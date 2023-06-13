# to be debug
# import os
# import sys
from glob import glob
import re
import argparse
import numpy as np
# import time
# import pickle
# import yaml

from tqdm import tqdm
# import pandas as pd
import h5py
# from pathlib import Path

from Xana import Xana

from functions.metadata import get_rep, get_scan_number
# from functions.slurm import submit_job

# run as for example
# python test_convert2hdf5.py --repsperspot 4 --proc ferritin_conc_gly_50_6 5
# or using sbatch job submission as
# sbatch convert_job.sbatch ferritin_conc_gly_50_6 5 4

parser = argparse.ArgumentParser()
parser.add_argument("datafolder", type=str, 
    help='the name of the measurement folder without scan number')
parser.add_argument("datasetnumber", type=int, help="the dataset number")
parser.add_argument('--repsperspot', '-reps', type=int, help="number of repetitions per spot", default=1)
parser.add_argument('--proc', action="store_true",  help="activate processing")

def convert2hdf5(xana, filename):
    xana.db = xana.db.sort_values(by='scannumber')
    xpcs_indices = xana.db[xana.db['analysis'] == 'xpcs'].index.values
    nxpcs = len(xpcs_indices)
    xpcs_scans = xana.db.loc[xpcs_indices, 'scannumber'].values
    saxs_indices = xana.db[xana.db['analysis'] == 'saxs'].index.values
    nsaxs = len(saxs_indices)
    saxs_scans = xana.db.loc[saxs_indices, 'scannumber'].values

    qI = xana.get_item(saxs_indices[0])['soq'][:,0]
    nqI = len(qI)
    tmp = xana.get_item(xpcs_indices[0])
    ttc_times = tmp['twotime_xy']
    ntimes = len(ttc_times)
    delay = tmp['corf'][1:,0]
    ndelay = len(delay)
    qv = tmp['qv']
    nq = len(qv)
    twotime_par = tmp['twotime_par']
    nttc = len(twotime_par)

    path = '/cfs/data/pg/sdaqs/esrf-ebs/id10/sc5275/20220614/processed/h5-files/'
    print("\n\npath+filename", path+filename, "\n\n")


    # -- write metadata from elog
    # load elog as pd dataframe
    elog = pd.read_pickle("../03-source/elog")
    condition = elog['measurement folder'].str.contains(filename[:-3], na=False) # select the right entry
    elog_entries = elog[condition] # select elog entry
    print(elog_entries)
    print("\nTemp: ", elog_entries["Temperature, K"].values, "\n")    
    # -- end of elog writing prep
    
    print((f"{nqI = }, {nq = }, {ntimes = }, {nxpcs = }, {nsaxs = }, {ndelay = }, {nttc = }"))

    with h5py.File( path + filename, 'a') as f: # w for read/write/create access
        # write metadata
        for col in elog.columns:
            f[f'/elog/{col}'] = elog_entries[col].values[0]
        
        # write data
        ttcs = f.create_dataset('/xpcs/ttcs/ttcs_all', shape=(nxpcs, nttc, ntimes, ntimes), dtype=np.float32, compression="gzip")
        ttcs_avg = f.create_dataset('/xpcs/ttcs/ttcs', shape=(nttc, ntimes, ntimes), dtype=np.float32, compression="gzip")
        ttcs_avg_f = f.create_dataset('/xpcs/ttcs/ttcs_f', shape=(nttc, ntimes, ntimes), dtype=np.float32, compression="gzip")
        g2s = f.create_dataset('/xpcs/g2s/g2s', shape=(nxpcs, nq, ndelay), dtype=np.float32, compression="gzip")
        I = f.create_dataset('/saxs/I', shape=(nsaxs, nqI), dtype=np.float32, compression="gzip")
        f['/xpcs/g2s/delay'] = delay
        f['/xpcs/g2s/q'] = qv
        f['/xpcs/ttcs/q'] = qv[twotime_par]
        f['/xpcs/ttcs/times'] = ttc_times
        f['/xpcs/scans'] = xpcs_scans
        f['/saxs/q'] = qI
        f['/saxs/scans'] = saxs_scans
        
        # updates
        print("defining ttcs and ttc_int")
        # ttcs = np.empty(shape=(nttc, ntimes, ntimes))
        ttc_int = np.empty(shape=nxpcs)

        print("\nentering for loop for writing xpcs results with ttcs")
        for i, xpcs_index in tqdm(enumerate(xpcs_indices[:10]), desc='writing XPCS results',  total=len(xpcs_indices)):
            tmp = xana.get_item(xpcs_index)
            # -- g2
            # g2s[i] = tmp['corf'][1:,1:].T
            # -- ttcs
            ttcs[i] = np.stack(list(tmp['twotime_corf'].values()), axis=0)
            # ttcs = list(tmp['twotime_corf'].values())
            # ttcs_avg[i] = = np.average(ttcs[i,:,:,:])
            # ttc_int[i] = np.average(ttcs)
            ttc_int[i] = np.average(list(tmp['twotime_corf'].values()))

        print("\nFor loop ended.\n")
        # filter based on ttc average intensity
        ttcs_avg_i = np.mean(ttc_int)
        print("Filtering ttcs based on average inteneisty which is: ", ttcs_avg_i)
        good_inds = ttc_int > ttcs_avg_i - 0.01*ttcs_avg_i # good_inds = ttc_int > 0.9 # or this ??
        print("good_inds: ", type(good_inds), len(good_inds), sum(good_inds), '\n')

        # issue here cayse ttcs is not defined anymore
        ttcs_avg = np.mean(ttcs,axis=0)
        ttcs_avg_f = np.mean(ttcs[good_inds],axis=0)
        print("written ttcs_avg and ttcs_avg_f ")

        for i, saxs_index in tqdm(enumerate(saxs_indices), desc='writing SAXS results',  total=len(saxs_indices)):
            tmp = xana.get_item(saxs_index)
            I[i] = tmp['soq'][:,1]
            # print(i, saxs_index)
        
if __name__ == "__main__":
    args = parser.parse_args()
    datafolder = args.datafolder
    datasetnumber = args.datasetnumber
    
    filename = f"{datafolder}_{datasetnumber:04}"
    print(filename)

    ana_db_files = glob(f'/cfs/data/pg/sdaqs/esrf-ebs/id10/sc5275/20220614/processed/results/{datafolder}_{datasetnumber:04}/p**/analysis_database.pkl', recursive=True)
    
    xana = Xana()
    for i, f in enumerate(ana_db_files):
        if i == 0:
            xana.load_db(f)
        else:
            xana.append_db(f)

    # find the error here
    # xana.db['rep'] = xana.db['master'].apply(get_rep, reps_per_spot=args.repsperspot)
    xana.db['rep'] = xana.db['master'].apply(lambda x: get_rep(x, reps_per_spot=args.repsperspot))
    xana.db['scannumber'] = xana.db['datdir'].apply(lambda x: get_scan_number(str(x)))
    # # filename = f"{datafolder}_{datasetnumber:04}"
    # xana.db = xana.db[xana.db['datdir'].apply(str).str.contains(filename)]
    # print(xana.db.head(10))

    if args.proc:
        convert2hdf5(xana, filename + '.h5')
