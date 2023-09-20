import os
from glob import glob
import re
import argparse
import numpy as np
from tqdm import tqdm
import pandas as pd
import h5py
from Xana import Xana
from functions.metadata import get_rep, get_scan_number

# run as for example
# python convert2hdf5.py --repsperspot 4 --proc ferritin_conc_gly_50_6 5
# or using sbatch job submission as
# sbatch job_convert.sbatch ferritin_conc_gly_50_5 4 4

parser = argparse.ArgumentParser()
parser.add_argument("datafolder", type=str, help='the name of the measurement folder without scan number')
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
    print("\nWriting metadata from elog\n")
    elog = pd.read_pickle("../03-source/elog")
    condition = elog['measurement folder'].str.contains(filename[:-3], na=False) # select the right entry
    elog_entries = elog[condition] # select elog entry
    print(elog_entries)
    print("\nTemp: ", elog_entries["Temperature, K"].values, "\n")    
    # -- end of elog writing prep
    

    print(f"Create h5 file: {path + filename}", '\n')
    nreps = args.repsperspot
    print((f"{nqI = }, {nq = }, {ntimes = }, {nxpcs = }, {nsaxs = }, {ndelay = }, {nttc = }, {nreps = }"))

    with h5py.File( path + filename, 'a') as f: # w for read/write/create access
        # write metadata
        for col in elog.columns:
            f[f'/elog/{col}'] = str(elog_entries[col].values[0])
        
        # write data
        f['/xpcs/g2s/delay'] = delay
        f['/xpcs/g2s/q'] = qv
        f['/xpcs/ttcs/q'] = qv[twotime_par]
        f['/xpcs/ttcs/times'] = ttc_times
        f['/xpcs/ttcs/twotime_par'] = twotime_par
        f['/xpcs/scans'] = xpcs_scans
        f['/saxs/q'] = qI
        f['/saxs/scans'] = saxs_scans
        
        good_indices = []

        # -- divide per repetition
        reps = np.arange(1, args.repsperspot+1)
        ttcs_filtered = f.create_dataset('/xpcs/ttcs/ttc_rep_qs_avg_filtered', shape=(args.repsperspot, nttc, ntimes, ntimes), dtype=np.float32, compression="gzip")
        g2s_filtered = f.create_dataset('/xpcs/g2s/g2s_filtered', shape=(args.repsperspot, nq, ndelay), dtype=np.float32, compression="gzip") # filtered g2s
        dg2s = f.create_dataset('/xpcs/g2s/g2s_err', shape=(args.repsperspot, nq, ndelay), dtype=np.float32, compression="gzip") # filtered g2s
        # not filtered datasets, which are not used anymore
        # ttcs = f.create_dataset('/xpcs/ttcs/ttc_rep_qs_avg', shape=(args.repsperspot, nttc, ntimes, ntimes), dtype=np.float32, compression="gzip")
        # g2s = f.create_dataset('/xpcs/g2s/g2s', shape=(args.repsperspot, nq, ndelay), dtype=np.float32, compression="gzip") # all g2s

        print("Number of repetitions: ", args.repsperspot, '\n')

        print("Writing xpcs data")

        for j,rep in enumerate(reps):
            # temporarily here, in the rep loop it will be overwritten
            ind_xpcs = xana.db[(xana.db['analysis'] == 'xpcs') & (xana.db['rep'] == rep)].index.values
            
            print("\trepetition: ", rep, '\t number of spots: ', len(ind_xpcs))
            ttcs_qs = np.empty(shape=(len(ind_xpcs), nttc, ntimes, ntimes))
            g2s_qs = np.empty(shape=(len(ind_xpcs), nq, ndelay))

            baseline = np.empty(shape=len(ind_xpcs)) # test

            for i,ind in tqdm(enumerate(ind_xpcs), total=len(ind_xpcs)):
                ttcs_qs[i] = np.stack(list(xana.get_item(ind)['twotime_corf'].values()), axis=0) # (nspots, nqs, t1, t2)
                g2s_qs[i] = xana.get_item(ind)['corf'][1:,1:].T # (nspots, nqs, delays)
                baseline[i] = np.average(ttcs_qs[i][1][:300,-300:]) # test

            # not filtered datasets, not saved anymore
            # g2s[j] = np.average(g2s_qs, axis=0)
            # ttcs[j] = np.average(ttcs_qs, axis=0)

            threshold = 1.03
            cond = (baseline < threshold)
            good_indices.append(cond) # good_runs
            print("\t\t\t good runs: ", np.sum(good_indices[j]))
            ttcs_filtered[j] = np.average(ttcs_qs[good_indices[j],:,:,:], axis=0)
            g2s_filtered[j] = np.average(g2s_qs[good_indices[j],:,:], axis=0)
            dg2s[j] = np.std(g2s_qs[good_indices[j],:,:], axis=0) / np.sqrt(np.sum(good_indices[j])-1)
            
            f.create_dataset(f'/xpcs/good_indices_{j}', data=good_indices[j]) # good indices used for filtering
            print("\t\t\t good indices saved: ", np.shape(good_indices[j]))

        # f.create_dataset('/xpcs/good_indices', data=np.asarray(good_indices)) # good indices used for filtering
        # print("\nSaved good indices with shape ", np.shape(good_indices))

        # -- saxs
        print("\nWriting saxs data")
        I = f.create_dataset('/saxs/I', shape=(nsaxs, nqI), dtype=np.float32, compression="gzip") # all Iqs
        # I_rep = f.create_dataset('/saxs/I_reps', shape=(args.repsperspot, nqI), dtype=np.float32, compression="gzip") # Iqs per repetition 
        I_rep_filtered = f.create_dataset('/saxs/I_reps_filtered', shape=(args.repsperspot, nqI), dtype=np.float32, compression="gzip") # filtered Iqs
        
        # save Iq per spot
        print("\tall spots")
        for i, saxs_index in tqdm(enumerate(saxs_indices), desc='writing SAXS results',  total=len(saxs_indices)):
            tmp = xana.get_item(saxs_index)
            I[i] = tmp['soq'][:,1]
        
        # save avg Iq per repetition
        for j,rep in enumerate(reps):
            print("\trepetition: ", rep, '\t number of spots: ', len(ind_xpcs))
            ind_saxs = xana.db[(xana.db['analysis'] == 'saxs') & (xana.db['rep'] == rep)].index.values
            Is = np.empty(shape=(len(ind_saxs), nqI))

            for i,ind in tqdm(enumerate(ind_saxs), total=len(ind_saxs)):
                Is[i] = xana.get_item(ind)['soq'][:,1]
                
            # I_rep[j] = np.average(Is, axis=0)
            I_rep_filtered[j] = np.average(Is[good_indices[j],:], axis=0)

        print("\n\nEnd of the scripts\n")


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
    xana.db['rep'] = xana.db['master'].apply(lambda x: get_rep(x, reps_per_spot=args.repsperspot))
    xana.db['scannumber'] = xana.db['datdir'].apply(lambda x: get_scan_number(str(x)))

    if args.proc:
        convert2hdf5(xana, filename + '.h5')
