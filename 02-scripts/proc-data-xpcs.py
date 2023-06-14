from Xana import Xana
import time
import sys
import os
import yaml
import numpy as np
import warnings

def run_analysis(proc_id, slurm_id, maskfile, setupfile, first, last, outdir, *args):
    warnings.filterwarnings("ignore")
    
    # watch out filelist.yml hard coded
    with open('filelist.yml', 'r') as f:
        folders = yaml.load(f, Loader=yaml.FullLoader)[proc_id]
    
    ana = Xana(fmtstr='ebs_id10_eiger500k',
                detector='eiger500k',
                maskfile=maskfile,
                setupfile=setupfile)

    outdir = f'/cfs/data/pg/sdaqs/esrf-ebs/id10/sc5275/20220614/processed/results/{outdir}/p{slurm_id:02d}'
    if not os.path.isdir(outdir):
        os.makedirs(outdir, exist_ok=True)
    ana.mksavdir(outdir)
    
    for folder in folders:
        ana.connect(folder)

    for index in ana.meta.index.values:
        ana.analyze(index, 'saxs',  first=first, last=last)
        if ana.meta.loc[index, 'nframes'] > 100:
            ana.analyze(index, 'xpcs', norm='symmetric', first=first, last=last, 
                        #twotime_par=[2,4,6], # twotime_par=list(np.arange(len(ana.setup.qroi))) 
                        # twotime_par=list(np.arange(1,len(ana.setup.qroi))),
                        twotime_par=[1, 2, 3, 4, 5, 6, 7, 8, 9, 11], # skip the roi where the detector is completely masked
			            saxs=None, nread_procs=1, nprocs=30, chunk_size=250, verbose=True)

if __name__ == '__main__':
    proc_id, slurm_id, maskfile, setupfile, first, last, outdir, *args = sys.argv[1:]
    run_analysis(int(proc_id), int(slurm_id), maskfile, setupfile, int(first), int(last), outdir)