#!/usr/bin/env python
# coding: utf-8

import os
import sys
import argparse
import numpy as np
import time
import pickle
import yaml
from glob import glob
import re
from pathlib import Path

from functions.helpers import dump_filelist
from functions.metadata import get_scan_number
from functions.slurm import submit_job

parser = argparse.ArgumentParser()
parser.add_argument("datafolder", type=str, help='the name of the measurement folder without scan number')
parser.add_argument("datasetnumber", type=int, help="the dataset number")
parser.add_argument('--nprocs', type=int, help="number of subprocesses", default=4)
parser.add_argument('--proc', action="store_true",  help="activate processing")
# parser.add_argument('--filelist-name', default='filelist.yml', help="name of the file to store the filelist.")

# examples:
# datafolder = 'ferritin_conc120_gly_50_2'
# datasetnumber = 2
# run for example with 
# python analysis-script.py ferritin_conc120_gly_50_2 2 --proc

def make_filelist(datafolder, datasetnumber):
    datafolder_ = f"/cfs/data/pg/sdaqs/esrf-ebs/id10/sc5275/{datafolder}/{datafolder}_{datasetnumber:04d}"
    filelist = glob(datafolder_ + '/**/*.h5',  recursive=True)
    filelist = sorted(list(filter(lambda x: 'scan' in x, filelist)))
    filelist = list(map(lambda x: str(Path(x).parent), filelist))

    if not len(filelist):
        raise FileNotFoundError(f"Folder '{datafolder_}' does not contain scans.")
    
    print(f"{len(filelist) = :d}")
    print(f"{filelist[0] = :s}")
    print(f"{filelist[-1] = :s}")

    return filelist

if __name__ == "__main__":
    args = parser.parse_args()
    datafolder = args.datafolder
    datasetnumber = args.datasetnumber

    filelist = make_filelist(datafolder, datasetnumber)
    scans = np.array(list(map(get_scan_number, filelist)))

    # dump_filelist(filelist, args.filelist_name, nprocs=args.nprocs)
    dump_filelist(filelist, f"../05-filelists/{datafolder}_{datasetnumber:04d}.yml", nprocs=args.nprocs)

    first = 10
    last = 10000
    func_file = 'proc-data-xpcs.py'
    test = not args.proc
    args_list = ['/cfs/data/pg/sdaqs/esrf-ebs/id10/sc5275/20220614/processed/mask-setup/cryo-mask-230417_03.npy',
            '/cfs/data/pg/sdaqs/esrf-ebs/id10/sc5275/20220614/processed/mask-setup/setup-fullmask-cryo-230417.pkl',
            first,
            last,
            f'{datafolder}_{datasetnumber:04d}',
           ]
    submit_job(func_file, args.nprocs, args=args_list, test=test)