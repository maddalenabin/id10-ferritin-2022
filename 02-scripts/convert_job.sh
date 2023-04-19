#!/bin/bash
#SBATCH -p cops
#SBATCH --output=/cfs/home/mabi3848/id10-ferritin-2022/04-jobs/convert-slurm-%j.out
#SBATCH --error=/cfs/home/mabi3848/id10-ferritin-2022/04-jobs/convert-slurm-%j.err
#SBATCH --time 24:00:00
#SBATCH --job-name=id10-convert
#SBATCH --mail-type=ALL
#SBATCH --mail-user=maddalenabin@hotmail.com

source /cfs/home/mabi3848/id10-ferritin-2022/.venv/bin/activate
echo "SLURM_JOB_ID           $SLURM_JOB_ID"
echo "SLURM_ARRAY_TASK_ID    $SLURM_ARRAY_TASK_ID"

python convert2hdf5.py --repsperspot $3 --proc $1 $2

exit