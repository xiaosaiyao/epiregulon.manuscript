#!/bin/bash
#SBATCH -n 20
#SBATCH --mem=256GB
#SBATCH --qos=long
#SBATCH --job-name="topics PBMC"
#SBATCH -p himem
ml Miniforge3
conda activate scenicplus2
python cistopic.py 20
python cistopic_2.py 20
echo $(scontrol show jobid -dd $SLURM_JOB_ID)
bash /gstore/project/epigen/benchmark/PBMC/save_stats.sh $SLURM_JOB_ID
