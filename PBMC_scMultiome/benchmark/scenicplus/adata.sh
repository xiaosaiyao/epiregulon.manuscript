#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=64GB
#SBATCH --qos=short
#SBATCH --job-name="adata PBMC"
ml Miniforge3
conda activate scenicplus2
python adata.py
echo $(scontrol show jobid -dd $SLURM_JOB_ID)
echo $(sacct --units=G --format=jobid,jobname%50,maxvmsizenode,maxvmsize,avevmsize,avecpu,CPUTime,consumedenergy,MaxRSS,alloccpus,elapsed,exitcode -j $SLURM_JOB_ID -P)

