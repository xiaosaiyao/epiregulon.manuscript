#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=64GB
#SBATCH --qos=short
#SBATCH --job-name="peak_matrix PBMC"
ml add R/dev
Rscript Scenic_plus_preparation.R
echo $(scontrol show jobid -dd $SLURM_JOB_ID)
bash /gstore/project/epigen/benchmark/PBMC/save_stats.sh $SLURM_JOB_ID
