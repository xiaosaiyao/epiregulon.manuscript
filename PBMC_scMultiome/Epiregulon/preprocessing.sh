#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=64GB
#SBATCH --qos=short
#SBATCH --job-name="Epiregulon preprocessing PBMC"
ml R/dev
Rscript preprocessing.R
echo $(scontrol show jobid -dd $SLURM_JOB_ID)
