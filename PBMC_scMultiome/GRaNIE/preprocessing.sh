#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=64GB
#SBATCH --qos=short
#SBATCH --job-name="GRaNIE preprocessing PBMC"
##SBATCH -p himem
ml R/dev
Rscript preprocessing.R 20
echo $(scontrol show jobid -dd $SLURM_JOB_ID)

