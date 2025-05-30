#!/bin/bash
#SBATCH -n 20
#SBATCH --mem=128GB
#SBATCH --qos=long
#SBATCH --job-name="GRaNIE main PBMC"
##SBATCH -p himem
ml R/dev
Rscript main.R 20
echo $(scontrol show jobid -dd $SLURM_JOB_ID)
