#!/bin/bash
#SBATCH -n 20
#SBATCH --mem=64GB
#SBATCH --qos=medium
#SBATCH --job-name="Epiregulon main PBMC corr"
##SBATCH -p himem
ml R/dev
Rscript main.R
echo $(scontrol show jobid -dd $SLURM_JOB_ID)
