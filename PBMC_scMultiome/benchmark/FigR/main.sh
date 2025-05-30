#!/bin/bash
#SBATCH -n 20
#SBATCH --mem=256GB
#SBATCH --qos=long
#SBATCH --job-name="FigR_main_PBMC_250k"
#SBATCH -p himem
ml R/dev
Rscript FigR.R 20
echo $(scontrol show jobid -dd $SLURM_JOB_ID)

