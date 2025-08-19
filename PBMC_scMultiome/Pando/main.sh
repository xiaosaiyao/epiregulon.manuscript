#!/bin/bash
#SBATCH -n 20
#SBATCH --mem=64GB
#SBATCH --qos=long
#SBATCH --job-name="Pando main PBMC SCT"
##SBATCH -p himem
ml R/r430-bioc317-20230522_prd
Rscript main.R
echo $(scontrol show jobid -dd $SLURM_JOB_ID)
