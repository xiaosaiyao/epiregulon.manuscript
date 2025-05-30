#!/bin/bash
#SBATCH -n 20
#SBATCH --mem=64GB
#SBATCH --qos=long
#SBATCH --job-name="scenic plus main PBMC"
##SBATCH -p himem
cd scplus_pipeline/Snakemake
snakemake --cores 20
echo $(scontrol show jobid -dd $SLURM_JOB_ID)
