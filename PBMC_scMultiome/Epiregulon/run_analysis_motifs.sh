#!/bin/bash
#SBATCH -n 10
#SBATCH --mem=128GB
#SBATCH --qos=medium
#SBATCH --job-name="PBMC Epiregulon motifs"
##SBATCH -p himem
ml R/dev
Rscript motifs.R
