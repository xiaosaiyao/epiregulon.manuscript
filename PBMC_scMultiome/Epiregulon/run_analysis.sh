#!/bin/bash
#SBATCH -n 10
#SBATCH --mem=128GB
#SBATCH --qos=medium
#SBATCH --job-name="PBMC Epiregulon blood"
##SBATCH -p himem
ml R/dev
Rscript generate_regulon_df_blood_grl.R
