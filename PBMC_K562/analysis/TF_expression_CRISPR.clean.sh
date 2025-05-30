#!/bin/sh
#SBATCH -n 4
#SBATCH -N 1
#SBATCH --mem=200G
#SBATCH --qos=long


ml R/dev
Rscript TF_expression_CRISPR.clean.R


