#!/bin/bash
#SBATCH -n 4
#SBATCH -N 1
#SBATCH --mem=200G
#SBATCH --qos=long

part_no=$1
ml R/dev
Rscript 01_rna_integration.R $part_no

