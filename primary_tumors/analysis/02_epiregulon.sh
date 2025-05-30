#!/bin/sh
#SBATCH -n 4
#SBATCH -N 1
#SBATCH --mem=200G
#SBATCH --qos=long

part_no=$1
ml R/dev
Rscript 02_epiregulon."$part_no".R $part_no


