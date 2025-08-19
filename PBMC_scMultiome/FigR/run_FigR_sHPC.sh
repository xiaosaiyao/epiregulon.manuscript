#!/bin/bash
#BSUB -n 12
#BSUB -R "rusage[mem=12GB]"
#BSUB -o my_job.o%J
#BSUB -e my_job.e%J
#BSUB -q long
ml add CEDAR
ml add R/cedar_r4.4_bioc3.19-release
Rscript FigR.R 12
