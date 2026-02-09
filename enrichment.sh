#!/bin/bash
#SBATCH -A MRC-BSU2-SL2-CPU
#SBATCH -p icelake
#SBATCH -N 1  
#SBATCH --cpus-per-task=24 
#SBATCH -t 24:00:00

Rscript enrichment.R "$@"
