#!/bin/bash
#SBATCH -A MRC-BSU-SL2-CPU
#SBATCH -p icelake-himem
#SBATCH -N 1  
#SBATCH --cpus-per-task=8 
#SBATCH -t 6:00:00

Rscript logistic.R "$@"
