#!/bin/bash
#SBATCH -A MRC-BSU-SL2-CPU
#SBATCH -p icelake
#SBATCH -N 1  
#SBATCH --cpus-per-task=24 
#SBATCH -t 12:00:00

Rscript logistic.R "$@"
