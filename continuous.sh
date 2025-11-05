#!/bin/bash
#SBATCH -A MRC-BSU-SL2-CPU
#SBATCH -p icelake
#SBATCH -N 1  
#SBATCH --cpus-per-task=4 
#SBATCH -t 3:00:00

#Previous time was 3:00:00
Rscript continuous.R "$@"  
