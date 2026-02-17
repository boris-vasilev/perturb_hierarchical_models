#!/usr/bin/env bash
#SBATCH -A MRC-BSU2-SL2-CPU
#SBATCH -p icelake-himem
#SBATCH -N 1
#SBATCH --job-name=dat_genome_wide
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=30

. /etc/profile.d/modules.sh

module load openjdk/17.0.11_9/gcc/szfinnph
source ~/.bashrc
conda activate env_nf

python -u dat_genome_wide.py