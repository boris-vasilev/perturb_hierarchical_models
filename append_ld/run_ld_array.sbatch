#!/usr/bin/env bash
#SBATCH -A MRC-BSU2-SL2-CPU
#SBATCH -p icelake-himem
#SBATCH --job-name=LD_array
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --array=1-192
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

. /etc/profile.d/modules.sh
module load openjdk/17.0.11_9/gcc/szfinnph
source ~/.bashrc
conda activate env_nf

# Get file corresponding to this array index
FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" gene_list.txt)

echo "Processing $FILE"

python compute_ld.py "$FILE"