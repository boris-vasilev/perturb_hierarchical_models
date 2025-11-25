#!/bin/bash
. /etc/profile.d/modules.sh
module load gcc/11

#module load miniconda/3
source ~/.bashrc
conda activate env_nf

nextflow fit_models.nf  \
    --perturbList /rds/project/rds-csoP2nj6Y6Y/biv22/perturb_hierarchical_models/cross_screen_models/perturblist.txt