#!/bin/bash
. /etc/profile.d/modules.sh
module load gcc/11

source ~/.bashrc

conda activate env_nf

python QC.py --cells Jurkat
python QC.py --cells HepG2
python QC.py --cells K562_GenomeWide
python QC.py --cells K562_essential
python QC.py --cells RPE1
