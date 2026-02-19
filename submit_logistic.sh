#!/bin/bash

# Jurkat
sbatch --job-name="l-J-vivs" logistic.sh --screen Jurkat --model vivs
sbatch --job-name="l-J-vivs_hs" logistic.sh --screen Jurkat --model vivs_horseshoe

# K562_essential
sbatch --job-name="l-KE-vivs" logistic.sh --screen K562_essential --model vivs
sbatch --job-name="l-KE-vivs_hs" logistic.sh --screen K562_essential --model vivs_horseshoe

# K562_GenomeWide
sbatch --job-name="l-KG-vivs" logistic.sh --screen K562_GenomeWide --model vivs
sbatch --job-name="l-KG-vivs_hs" logistic.sh --screen K562_GenomeWide --model vivs_horseshoe

# HepG2
sbatch --job-name="l-H-vivs" logistic.sh --screen HepG2 --model vivs
sbatch --job-name="l-H-vivs_hs" logistic.sh --screen HepG2 --model vivs_horseshoe

# RPE1
sbatch --job-name="l-R-vivs" logistic.sh --screen RPE1 --model vivs
sbatch --job-name="l-R-vivs_hs" logistic.sh --screen RPE1 --model vivs_horseshoe