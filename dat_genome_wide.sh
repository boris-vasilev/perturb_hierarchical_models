#!/bin/bash

sbatch --job-name="dat_K562GW" dat_genome_wide.sbatch --screen K562_GenomeWide --output dat_K562GW.csv
sbatch --job-name="dat_K562E" dat_genome_wide.sbatch --screen K562_Essential --output dat_K562E.csv
sbatch --job-name="dat_Jurkat" dat_genome_wide.sbatch --screen Jurkat --output dat_Jurkat.csv
sbatch --job-name="dat_HepG2" dat_genome_wide.sbatch --screen HepG2 --output dat_HepG2.csv
sbatch --job-name="dat_RPE1" dat_genome_wide.sbatch --screen RPE1 --output dat_RPE1.csv