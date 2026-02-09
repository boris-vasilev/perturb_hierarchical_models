# Jurkat
sbatch --job-name="l-J-vivs" logistic.sh --cells Jurkat --model vivs
sbatch --job-name="l-J-vivs_hs" logistic.sh --cells Jurkat --model vivs_horseshoe

# K562_essential
sbatch --job-name="l-KE-vivs" logistic.sh --cells K562_essential --model vivs
sbatch --job-name="l-KE-vivs_hs" logistic.sh --cells K562_essential --model vivs_horseshoe

# K562_GenomeWide
# sbatch --job-name="l-KG-vivs" logistic.sh --cells K562_GenomeWide --model vivs
# sbatch --job-name="l-KG-vivs_hs" logistic.sh --cells K562_GenomeWide --model vivs_horseshoe

# HepG2
sbatch --job-name="l-H-vivs" logistic.sh --cells HepG2 --model vivs
sbatch --job-name="l-H-vivs_hs" logistic.sh --cells HepG2 --model vivs_horseshoe

# RPE1
sbatch --job-name="l-R-vivs" logistic.sh --cells RPE1 --model vivs
sbatch --job-name="l-R-vivs_hs" logistic.sh --cells RPE1 --model vivs_horseshoe
