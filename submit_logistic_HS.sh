# Jurkat
#  sbatch --job-name="l-J-vivs" logistic_HS.sh --cells Jurkat --model varying_intercept_varying_slope --efficient

# K562_essential
#  sbatch --job-name="l-K-vivs" logistic_HS.sh --cells K562_essential --model varying_intercept_varying_slope --efficient

# K562_GenomeWide
sbatch --job-name="l-KG-vivs" logistic_HS.sh --cells K562_GenomeWide --model varying_intercept_varying_slope

# HepG2
# sbatch --job-name="l-H-vivs" logistic_HS.sh --cells HepG2 --model varying_intercept_varying_slope --efficient

# RPE1
# sbatch --job-name="l-R-vivs" logistic_HS.sh --cells RPE1 --model varying_intercept_varying_slope --efficient
