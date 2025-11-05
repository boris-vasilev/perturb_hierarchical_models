# Jurkat
#  sbatch --job-name="l-J-vifs" logistic.sh  --cells Jurkat --model varying_intercept_fixed_slope --efficient
# sbatch --job-name="l-J-student" logistic.sh --cells Jurkat --model vivs_student --efficient
 sbatch --job-name="l-J-vivs" logistic.sh --cells Jurkat --model varying_intercept_varying_slope --efficient

# K562_essential
#  sbatch --job-name="l-K-vifs" logistic.sh --cells K562_essential --model varying_intercept_fixed_slope --efficient
# sbatch --job-name="l-K-student" logistic.sh --cells K562_essential --model vivs_student --efficient
 sbatch --job-name="l-K-vivs" logistic.sh --cells K562_essential --model varying_intercept_varying_slope --efficient

# K562_GenomeWide
# sbatch --job-name="l-KG-vifs" logistic.sh --cells K562_GenomeWide --model varying_intercept_fixed_slope --efficient
# sbatch --job-name="l-KG-student" logistic.sh --cells K562_GenomeWide --model vivs_student --efficient
# sbatch --job-name="l-KG-vivs" logistic.sh --cells K562_GenomeWide --model varying_intercept_varying_slope --efficient

# HepG2
# sbatch --job-name="l-H-vifs" logistic.sh --cells HepG2 --model varying_intercept_fixed_slope --efficient
# sbatch --job-name="l-H-student" logistic.sh --cells HepG2 --model vivs_student --efficient
sbatch --job-name="l-H-vivs" logistic.sh --cells HepG2 --model varying_intercept_varying_slope --efficient

# RPE1
# sbatch --job-name="l-R-vifs" logistic.sh --cells RPE1 --model varying_intercept_fixed_slope --efficient
# sbatch --job-name="l-R-student" logistic.sh --cells RPE1 --model vivs_student --efficient
sbatch --job-name="l-R-vivs" logistic.sh --cells RPE1 --model varying_intercept_varying_slope --efficient
