# Jurkat
sbatch --job-name="eff-l-J-vifs" logistic.sh  --cells Jurkat --model varying_intercept_fixed_slope --efficient --bulk
#sbatch --job-name="eff-l-J-student" logistic.sh --cells Jurkat --model vivs_student --efficient --bulk
sbatch --job-name="eff-l-J-vivs" logistic.sh --cells Jurkat --model varying_intercept_varying_slope --efficient --bulk

# K562_essential
sbatch --job-name="eff-l-K-vifs" logistic.sh --cells K562_essential --model varying_intercept_fixed_slope --efficient --bulk
#sbatch --job-name="eff-l-K-student" logistic.sh --cells K562_essential --model vivs_student --efficient --bulk
sbatch --job-name="eff-l-K-vivs" logistic.sh --cells K562_essential --model varying_intercept_varying_slope --efficient --bulk
