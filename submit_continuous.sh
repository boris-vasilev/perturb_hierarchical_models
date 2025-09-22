# Jurkat
#sbatch --job-name="eff-c-J-vifs-ratio" continuous.sh --cells Jurkat --model varying_intercept_fixed_slope --response ratio --efficient --bulk
#sbatch --job-name="eff-c-J-vivs-student-ratio" continuous.sh --cells Jurkat --model vivs_student --response ratio --efficient --bulk
#sbatch continuous.sh --cells Jurkat --model no_pooling_intercept_varying_slope
#sbatch --job-name="eff-c-J-vivs-ratio" continuous.sh --cells Jurkat --model varying_intercept_varying_slope --response ratio --efficient --bulk
sbatch --job-name="eff-concord-J-vivs-ratio" continuous.sh --cells Jurkat --model vivs_concordance --response ratio --efficient --bulk
sbatch --job-name="eff-concord-J-vifs-ratio" continuous.sh --cells Jurkat --model vifs_concordance --response ratio --efficient --bulk
#sbatch --job-name="eff-concord-J-vivs-student-ratio" continuous.sh --cells Jurkat --model vivs_student_concordance --response ratio --efficient --bulk



# K562_essential
#sbatch --job-name="eff-c-K-vifs-ratio" continuous.sh --cells K562_essential --model varying_intercept_fixed_slope --response ratio --efficient --bulk
#sbatch --job-name="eff-c-K-vivs-student-ratio" continuous.sh --cells K562_essential --model vivs_student --response ratio --efficient --bulk
#sbatch continuous.sh --cells K562_essential --model no_pooling_intercept_varying_slope
#sbatch --job-name="eff-c-K-vivs-ratio" continuous.sh --cells K562_essential --model varying_intercept_varying_slope --response ratio --efficient --bulk
sbatch --job-name="eff-concord-K-vivs-ratio" continuous.sh --cells K562_essential --model vivs_concordance --response ratio --efficient --bulk
sbatch --job-name="eff-concord-K-vifs-ratio" continuous.sh --cells K562_essential --model vifs_concordance --response ratio --efficient --bulk
#sbatch --job-name="eff-concord-K-vivs-student-ratio" continuous.sh --cells K562_essential --model vivs_student_concordance --response ratio --efficient --bulk

