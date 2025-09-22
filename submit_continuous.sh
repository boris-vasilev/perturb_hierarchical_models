# Jurkat
# sbatch --job-name="c-J-vifs" continuous.sh --cells Jurkat --model varying_intercept_fixed_slope --response ratio --efficient
# sbatch --job-name="c-J-vivs-student" continuous.sh --cells Jurkat --model vivs_student --response ratio --efficient
# sbatch continuous.sh --cells Jurkat --model no_pooling_intercept_varying_slope
# sbatch --job-name="c-J-vivs" continuous.sh --cells Jurkat --model varying_intercept_varying_slope --response ratio --efficient
# sbatch --job-name="concord-J-vivs" continuous.sh --cells Jurkat --model vivs_concordance --response ratio --efficient
# sbatch --job-name="concord-J-vifs" continuous.sh --cells Jurkat --model vifs_concordance --response ratio --efficient
# sbatch --job-name="concord-J-vivs-student" continuous.sh --cells Jurkat --model vivs_student_concordance --response ratio --efficient



# K562_essential
# sbatch --job-name="c-K-vifs" continuous.sh --cells K562_essential --model varying_intercept_fixed_slope --response ratio --efficient
# sbatch --job-name="c-K-vivs-student" continuous.sh --cells K562_essential --model vivs_student --response ratio --efficient
# sbatch continuous.sh --cells K562_essential --model no_pooling_intercept_varying_slope
# sbatch --job-name="c-K-vivs" continuous.sh --cells K562_essential --model varying_intercept_varying_slope --response ratio --efficient
# sbatch --job-name="concord-K-vivs" continuous.sh --cells K562_essential --model vivs_concordance --response ratio --efficient
# sbatch --job-name="concord-K-vifs" continuous.sh --cells K562_essential --model vifs_concordance --response ratio --efficient
# sbatch --job-name="concord-K-vivs-student" continuous.sh --cells K562_essential --model vivs_student_concordance --response ratio --efficient

# K562_GenomeWide
# sbatch --job-name="c-KG-vifs" continuous.sh --cells K562_GenomeWide --model varying_intercept_fixed_slope --response ratio --efficient
# sbatch --job-name="c-KG-vivs-student" continuous.sh --cells K562_GenomeWide --model vivs_student --response ratio --efficient
# sbatch --job-name="c-KG-vivs" continuous.sh --cells K562_GenomeWide --model varying_intercept_varying_slope --response ratio --efficient
sbatch --job-name="concord-KG-vivs" continuous.sh --cells K562_GenomeWide --model vivs_concordance --response ratio --efficient
sbatch --job-name="concord-KG-vifs" continuous.sh --cells K562_GenomeWide --model vifs_concordance --response ratio --efficient
# sbatch --job-name="concord-KG-vivs-student" continuous.sh --cells K562_GenomeWide --model vivs_student_concordance --response ratio --efficient

# HepG2
# sbatch --job-name="c-H-vifs" continuous.sh --cells HepG2 --model varying_intercept_fixed_slope --response ratio --efficient
# sbatch --job-name="c-H-vivs-student" continuous.sh --cells HepG2 --model vivs_student --response ratio --efficient
# sbatch --job-name="c-H-vivs" continuous.sh --cells HepG2 --model varying_intercept_varying_slope --response ratio --efficient
sbatch --job-name="concord-H-vivs" continuous.sh --cells HepG2 --model vivs_concordance --response ratio --efficient
sbatch --job-name="concord-H-vifs" continuous.sh --cells HepG2 --model vifs_concordance --response ratio --efficient
# sbatch --job-name="concord-H-vivs-student" continuous.sh --cells HepG2 --model vivs_student_concordance --response ratio --efficient

# RPE1
# sbatch --job-name="c-R-vifs" continuous.sh --cells RPE1 --model varying_intercept_fixed_slope --response ratio --efficient
# sbatch --job-name="c-R-vivs-student" continuous.sh --cells HepG2 --model vivs_student --response ratio --efficient
# sbatch --job-name="c-R-vivs" continuous.sh --cells RPE1 --model varying_intercept_varying_slope --response ratio --efficient
sbatch --job-name="concord-R-vivs" continuous.sh --cells RPE1 --model vivs_concordance --response ratio --efficient
sbatch --job-name="concord-R-vifs" continuous.sh --cells RPE1 --model vifs_concordance --response ratio --efficient
# sbatch --job-name="concord-R-vivs-student" continuous.sh --cells RPE1 --model vivs_student_concordance --response ratio --efficient
