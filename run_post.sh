# K562_essential
# Rscript post_gaussian.R --cells K562_essential --response ratio --model varying_intercept_fixed_slope --efficient &
# Rscript post_gaussian.R --cells K562_essential --response ratio --model varying_intercept_varying_slope --efficient &
# Rscript post_gaussian.R --cells K562_essential --response ratio --model vivs_student --efficient &

# Rscript post_concordance.R --cells K562_essential --response ratio --model vifs_concordance --efficient &
# Rscript post_concordance.R --cells K562_essential --response ratio --model vivs_concordance --efficient &
# Rscript post_concordance.R --cells K562_essential --response ratio --model vivs_student_concordance --efficient &

Rscript post_logistic.R --cells K562_essential --model varying_intercept_fixed_slope --efficient &
Rscript post_logistic.R --cells K562_essential --model varying_intercept_varying_slope --efficient &
Rscript post_logistic.R --cells K562_essential --model vivs_student --efficient &

# K562_GenomeWide
# Rscript post_gaussian.R --cells K562_GenomeWide --response ratio --model varying_intercept_fixed_slope --efficient &
# Rscript post_gaussian.R --cells K562_GenomeWide --response ratio --model varying_intercept_varying_slope --efficient &
# Rscript post_gaussian.R --cells K562_GenomeWide --response ratio --model vivs_student --efficient &

# Rscript post_concordance.R --cells K562_GenomeWide --response ratio --model vifs_concordance --efficient &
# Rscript post_concordance.R --cells K562_GenomeWide --response ratio --model vivs_concordance --efficient &
# Rscript post_concordance.R --cells K562_GenomeWide --response ratio --model vivs_student_concordance --efficient &

Rscript post_logistic.R --cells K562_GenomeWide --model varying_intercept_fixed_slope --efficient &
Rscript post_logistic.R --cells K562_GenomeWide --model varying_intercept_varying_slope --efficient &
Rscript post_logistic.R --cells K562_GenomeWide --model vivs_student --efficient &

# Jurkat
# Rscript post_gaussian.R --cells Jurkat --response ratio --model varying_intercept_fixed_slope --efficient &
# Rscript post_gaussian.R --cells Jurkat --response ratio --model varying_intercept_varying_slope --efficient &
# Rscript post_gaussian.R --cells Jurkat --response ratio --model vivs_student --efficient &

# Rscript post_concordance.R --cells Jurkat --response ratio --model vifs_concordance --efficient &
# Rscript post_concordance.R --cells Jurkat --response ratio --model vivs_concordance --efficient &
# Rscript post_concordance.R --cells Jurkat --response ratio --model vivs_student_concordance --efficient &

# Rscript post_logistic.R --cells Jurkat --model varying_intercept_fixed_slope --efficient &
# Rscript post_logistic.R --cells Jurkat --model varying_intercept_varying_slope --efficient &
# Rscript post_logistic.R --cells Jurkat --model vivs_student --efficient &
