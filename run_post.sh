# K562_essential
# Rscript post_gaussian.R --cells K562_essential --response ratio --bulk --model varying_intercept_fixed_slope --efficient &
# Rscript post_gaussian.R --cells K562_essential --response ratio --bulk --model varying_intercept_varying_slope --efficient &
# Rscript post_gaussian.R --cells K562_essential --response ratio --bulk --model vivs_student --efficient &

# Rscript post_concordance.R --cells K562_essential --response ratio --bulk --model vifs_concordance --efficient &
# Rscript post_concordance.R --cells K562_essential --response ratio --bulk --model vivs_concordance --efficient &
# Rscript post_concordance.R --cells K562_essential --response ratio --bulk --model vivs_student_concordance --efficient &

Rscript post_logistic.R --cells K562_essential --bulk --model varying_intercept_fixed_slope --efficient &
Rscript post_logistic.R --cells K562_essential --bulk --model varying_intercept_varying_slope --efficient &
Rscript post_logistic.R --cells K562_essential --bulk --model vivs_student --efficient &

# K562_GenomeWide
# Rscript post_gaussian.R --cells K562_GenomeWide --response ratio --bulk --model varying_intercept_fixed_slope --efficient &
# Rscript post_gaussian.R --cells K562_GenomeWide --response ratio --bulk --model varying_intercept_varying_slope --efficient &
# Rscript post_gaussian.R --cells K562_GenomeWide --response ratio --bulk --model vivs_student --efficient &

# Rscript post_concordance.R --cells K562_GenomeWide --response ratio --bulk --model vifs_concordance --efficient &
# Rscript post_concordance.R --cells K562_GenomeWide --response ratio --bulk --model vivs_concordance --efficient &
# Rscript post_concordance.R --cells K562_GenomeWide --response ratio --bulk --model vivs_student_concordance --efficient &

Rscript post_logistic.R --cells K562_GenomeWide --bulk --model varying_intercept_fixed_slope --efficient &
Rscript post_logistic.R --cells K562_GenomeWide --bulk --model varying_intercept_varying_slope --efficient &
Rscript post_logistic.R --cells K562_GenomeWide --bulk --model vivs_student --efficient &

# Jurkat
# Rscript post_gaussian.R --cells Jurkat --response ratio --bulk --model varying_intercept_fixed_slope --efficient &
# Rscript post_gaussian.R --cells Jurkat --response ratio --bulk --model varying_intercept_varying_slope --efficient &
# Rscript post_gaussian.R --cells Jurkat --response ratio --bulk --model vivs_student --efficient &

# Rscript post_concordance.R --cells Jurkat --response ratio --bulk --model vifs_concordance --efficient &
# Rscript post_concordance.R --cells Jurkat --response ratio --bulk --model vivs_concordance --efficient &
# Rscript post_concordance.R --cells Jurkat --response ratio --bulk --model vivs_student_concordance --efficient &

# Rscript post_logistic.R --cells Jurkat --bulk --model varying_intercept_fixed_slope --efficient &
# Rscript post_logistic.R --cells Jurkat --bulk --model varying_intercept_varying_slope --efficient &
# Rscript post_logistic.R --cells Jurkat --bulk --model vivs_student --efficient &
