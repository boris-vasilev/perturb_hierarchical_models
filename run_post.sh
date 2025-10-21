# K562_essential

# Rscript post_concordance.R --cells K562_essential --response ratio --model vifs_concordance --efficient &
# Rscript post_concordance.R --cells K562_essential --response ratio --model vivs_concordance --efficient &

Rscript post_logistic.R --cells K562_essential --model varying_intercept_fixed_slope --efficient &
Rscript post_logistic.R --cells K562_essential --model varying_intercept_varying_slope --efficient &

# K562_GenomeWide

# Rscript post_concordance.R --cells K562_GenomeWide --response ratio --model vifs_concordance --efficient &
# Rscript post_concordance.R --cells K562_GenomeWide --response ratio --model vivs_concordance --efficient &

# Rscript post_logistic.R --cells K562_GenomeWide --model varying_intercept_fixed_slope --efficient &
# Rscript post_logistic.R --cells K562_GenomeWide --model varying_intercept_varying_slope --efficient &

# Jurkat

# Rscript post_concordance.R --cells Jurkat --response ratio --model vifs_concordance --efficient &
# Rscript post_concordance.R --cells Jurkat --response ratio --model vivs_concordance --efficient &

Rscript post_logistic.R --cells Jurkat --model varying_intercept_fixed_slope --efficient &
Rscript post_logistic.R --cells Jurkat --model varying_intercept_varying_slope --efficient &

# RPE1

# Rscript post_concordance.R --cells RPE1 --response ratio --model vifs_concordance --efficient &
# Rscript post_concordance.R --cells RPE1 --response ratio --model vivs_concordance --efficient &

Rscript post_logistic.R --cells RPE1 --model varying_intercept_fixed_slope --efficient &
Rscript post_logistic.R --cells RPE1 --model varying_intercept_varying_slope --efficient &

# HepG2

# Rscript post_concordance.R --cells HepG2 --response ratio --model vifs_concordance --efficient &
# Rscript post_concordance.R --cells HepG2 --response ratio --model vivs_concordance --efficient &

Rscript post_logistic.R --cells HepG2 --model varying_intercept_fixed_slope --efficient &
Rscript post_logistic.R --cells HepG2 --model varying_intercept_varying_slope --efficient &
