#!/bin/bash

# Jurkat
sbatch --job-name="e-J-fixed" enrichment.sh --cells Jurkat --model binom_enrichment_fixed --chains 6 --threads 4
sbatch --job-name="e-J-varying" enrichment.sh --cells Jurkat --model binom_enrichment_varying_intercept --chains 6 --threads 4
sbatch --job-name="e-J-cov" enrichment.sh --cells Jurkat --model binom_enrichment_with_covariates --chains 6 --threads 4

# K562_essential
sbatch --job-name="e-KE-fixed" enrichment.sh --cells K562_essential --model binom_enrichment_fixed --chains 6 --threads 4
sbatch --job-name="e-KE-varying" enrichment.sh --cells K562_essential --model binom_enrichment_varying_intercept --chains 6 --threads 4
sbatch --job-name="e-KE-cov" enrichment.sh --cells K562_essential --model binom_enrichment_with_covariates --chains 6 --threads 4

# K562_GenomeWide
sbatch --job-name="e-KG-fixed" enrichment.sh --cells K562_GenomeWide --model binom_enrichment_fixed --chains 6 --threads 4
sbatch --job-name="e-KG-varying" enrichment.sh --cells K562_GenomeWide --model binom_enrichment_varying_intercept --chains 6 --threads 4
sbatch --job-name="e-KG-cov" enrichment.sh --cells K562_GenomeWide --model binom_enrichment_with_covariates --chains 6 --threads 4

# HepG2
sbatch --job-name="e-H-fixed" enrichment.sh --cells HepG2 --model binom_enrichment_fixed --chains 6 --threads 4
sbatch --job-name="e-H-varying" enrichment.sh --cells HepG2 --model binom_enrichment_varying_intercept --chains 6 --threads 4
sbatch --job-name="e-H-cov" enrichment.sh --cells HepG2 --model binom_enrichment_with_covariates --chains 6 --threads 4

# RPE1
sbatch --job-name="e-R-fixed" enrichment.sh --cells RPE1 --model binom_enrichment_fixed --chains 6 --threads 4
sbatch --job-name="e-R-varying" enrichment.sh --cells RPE1 --model binom_enrichment_varying_intercept --chains 6 --threads 4
sbatch --job-name="e-R-cov" enrichment.sh --cells RPE1 --model binom_enrichment_with_covariates --chains 6 --threads 4
