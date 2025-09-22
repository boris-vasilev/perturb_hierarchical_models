library(tidyverse)
library(here)
library(data.table)

# Read MAF from 1000 Genomes European ancestry
chr_frq_list <- list.files(path = ("/rds/project/rds-csoP2nj6Y6Y/biv22/data/1000G_Phase3_frq"), full.names = TRUE)
snp_frq_dt <- lapply(chr_frq_list, fread) %>% rbindlist

# Adjust MAF based on allele matching between eQTLgen and 1000G
adjust_MAF <- function(eqtl_dt) {
  # Inner join on SNP
  eqtl_dt <- merge(eqtl_dt, snp_frq_dt, by = "SNP")

  # Match alleles
  eqtl_dt[, MAF := fifelse(
    AssessedAllele == A1 & OtherAllele == A2, 
    MAF, 
    fifelse(AssessedAllele == A2 & OtherAllele == A1, 
            1 - MAF, 
            NA_real_)
  )]

  # Remove mismatched alleles
  eqtl_dt <- eqtl_dt[!is.na(MAF)]

  return(eqtl_dt)
}

# Compute effect sizes (beta) from Z-score using adjusted MAF
calculate_eQTL_beta <- function(eqtl_dt) {
  eqtl_dt <- adjust_MAF(eqtl_dt)

  # Compute beta
  eqtl_dt[, Beta := Zscore / sqrt(2 * MAF * (1 - MAF) * (NrSamples + Zscore^2))]

  # Return with useful columns
  return(eqtl_dt[, .(SNP, GeneSymbol, Pvalue, FDR, Beta)])
}

