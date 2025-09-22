library(tidyverse)
library(here)
library(data.table)

# Read MAF from 1000 Genomes European ancestry
chr_frq_list <- list.files(path = here("LDSC/LDSCORE/1000G_Phase3_frq"), full.names = TRUE)
snp_frq_dt <- lapply(chr_frq_list, fread) %>% rbindlist

## Function to fix mismatches between eQTLgen and 1000G effect allele (AssessedAllele != A1)
#adjust_MAF <- function(eqtl_df) {
#  # Inner join the (cis/trans)-eQTL data from eQTLgen with the 1000G data
#  eqtl_df <- eqtl_df %>% inner_join(snp_frq, by="SNP")
#  # Loop over the rows of the dataframe to adjust MAF in place
#  eqtl_df <- eqtl_df %>%
#    mutate(
#      # If the alleles match, keep the MAF as is
#      MAF = case_when(
#        AssessedAllele == A1 & OtherAllele == A2 ~ MAF,                      # matching alleles
#        AssessedAllele == A2 & OtherAllele == A1 ~ 1 - MAF,                    # swapped alleles
#        TRUE ~ NA_real_                                                        # completely mismatched alleles
#      )
#    ) %>%
#    # Remove rows with completely mismatched alleles (MAF = NA)
#    filter(!is.na(MAF))
#
#  return(eqtl_df)
#}
#
## Function to calculate the effect size of eQTLs from the Z-score
#calculate_eQTL_beta <- function(eqtl_df) {
#  eqtl_df <- adjust_MAF(eqtl_df)
#  # As calculated in the SNP/HEIDI paper
#  eqtl_df <- eqtl_df %>%
#    mutate(Beta = Zscore / sqrt(2*MAF*(1-MAF)*(NrSamples + Zscore^2))) 
#
#  eqtl_df <- eqtl_df[, .(SNP, GeneSymbol, Pvalue, FDR)]
#
#  return(eqtl_df)
#}

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

