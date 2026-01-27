library(tidyverse)
library(data.table)
library(parallel)
library(ashr)
library(argparse)
source("functions_create_perturb_df.R")

parser <- ArgumentParser()
parser$add_argument("--cores", type = "numeric", help = "Number of cores", required = TRUE)

args <- parser$parse_args()

n_cores <- args$cores

data_dir <- "/rds/project/rds-csoP2nj6Y6Y/biv22/data"
pairs_dir <- file.path(data_dir, "pairs")

output_file <- file.path(pairs_dir, "full_eqtl_dat.csv")

####### PROCESS EQTLS
### -- Load and filter cis-eQTLs -- ###
message("[1/7] Reading cis-eQTLs")
cis_eQTL <- fread("/rds/project/rds-csoP2nj6Y6Y/biv22/data/eqtl/cis_eQTLs_eQTLgen.tsv",
                  sep = "\t",
                  select = c("SNP", "Pvalue", "FDR", "GeneSymbol", "Zscore", "NrSamples", "AssessedAllele", "OtherAllele"))
message("  ✔ cis-eQTLs loaded: ", nrow(cis_eQTL))

message("[2/7] Selecting significant cis-eSNPs and calculating beta")
setkey(cis_eQTL, GeneSymbol)
cis_eSNPs <- cis_eQTL %>%
  as.data.table %>%
  calculate_eQTL_beta %>%
  # group_by(GeneSymbol) %>%
  # filter(Pvalue == min(Pvalue)) %>%
  # filter(abs(Beta) == max(abs(Beta))) %>%
  # ungroup() %>%
  filter(FDR < 0.05)

message("  ✔ Significant perturbation cis-eSNPs found: ", nrow(cis_eSNPs))

message("[3/7] Reading trans-eQTLs")
trans_eQTL <- fread("/rds/project/rds-csoP2nj6Y6Y/biv22/data/eqtl/trans_eQTLs_eQTLgen.txt",
                    sep = "\t",
                    select = c("SNP", "Pvalue", "FDR", "GeneSymbol", "Zscore", "NrSamples", "AssessedAllele", "OtherAllele"))
message("  ✔ trans-eQTLs loaded: ", nrow(trans_eQTL))

message("[4/7] Selecting expressed genes trans-eQTLs and calculating beta")
setkey(trans_eQTL, GeneSymbol)
trans_eQTL <- trans_eQTL[SNP %in% cis_eSNPs$SNP] %>% calculate_eQTL_beta
message("  ✔ Matching trans-eQTLs: ", nrow(trans_eQTL))

message("[5/7] Merging cis and trans-eQTLs of the lead cis-eSNPs of perturbed genes")
merged.QTL <- merge(
  cis_eSNPs,
  trans_eQTL,
  allow.cartesian = TRUE,
  by = "SNP", suffixes = c(".perturb", ".effect")
) %>% rename(
  perturb = "GeneSymbol.perturb",
  effect = "GeneSymbol.effect",
) %>% mutate(y = ifelse(FDR.effect < 0.05, 1, 0))

message("  ✔ eQTL pairs: ", nrow(merged.QTL))

message("[6/7] Selecting lead cis-eSNP of perturbed gene")

merged.QTL <- merged.QTL %>%
  group_by(perturb) %>%
  filter(Pvalue.perturb == min(Pvalue.perturb)) %>%
  filter(abs(Beta.perturb) == max(abs(Beta.perturb))) %>%
  ungroup()

message("✔ Writing final outputs...")
fwrite(merged.QTL, file = output_file)
message("🎉 Done.")
