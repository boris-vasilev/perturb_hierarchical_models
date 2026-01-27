library(tidyverse)
library(data.table)
library(parallel)

library(argparse)

parser <- ArgumentParser()
parser$add_argument("--cis_genes", type = "character", help = "File containing cis genes", required = TRUE)
parser$add_argument("--trans_genes", type = "character", help = "File containing trans genes", required = TRUE)

data_dir <- "/rds/project/rds-csoP2nj6Y6Y/biv22/data"
pairs_dir <- file.path(data_dir, "coloc")

cis_output_file <- file.path(pairs_dir, "x1y1_cis_eqtls.csv")
trans_output_file <- file.path(pairs_dir, "x1y1_trans_eqtls.csv")
merged_output_file <- file.path(pairs_dir, "x1y1_merged_eqtls.csv")

args <- parser$parse_args()
file_cis_genes <- args$cis_genes
file_trans_genes <- args$trans_genes


cis_genes <- readLines(file_cis_genes)
trans_genes <- readLines(file_trans_genes)

message("[1/7] Reading cis-eQTLs")
cis_eQTL <- fread("/rds/project/rds-csoP2nj6Y6Y/biv22/data/eqtl/cis_eQTLs_eQTLgen.tsv",
                  sep = "\t",
                  select = c("SNP", "Pvalue", "FDR", "GeneSymbol", "Zscore", "NrSamples", "AssessedAllele", "OtherAllele"))
message("  ✔ cis-eQTLs loaded: ", nrow(cis_eQTL))

cis_eQTL <- cis_eQTL %>% filter(GeneSymbol %in% cis_genes)

message("[3/7] Reading trans-eQTLs")
trans_eQTL <- fread("/rds/project/rds-csoP2nj6Y6Y/biv22/data/eqtl/trans_eQTLs_eQTLgen.txt",
                    sep = "\t",
                    select = c("SNP", "Pvalue", "FDR", "GeneSymbol", "Zscore", "NrSamples", "AssessedAllele", "OtherAllele"))
message("  ✔ trans-eQTLs loaded: ", nrow(trans_eQTL))

trans_eQTL <- trans_eQTL %>% filter(GeneSymbol %in% trans_genes)

merged_eQTL <- merge(cis_eQTL, trans_eQTL,
by = "SNP",
allow.cartesian = TRUE,
suffixes = c(".cis", ".trans"))

fwrite(cis_eQTL, file = cis_output_file)
fwrite(trans_eQTL, file = trans_output_file)
fwrite(merged_eQTL, file = merged_output_file)