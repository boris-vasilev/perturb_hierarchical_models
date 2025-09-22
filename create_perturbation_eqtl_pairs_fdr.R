library(tidyverse)
library(glue)
library(here)
library(org.Hs.eg.db)
library(data.table)
library(argparse)
source(here("src/functions_create_perturb_df.R"))
source(here("src/functions_ashr.R"))

parser <- ArgumentParser()
parser$add_argument("--cells", type = "character", help="Cell type: K562_essential, K562, Jurkat, RPE1", required=T)
parser$add_argument("--ash", action = "store_true", help="Run adaptive shrinkage (ashr)")
parser$add_argument("--bulk", action = "store_true", help="Pseudobulk?")
parser$add_argument("--efficient", action = "store_true", help="Filter efficient perturbations only (>70% perturb)")
parser$add_argument("--unfiltered", action = "store_true", help="Don't filter perturbations-effect pairs for significance")
args <- parser$parse_args()

cells <- args$cells
ash <- args$ash
bulk <- args$bulk
unfiltered <- args$unfiltered
efficient <- args$efficient

###### PROCESS PERTURB-SEQ DEGS

bulk_prefix <- if(bulk) "bulk_" else ""
perturb_DEG_dir <- here(glue("data/perturb/DEGs/{bulk_prefix}{cells}"))

DEG_files <- list.files(perturb_DEG_dir, full.names = TRUE)

if(ash) {
  perturbation_effect_df <- apply_ash(DEG_files, cores=16, significant=!unfiltered, bulk=bulk)
} else {
  results <- extract_DEGs(DEG_files,
                          cores=16,
                          significant=!unfiltered,
                          bulk=bulk,
                          efficient = efficient)
}

perturbation_effect_df <- results$perturbation_effect_df
filter_summary <- results$filter_summary

ash_prefix = if(ash) {"ash_"} else {""}
efficient_suffix = if(efficient) {"_eff"} else {""}
unfiltered_suffix = if(unfiltered) {"_all"} else {""}

# Write valid perturbation-effect pairs
fwrite(perturbation_effect_df, here(glue("data/perturb/pairs/{bulk_prefix}{cells}/{ash_prefix}{bulk_prefix}perturbation_pairs{unfiltered_suffix}{efficient_suffix}.csv")))

####### PROCESS EQTLS
#### CIS-EQTLS
print("Reading cis-eQTLs")
cis_eQTL <- fread(here("data/summary_stats/cis_eQTLs_eQTLgen_significant.txt.gz"), sep = "\t")
cis_eQTL <- cis_eQTL[, .(SNP, FDR, Pvalue, GeneSymbol, Zscore, NrSamples, AssessedAllele, OtherAllele)]

filter_summary$total_cis_eQTLs <- nrow(cis_eQTL)
filter_summary$total_cis_eQTL_genes <- cis_eQTL$GeneSymbol %>% unique %>% length

cis_eQTL.perturbations.all <- cis_eQTL %>%
  filter(GeneSymbol %in% perturbation_effect_df$perturbation) # get the cis-eQTLs of perturbation genes

filter_summary$cis_eQTLs_perturbed <- nrow(cis_eQTL.perturbations.all)
filter_summary$cis_eQTL_perturbed_genes <- cis_eQTL.perturbations.all$GeneSymbol %>% unique %>% length

cis_eQTL.perturbations.all <- cis_eQTL.perturbations.all %>%
  calculate_eQTL_beta  # calculate beta from z-score

# aggregate the eQTLs by SNP. So each row has a unique SNP ID along with lists of eQTLs (p-vals and gene symbols)
cis_eQTL.perturbations.aggregated <- cis_eQTL.perturbations.all[, lapply(.SD, function(x) list(x)), by = SNP]

#### TRANS-EQTLS
print("Reading trans-eQTLs")
trans_eQTL <- fread(here("data/summary_stats/trans_eQTLs_eQTLgen_significant.txt.gz"), sep = "\t")
trans_eQTL <- trans_eQTL[, .(SNP, FDR, Pvalue, GeneSymbol, Zscore, NrSamples, AssessedAllele, OtherAllele)]

filter_summary$total_trans_eQTLs <- nrow(trans_eQTL)
filter_summary$total_trans_eQTL_genes <- trans_eQTL$GeneSymbol %>% unique %>% length

# Get all trans-eQTLs (for SNPs that are significant cis-eQTLs, and their target gene is part of the affected genes in perturb-seq)
trans_eQTL.effect.unfiltered <- trans_eQTL %>%
  filter(GeneSymbol %in% perturbation_effect_df$effect)

filter_summary$trans_eQTL_DEG <- nrow(trans_eQTL.effect.unfiltered)
filter_summary$trans_eQTL_DEG_genes <- trans_eQTL.effect.unfiltered$GeneSymbol %>% unique %>% length

trans_eQTL.effect.unfiltered <- trans_eQTL.effect.unfiltered %>%
  filter(SNP %in% cis_eQTL.perturbations.all$SNP) %>%
  calculate_eQTL_beta

filter_summary$trans_eQTL_DEG_matching_cis <- nrow(trans_eQTL.effect.unfiltered)
filter_summary$trans_eQTL_DEG_matching_cis_genes <- trans_eQTL.effect.unfiltered$GeneSymbol %>% unique %>% length

trans_eQTL.effect.unfiltered.aggregated <- trans_eQTL.effect.unfiltered[, lapply(.SD, function(x) list(x)), by = SNP]

merged.QTL.unfiltered <- merge(cis_eQTL.perturbations.aggregated, trans_eQTL.effect.unfiltered.aggregated, by = "SNP", suffixes = c(".perturb", ".effect"))

# this effectively filters out the cis/trans-only SNPs due to unnest_longer
# e.g. unnest_longer(c(Pvalue.perturb, GeneSymbol.perturb)) will drop rows where Pvalue.perturb is NA because there is nothing to expand on that row
QTL.pairs.unfiltered <- merged.QTL.unfiltered %>%
  unnest_longer(c(Beta.perturb, Pvalue.perturb, FDR.perturb, GeneSymbol.perturb)) %>%
  unnest_longer(c(Beta.effect, Pvalue.effect, FDR.effect, GeneSymbol.effect)) %>%
  rename_with(~ c("perturbation", "effect"), c(GeneSymbol.perturb, GeneSymbol.effect))

filter_summary$eQTL_pairs <- nrow(QTL.pairs.unfiltered)
filter_summary$eQTL_perturbs <- QTL.pairs.unfiltered$perturbation %>% unique %>% length
filter_summary$eQTL_effects <- QTL.pairs.unfiltered$perturbation %>% unique %>% length

fwrite(QTL.pairs.unfiltered, file=here(glue("data/perturb/pairs/{bulk_prefix}{cells}/{ash_prefix}{bulk_prefix}eQTL_pairs{unfiltered_suffix}{efficient_suffix}.csv")))

saveRDS(filter_summary, file=here(glue("data/perturb/pairs/{bulk_prefix}{cells}/filter_summary.RDS")))



