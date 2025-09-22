library(tidyverse)
library(glue)
library(here)
library(org.Hs.eg.db)
library(data.table)
library(argparse)
source(here("src/functions_create_perturb_df.R"))
source(here("src/functions_ashr.R"))

parser <- ArgumentParser()
parser$add_argument("--ash", action = "store_true", help="Run adaptive shrinkage (ashr)")
parser$add_argument("--efficient", action = "store_true", help="Filter efficient perturbations only (>70% perturb)")
parser$add_argument("--unfiltered", action = "store_true", help="Don't filter perturbations-effect pairs for significance")
args <- parser$parse_args()

ash <- args$ash
unfiltered <- args$unfiltered
efficient <- args$efficient

###### PROCESS PERTURB-SEQ DEGS
perturb_DEG_dir <- here("data/perturb/DEGs/Jurkat_pseudo")

DEG_files <- list.files(perturb_DEG_dir, full.names = TRUE)

if(ash) {
  perturbation_effect_df <- apply_ash(DEG_files, cores=16, significant=!unfiltered)
} else {
  perturbation_effect_df <- extract_DEGs(DEG_files, cores=16, significant=!unfiltered)
}

effects.ensg <- perturbation_effect_df$effect %>% unique

effects.symbols <- mapIds(org.Hs.eg.db, keys = effects.ensg,
			  column = "SYMBOL", keytype = "ENSEMBL",
			  multiVals = "first")

effects.symbols <- effects.symbols %>% na.omit


perturbation_effect_df <- perturbation_effect_df %>%
  mutate(effect = effects.symbols[effect])

# Identify perturbations that appear as effects (i.e. they successfully perturb their targets)
if(efficient) {
  valid_perturbations <- perturbation_effect_df %>%
    filter(effect == perturbation) %>%
    mutate(percent_change = (2^avg_log2FC - 1)) %>%
    filter(percent_change < -0.7) %>% 
    pull(perturbation)
} else {
  valid_perturbations <- perturbation_effect_df %>%
    filter(effect == perturbation) %>%
    pull(perturbation)
}
# Keep only rows where perturbation is in the valid list
perturbation_effect_df <- perturbation_effect_df %>%
  filter(perturbation %in% valid_perturbations)

# filter out recursive edges
#perturbation_effect_df <- perturbation_effect_df %>% filter(perturbation != effect)

ash_prefix = if(ash) {"ash_"} else {""}
efficient_suffix = if(efficient) {"_eff"} else {""}
unfiltered_suffix = if(unfiltered) {"_all"} else {""}

# Write valid perturbation-effect pairs
fwrite(perturbation_effect_df, here(glue("data/perturb/pairs/Jurkat/{ash_prefix}bulk_pertubration_pairs{unfiltered_suffix}{efficient_suffix}.csv")))

####### PROCESS EQTLS
#### CIS-EQTLS
cis_eQTL <- fread(here("data/summary_stats/cis_eQTLs_eQTLgen.txt.gz"), sep = "\t")

cis_eQTL.perturbations.all <- cis_eQTL %>%
  filter(GeneSymbol %in% perturbation_effect_df$perturbation) %>%  # get the cis-eQTLs of perturbation genes
  filter(Pvalue < 5e-8) %>%  # only significant
  calculate_eQTL_beta  # calculate beta from z-score

# aggregate the eQTLs by SNP. So each row has a unique SNP ID along with lists of eQTLs (p-vals and gene symbols)
cis_eQTL.perturbations.aggregated <- cis_eQTL.perturbations.all[, lapply(.SD, function(x) list(x)), by = SNP]

#### TRANS-EQTLS
trans_eQTL <- fread(here("data/summary_stats/trans_eQTLs_eQTLgen.txt.gz"), sep = "\t")

# Get all trans-eQTLs (for SNPs that are significant cis-eQTLs, and their target gene is part of the affected genes in perturb-seq)
trans_eQTL.effect.unfiltered <- trans_eQTL %>%
  filter(GeneSymbol %in% perturbation_effect_df$effect & SNP %in% cis_eQTL.perturbations.all$SNP) %>%
  calculate_eQTL_beta

trans_eQTL.effect.unfiltered.aggregated <- trans_eQTL.effect.unfiltered[, lapply(.SD, function(x) list(x)), by = SNP]

merged.QTL.unfiltered <- merge(cis_eQTL.perturbations.aggregated, trans_eQTL.effect.unfiltered.aggregated, by = "SNP", suffixes = c(".perturb", ".effect"))

# this effectively filters out the cis/trans-only SNPs due to unnest_longer
# e.g. unnest_longer(c(Pvalue.perturb, GeneSymbol.perturb)) will drop rows where Pvalue.perturb is NA because there is nothing to expand on that row
QTL.pairs.unfiltered <- merged.QTL.unfiltered %>%
  unnest_longer(c(Beta.perturb, Pvalue.perturb, FDR.perturb, GeneSymbol.perturb)) %>%
  unnest_longer(c(Beta.effect, Pvalue.effect, FDR.effect, GeneSymbol.effect)) %>%
  rename_with(~ c("perturbation", "effect"), c(GeneSymbol.perturb, GeneSymbol.effect))

fwrite(QTL.pairs.unfiltered, file=here(glue("data/perturb/pairs/Jurkat/{ash_prefix}bulk_eQTL_pairs{unfiltered_suffix}{efficient_suffix}.csv")))

