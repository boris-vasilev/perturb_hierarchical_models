library(tidyverse)
library(here)
library(org.Hs.eg.db)
library(data.table)
library(argparse)

source(here("src/functions_create_perturb_df.R"))
source(here("src/functions_ashr.R"))

parser <- ArgumentParser()
parser$add_argument("--unfiltered", action = "store_true", help="Don't filter perturbations-effect pairs for significance")
parser$add_argument("--efficient", action = "store_true", help="Efficient perturbations >70% only")
args <- parser$parse_args()

unfiltered <- args$unfiltered

unfiltered_suffix = if(unfiltered) {"_all"} else {""}

###### PROCESS PERTURB-SEQ DEGS
perturb_DEG_dir <- here("data/perturb/DEGs/Jurkat")

DEG_files <- list.files(perturb_DEG_dir, full.names = TRUE)

significant_perturb_effects <- list()

# Loop through each file
for(file in DEG_files) {
  # Read the file
  DEGs <- read.table(file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  
  # Extract perturbed gene from filename
  perturbed_gene <- basename(file) %>% sub("\\.tsv$", "", .)
  
  if(unfiltered) {
    DEGs <- DEGs %>%
      rownames_to_column(var = "effect") %>%
      mutate(perturbation = perturbed_gene) %>%
      dplyr::select(effect, perturbation, p_val_adj, avg_log2FC)  # Keep relevant columns
  } else {
    DEGs <- DEGs %>%
      filter(p_val_adj < 0.05) %>%
      rownames_to_column(var = "effect") %>%
      mutate(perturbation = perturbed_gene) %>%
      dplyr::select(effect, perturbation, p_val_adj, avg_log2FC)  # Keep relevant columns
  }
  # Store dataframe in list
  significant_perturb_effects[[perturbed_gene]] <- DEGs
}

# Combine all dataframes into one
perturbation_effect_df <- bind_rows(significant_perturb_effects)

effects.ensg <- perturbation_effect_df$effect %>% unique

effects.symbols <- mapIds(org.Hs.eg.db, keys = effects.ensg,
			  column = "SYMBOL", keytype = "ENSEMBL",
			  multiVals = "first")

effects.symbols <- effects.symbols %>% na.omit


perturbation_effect_df <- perturbation_effect_df %>%
  mutate(effect = effects.symbols[effect])

# Identify perturbations that appear as effects (i.e. they successfully perturb their targets)
valid_perturbations <- perturbation_effect_df %>%
  filter(effect == perturbation) %>%
  pull(perturbation)

# Keep only rows where perturbation is in the valid list
perturbation_effect_df <- perturbation_effect_df %>%
  filter(perturbation %in% valid_perturbations)

# filter out recursive edges
#perturbation_effect_df <- perturbation_effect_df %>% filter(perturbation != effect)

# Write valid perturbation-effect pairs
fwrite(perturbation_effect_df, here("data/perturb/pairs/Jurkat/pertubration_pairs{unfiltered_suffix}.csv"))

quit(save = "no", status = 0)  # or q("no", status = 0)

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

fwrite(QTL.pairs.unfiltered, file=here("data/perturb/pairs/Jurkat/eQTL_pairs{unfiltered_suffix}.csv"))

###### OVERLAP OF QLT PAIRS (same SNP) AND PERTURB PAIRS

# all possible pairs of perturbations and affected genes. Cartesian product of perturbations and sequenced genes
all_pairs <- expand.grid(perturbation=names(significant_perturb_effects), effect=effects.symbols$hgnc_symbol)

# check if the pairs are true perturbation pairs, significant eQTL pair, or any eQTL pair (signficant cis-eQTL, any trans-eQTL)
all_pairs <- all_pairs %>%
  mutate(perturbseq = if_else(paste(perturbation, effect) %in% paste(perturbation_effect_df$perturbation, perturbation_effect_df$effect), TRUE, FALSE)) %>%
  mutate(eqtl_significant = if_else(paste(perturbation, effect) %in% paste(QTL.pairs$perturbation, QTL.pairs$effect), TRUE, FALSE)) %>%
  mutate(eqtl_any = if_else(paste(perturbation, effect) %in% paste(QTL.pairs.unfiltered$perturbation, QTL.pairs.unfiltered$effect), TRUE, FALSE))

all_pairs$eqtl <- case_when(
  all_pairs$eqtl_significant == TRUE & all_pairs$eqtl_any == TRUE ~ TRUE,   # Significant trans-eQTL
  all_pairs$eqtl_significant == FALSE & all_pairs$eqtl_any == TRUE ~ FALSE,  # Not significant trans-eQTL
  all_pairs$eqtl_significant == FALSE & all_pairs$eqtl_any == FALSE ~ NA    # No trans-eQTL at all
)

confusion_matrix <- table(data.frame(perturbseq=all_pairs$perturbseq,
                                     eqtl =  all_pairs$eqtl))

print(confusion_matrix)

fisher.test(confusion_matrix)

# A few perturbations affect a huge number of genes
# Might skew the results of the Fisher's test
number_of_effects_per_perturbation <-lapply(significant_perturb_effects, length) %>% unlist %>% sort

png("perturbation_effects_count.png")
number_of_effects_per_perturbation %>% hist(., xlab="Number of affected genes", main="Number of affected genes/perturbation")
dev.off()

png("perturbation_effects_count_log.png")
log10(number_of_effects_per_perturbation) %>% hist(., xlab="Log10(number of affected genes per perturbation)", main="Number of affected genes/perturbation (log-scale)")
dev.off()

##### JACCARD OVERLAP OF THE AFFECTED GENE SETS BETWEEN PERTURB AND EQTL

grouped_eQTL_sets <- QTL.pairs %>% group_by(perturbation) %>% summarize(effect = list(as.character(unique(effect))))
grouped_perturb_sets <- perturbation_effect_df %>% group_by(perturbation) %>% summarize(effect = list(as.character(unique(effect))))

grouped_eQTL_sets <- grouped_eQTL_sets %>% mutate(perturbation = as.character(perturbation))
grouped_perturb_sets <- grouped_perturb_sets %>% mutate(perturbation = as.character(perturbation))

library(purrr)

# Convert tibbles to named lists
eQTL_sets <- grouped_eQTL_sets %>% 
  deframe()  # Creates a named list: perturbation -> gene set

perturb_sets <- grouped_perturb_sets %>% 
  deframe()  # Creates a named list: perturbation -> gene set

# Function to compute Jaccard similarity
jaccard_index <- function(eQTL_set, perturb_set) {
  if (is.null(eQTL_set) || is.null(perturb_set)) return(NA) # Handle missing perturbations
  intersection <- length(intersect(eQTL_set, perturb_set))
  union <- length(unique(c(eQTL_set, perturb_set)))
  return(intersection / union)
}

eQTL_recoverty_rate <- function(eQTL_set, perturb_set) {
  if (is.null(eQTL_set) || is.null(perturb_set)) return(NA) # Handle missing perturbations
  intersection <- length(intersect(eQTL_set, perturb_set))
  perturb_set_size <- length(unique(perturb_set))
  return(intersection / perturb_set_size)
}

# Compute Jaccard index for perturbations found in both datasets
jaccard_scores <- tibble(
  perturbation = intersect(names(eQTL_sets), names(perturb_sets)),  # Only common perturbations
  jaccard = map2_dbl(eQTL_sets[perturbation], perturb_sets[perturbation], jaccard_index),
  recovery_rate = map2_dbl(eQTL_sets[perturbation], perturb_sets[perturbation], eQTL_recoverty_rate)
)

# View results
print(jaccard_scores, n = 20)


###### QQ PLOT OF PERTURBATION PAIR TRANS-EQTL P-VALUES

pvals.true <- QTL.pairs.unfiltered %>%
  inner_join(perturbation_effect_df, by = c("perturbation", "effect")) %>%
  pull(Pvalue.effect)

pvals.false <- QTL.pairs.unfiltered %>%
  anti_join(perturbation_effect_df, by = c("perturbation", "effect")) %>%
  pull(Pvalue.effect)

library(qqman)

png("qqplot_true.png")
pvals.true %>% qq(., main="Trans-eQTL of true perturb-seq affected genes")
dev.off()

png("qqplot_false.png")
pvals.false  %>% qq(., main="Trans-eQTL of false perturb-seq affected genes")
dev.off()

df <- data.frame(
  pvalue = c(pvals.true, pvals.false),
  group = rep(c("True Pairs", "False Pairs"), c(length(pvals.true), length(pvals.false)))
)

png("perturb_qtl_ecdf.png")
ggplot(df, aes(x = pvalue, color = group)) +
  stat_ecdf(geom = "step") +
  theme_minimal() +
  labs(title = "ECDF of p-values: True vs. False Pairs",
       x = "p-value", y = "Cumulative Probability")
dev.off()

beta.true <- QTL.pairs.unfiltered %>%
  inner_join(perturbation_effect_df, by = c("perturbation", "effect")) %>%
  pull(Beta)

beta.false <- QTL.pairs.unfiltered %>%
  anti_join(perturbation_effect_df, by = c("perturbation", "effect")) %>%
  pull(Beta)

png("trans_eqtl_betas.png")
ggplot() +
  geom_density(aes(x = log10(abs(beta.true))), color = "blue") +
  geom_density(aes(x = log10(abs(beta.false))), color = "red")
dev.off()
