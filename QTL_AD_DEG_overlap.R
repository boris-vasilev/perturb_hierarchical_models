# Overlap of cis/trans pairs and DEGs based on the Anderson-Darling test from Replogle et al.
# in K562 and RPE1

library(tidyverse)
library(here)
library(biomaRt)
library(data.table)
library(rhdf5)

# Read GWPS pseudobulk expression data
# Not used directly for analysis just to rename the perturbations (row names) of the AD test statistic file from ENSG to symbol
perturb.bulk.k562 <- H5Fopen(here('data/perturb/seurat/K562_gwps_raw_singlecell.h5ad'))

gene_name <- factor(h5read(perturb.bulk.k562, "/var/gene_name"), labels=h5read(perturb.bulk.k562, "/var/__categories/gene_name")) %>% as.character
gene_id <- h5read(perturb.bulk.k562, "/var/gene_id")

H5Fclose(perturb.bulk.k562)

expressed_gene_name2symbol_mapping <- data.frame(row.names = gene_id, symbol = gene_name)
expressed_gene_name2symbol_mapping["ENSG00000284024", "symbol"] <- "MSANTD7"

# Read GWPS AD values and fix gene names for row and column names
# There are 11258 perturbations (columns)
# and 5529 sequenced genes (rows)
# The format of the data is a (affected gene X perturbation gene) matrix

perturb_AD <- fread(here('data/summary_stats/andersondarlingpvaluesBHcorrected.csv.gz'))

# Ensure the symbol column is excluded from the matrix portion
gene_names <- expressed_gene_name2symbol_mapping[perturb_AD$V1, "symbol"]
perturb_AD <- perturb_AD %>% dplyr::select(-V1)
colnames(perturb_AD) <- sapply(strsplit(colnames(perturb_AD), "_"), function(x) x[2])
rownames(perturb_AD) <- gene_names

# Create a named list of affected genes for each perturbation gene
perturbations.affected_genes <- lapply(seq_len(ncol(perturb_AD)), function(j) {
  as.character(gene_names[perturb_AD[[j]] < 0.05])  # Use the gene names vector and filter by < 0.05 for the AD statistic
})
names(perturbations.affected_genes) <- colnames(perturb_AD)  # Assign column (perturbation) names to the list names

perturbation_effect_df <- stack(perturbations.affected_genes)
colnames(perturbation_effect_df) <- c("effect", "perturbation")

# Identify perturbations that appear as effects (i.e. they successfully perturb their targets)
valid_perturbations <- perturbation_effect_df %>%
  filter(effect == perturbation) %>%
  pull(perturbation)

# Keep only rows where perturbation is in the valid list
perturbation_effect_df <- perturbation_effect_df %>%
  filter(perturbation %in% valid_perturbations)

# filter out recursive edges
perturbation_effect_df <- perturbation_effect_df %>% filter(perturbation != effect)

####### PROCESS EQTLS
cis_eQTL <- fread(here("data/summary_stats/cis_eQTLs_eQTLgen.txt.gz"), sep = "\t")

# get the cis-eQTLs of perturbation genes
cis_eQTL.perturbations.all <- cis_eQTL %>% filter(GeneSymbol %in% perturbation_effect_df$perturbation)

# get only the significant cis-eQTLs of perturbations
cis_eQTL.perturbations.all <- cis_eQTL.perturbations.all %>%
  filter(Pvalue < 5e-8)

cis_eQTL.perturbations.all <- cis_eQTL.perturbations.all %>%
  mutate(Beta = Zscore / sqrt(NrSamples)) %>%
  dplyr::select(SNP, Beta, Pvalue, GeneSymbol)

cis_eQTL.perturbations.beta.filter <- cis_eQTL.perturbations.all %>% filter(abs(Beta) > 0.2) 

trans_eQTL <- fread(here("data/summary_stats/trans_eQTLs_eQTLgen.txt.gz"), sep = "\t")

# get the trans-eQTLs of affected genes
trans_eQTL.effect.all <- trans_eQTL %>% filter(GeneSymbol %in% perturbation_effect_df$effect)

# get only the significant trans-eQTLs of affected genes
trans_eQTL.effect.all <- trans_eQTL.effect.all %>%
  filter(Pvalue < 5e-8)

trans_eQTL.effect.all <- trans_eQTL.effect.all %>%
  mutate(Beta = Zscore / sqrt(NrSamples)) %>%
  dplyr::select(SNP, Beta, Pvalue, GeneSymbol)

# aggregate the eQTLs by SNP. So each row has a unique SNP ID along with lists of eQTLs (p-vals and gene symbols) 
cis_eQTL.perturbations.aggregated <- cis_eQTL.perturbations.all[, lapply(.SD, function(x) list(x)), by = SNP]
cis_eQTL.perturbations.beta.filter.aggregated <- cis_eQTL.perturbations.beta.filter[, lapply(.SD, function(x) list(x)), by = SNP]
trans_eQTL.effect.aggregated <- trans_eQTL.effect.all[, lapply(.SD, function(x) list(x)), by = SNP]


# full outer join of cis and trans-eQTLs. Will contain SNPs that act only in cis and only in trans
merged.QTL <- merge(cis_eQTL.perturbations.aggregated, trans_eQTL.effect.aggregated, by = "SNP", all=T, suffixes = c(".perturb", ".effect"))

# this effectively filters out the cis/trans-only SNPs due to unnest_longer
# e.g. unnest_longer(c(Pvalue.perturb, GeneSymbol.perturb)) will drop rows where Pvalue.perturb is NA because there is nothing to expand on that row
QTL.pairs <- merged.QTL %>%
  unnest_longer(c(Beta.perturb, Pvalue.perturb, GeneSymbol.perturb)) %>%
  unnest_longer(c(Beta.effect, Pvalue.effect, GeneSymbol.effect)) %>%
  dplyr::select(GeneSymbol.perturb, GeneSymbol.effect) %>%
  rename("perturbation" = "GeneSymbol.perturb", "effect" = "GeneSymbol.effect") %>%
  unique()

# Get all trans-eQTLs (for SNPs that are significant cis-eQTLs, and their target gene is part of the affected genes in perturb-seq)
trans_eQTL.effect.unfiltered <- trans_eQTL %>%
  filter(GeneSymbol %in% perturbation_effect_df$effect & SNP %in% cis_eQTL.perturbations.all$SNP) %>%
  mutate(Beta = Zscore / sqrt(NrSamples)) %>%
  dplyr::select(SNP, Beta, Pvalue, GeneSymbol)

# Note: This trans_eQTL.effect.unfiltered is a mix of both true perturbation pairs, and pairs that are not matched
# Here we split them into true perturbation pairs, and false perturbation pairs

trans_eQTL.effect.unfiltered.aggregated <- trans_eQTL.effect.unfiltered[, lapply(.SD, function(x) list(x)), by = SNP]

merged.QTL.unfiltered <- merge(cis_eQTL.perturbations.aggregated, trans_eQTL.effect.unfiltered.aggregated, by = "SNP", suffixes = c(".perturb", ".effect"))
merged.QTL.beta.filtered <- merge(cis_eQTL.perturbations.beta.filter.aggregated, trans_eQTL.effect.unfiltered.aggregated, by = "SNP", suffixes = c(".perturb", ".effect"))

QTL.pairs.unfiltered <- merged.QTL.unfiltered %>%
  unnest_longer(c(Beta.perturb, Pvalue.perturb, GeneSymbol.perturb)) %>%
  unnest_longer(c(Beta.effect, Pvalue.effect, GeneSymbol.effect)) %>%
  rename("perturbation" = "GeneSymbol.perturb", "effect" = "GeneSymbol.effect")

QTL.pairs.beta.filtered <- merged.QTL.beta.filtered %>%
  unnest_longer(c(Beta.perturb, Pvalue.perturb, GeneSymbol.perturb)) %>%
  unnest_longer(c(Beta.effect, Pvalue.effect, GeneSymbol.effect)) %>%
  rename("perturbation" = "GeneSymbol.perturb", "effect" = "GeneSymbol.effect")

###### OVERLAP OF QLT PAIRS (same SNP) AND PERTURB PAIRS

# all possible pairs of perturbations and affected genes. Cartesian product of perturbations and sequenced genes
all_pairs <- expand.grid(perturbation=perturbation_effect_df$perturbation %>% unique, effect=perturbation_effect_df$effect %>% unique)

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

# Histogram of Jaccard indexes
png("k562_overlap_jaccard_hist.png")
jaccard_scores$jaccard %>% hist(., main="K562 Histogram of Jaccard indexes")
dev.off()

# Histogram of recovery rate
png("k562_overlap_recovery_rate_hist.png")
jaccard_scores$recovery_rate %>% hist(., main="K562 Histogram of recovery rate")
dev.off()

###### QQ PLOT OF PERTURBATION PAIR TRANS-EQTL P-VALUES

pvals.true <- QTL.pairs.unfiltered %>%
  inner_join(perturbation_effect_df, by = c("perturbation", "effect")) %>%
  pull(Pvalue.effect)

pvals.false <- QTL.pairs.unfiltered %>%
  anti_join(perturbation_effect_df, by = c("perturbation", "effect")) %>%
  pull(Pvalue.effect)

library(qqman)

png("ad_qqplot_true.png")
pvals.true %>% qq(., main="Trans-eQTL of true perturb-seq affected genes")
dev.off()

png("ad_qqplot_false.png")
pvals.false %>% qq(., main="Trans-eQTL of false perturb-seq affected genes")
dev.off()

df <- data.frame(
  pvalue = c(pvals.true, pvals.false),
  group = rep(c("True Pairs", "False Pairs"), c(length(pvals.true), length(pvals.false)))
)

png("ad_perturb_qtl_ecdf.png")
ggplot(df, aes(x = -log10(pvalue), color = group)) +
  stat_ecdf(geom = "step") +
  theme_minimal() +
  labs(title = "ECDF of p-values: True vs. False Pairs",
       x = "p-value", y = "Cumulative Probability")
dev.off()

beta.true <- QTL.pairs.unfiltered %>%
  inner_join(perturbation_effect_df, by = c("perturbation", "effect")) %>%
  mutate(Beta = Beta.effect / Beta.perturb) %>%
  pull(Beta)

beta.false <- QTL.pairs.unfiltered %>%
  anti_join(perturbation_effect_df, by = c("perturbation", "effect")) %>%
  mutate(Beta = Beta.effect / Beta.perturb) %>%
  pull(Beta)

beta.true.filtered <- QTL.pairs.beta.filtered %>%
  inner_join(perturbation_effect_df, by = c("perturbation", "effect")) %>%
  mutate(Beta = Beta.effect / Beta.perturb) %>%
  pull(Beta)

beta.false.filtered <- QTL.pairs.beta.filtered %>%
  anti_join(perturbation_effect_df, by = c("perturbation", "effect")) %>%
  mutate(Beta = Beta.effect / Beta.perturb) %>%
  pull(Beta)

png("ad_trans_eqtl_betas.png")
ggplot() +
  geom_density(aes(x = log10(abs(beta.true)), color = "True pair")) +
  geom_density(aes(x = log10(abs(beta.false)), color = "False pair")) +
  xlab("log10(abs(beta))") +
  scale_color_manual(name = "Pair Type", values = c("True pair" = "blue", "False pair" = "red")) +
  theme_classic()
dev.off()


signed_log <- function(x, c = 1e-2) {
  sign(x) * log10(abs(x) + c)
}


png("ad_trans_eqtl_betas_filtered.png")
ggplot() +
  geom_density(aes(x = log10(abs(beta.true.filtered)), color = "True pair")) +
  geom_density(aes(x = log10(abs(beta.false.filtered)), color = "False pair")) +
  xlab("log10(abs(beta))") +
  scale_color_manual(name = "Pair Type", values = c("True pair" = "blue", "False pair" = "red")) +
  theme_classic()
dev.off()

png("ad_trans_eqtl_betas_filtered_signed.png")
ggplot() +
  geom_density(aes(x = signed_log(beta.true.filtered), color = "True pair")) +
  geom_density(aes(x = signed_log(beta.false.filtered), color = "False pair")) +
  xlab("signed log10(abs(beta))") +
  scale_color_manual(name = "Pair Type", values = c("True pair" = "blue", "False pair" = "red")) +
  theme_classic()
dev.off()
