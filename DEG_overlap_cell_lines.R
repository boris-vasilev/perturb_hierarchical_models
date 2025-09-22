library(tidyverse)
library(here)
library(biomaRt)
library(data.table)

###### PROCESS PERTURB-SEQ DEGS
perturb_DEG_dir <- here("data/perturb/DEGs/Jurkat")

DEG_files <- list.files(perturb_DEG_dir, full.names = TRUE)

significant_perturb_effects <- list()

for(file in DEG_files) {
  DEGs.ensg <- read.table(file ,header = TRUE, sep = "\t", row.names = 1, check.names = FALSE) %>%
	  filter(p_val_adj < 0.05) %>% rownames
  perturbed_gene <- basename(file) %>% sub("\\.tsv$", "", .)
  significant_perturb_effects[[perturbed_gene]] <- DEGs.ensg
}

perturbation_effect_df <- stack(significant_perturb_effects)
colnames(perturbation_effect_df) <- c("effect", "perturbation")

effects.ensg <- perturbation_effect_df$effect %>% unique

mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
#mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="https://jan2024.archive.ensembl.org")

attributes <- c("ensembl_gene_id", "hgnc_symbol")

effects.symbols <- getBM(
  attributes = attributes,
  filters = "ensembl_gene_id",
  values = effects.ensg,
  mart = mart
)

# TODO: use external_synonym where no symbol is available and use lates ensembl when it comes back to life
effects.symbols <- effects.symbols %>% filter(hgnc_symbol != "")

perturbation_effect_df <- perturbation_effect_df %>%
  left_join(effects.symbols, by = c("effect" = "ensembl_gene_id")) %>%
  filter(!is.na(hgnc_symbol)) %>%  # Remove unmatched rows
  mutate(effect = hgnc_symbol) %>% # Replace effect with hgnc_symbol
  dplyr::select(-hgnc_symbol)  # Remove the extra column

# Identify perturbations that appear as effects (i.e. they successfully perturb their targets)
valid_perturbations <- perturbation_effect_df %>%
  filter(effect == perturbation) %>%
  pull(perturbation)

# Keep only rows where perturbation is in the valid list
perturbation_effect_df <- perturbation_effect_df %>%
  filter(perturbation %in% valid_perturbations)

# filter out recursive edges
perturbation_effect_df.Jurkat <- perturbation_effect_df %>% filter(perturbation != effect)


###### PROCESS PERTURB-SEQ DEGS (RPE-1)
perturb_DEG_dir <- here("data/perturb/DEGs/RPE1")

DEG_files <- list.files(perturb_DEG_dir, full.names = TRUE)

significant_perturb_effects <- list()

for(file in DEG_files) {
  DEGs.ensg <- read.table(file ,header = TRUE, sep = "\t", row.names = 1, check.names = FALSE) %>%
	  filter(p_val_adj < 0.05) %>% rownames
  perturbed_gene <- basename(file) %>% sub("\\.tsv$", "", .)
  significant_perturb_effects[[perturbed_gene]] <- DEGs.ensg
}

perturbation_effect_df <- stack(significant_perturb_effects)
colnames(perturbation_effect_df) <- c("effect", "perturbation")

effects.ensg <- perturbation_effect_df$effect %>% unique

effects.symbols <- getBM(
  attributes = attributes,
  filters = "ensembl_gene_id",
  values = effects.ensg,
  mart = mart
)

# TODO: use external_synonym where no symbol is available and use lates ensembl when it comes back to life
effects.symbols <- effects.symbols %>% filter(hgnc_symbol != "")

perturbation_effect_df <- perturbation_effect_df %>%
  left_join(effects.symbols, by = c("effect" = "ensembl_gene_id")) %>%
  filter(!is.na(hgnc_symbol)) %>%  # Remove unmatched rows
  mutate(effect = hgnc_symbol) %>% # Replace effect with hgnc_symbol
  dplyr::select(-hgnc_symbol)  # Remove the extra column

# Identify perturbations that appear as effects (i.e. they successfully perturb their targets)
valid_perturbations <- perturbation_effect_df %>%
  filter(effect == perturbation) %>%
  pull(perturbation)

# Keep only rows where perturbation is in the valid list
perturbation_effect_df <- perturbation_effect_df %>%
  filter(perturbation %in% valid_perturbations)

# filter out recursive edges
perturbation_effect_df.RPE1 <- perturbation_effect_df %>% filter(perturbation != effect)


# JACCARD OVERLAP OF THE AFFECTED GENE SETS BETWEEN PERTURB AND EQTL

grouped_perturb_sets.Jurkat <- perturbation_effect_df.Jurkat %>% group_by(perturbation) %>% summarize(effect = list(as.character(unique(effect))))
grouped_perturb_sets.RPE1 <- perturbation_effect_df.RPE1 %>% group_by(perturbation) %>% summarize(effect = list(as.character(unique(effect))))

grouped_perturb_sets.Jurkat <- grouped_perturb_sets.Jurkat %>% mutate(perturbation = as.character(perturbation))
grouped_perturb_sets.RPE1 <- grouped_perturb_sets.RPE1 %>% mutate(perturbation = as.character(perturbation))

library(purrr)

# Convert tibbles to named lists
perturb_sets.Jurkat <- grouped_perturb_sets.Jurkat %>%
  deframe()  # Creates a named list: perturbation -> gene set

perturb_sets.RPE1 <- grouped_perturb_sets.RPE1 %>%
  deframe()  # Creates a named list: perturbation -> gene set

# Function to compute Jaccard similarity
jaccard_index <- function(cell1_set, cell2_set) {
  if (is.null(cell1_set) || is.null(cell2_set)) return(NA) # Handle missing perturbations
  intersection <- length(intersect(cell1_set, cell2_set))
  union <- length(unique(c(cell1_set, cell2_set)))
  return(intersection / union)
}

# Compute Jaccard index for perturbations found in both datasets
jaccard_scores <- tibble(
  perturbation = intersect(names(perturb_sets.Jurkat), names(perturb_sets.RPE1)),  # Only common perturbations
  jaccard = map2_dbl(perturb_sets.Jurkat[perturbation], perturb_sets.RPE1[perturbation], jaccard_index),
)

# View results
print(jaccard_scores, n = 20)
