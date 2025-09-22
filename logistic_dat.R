library(tidyverse)
library(argparse)
library(glue)
library(data.table)


parser <- ArgumentParser()
parser$add_argument("--cells", type = "character", help = "Cell type/line", required = TRUE)
parser$add_argument("--ash", action = "store_true", help="Run adaptive shrinkage (ashr)")
parser$add_argument("--efficient", action = "store_true", help="Only efficient >70% perturbations")

args <- parser$parse_args()
cells <- args$cells


print(glue("Generating data file: Cells: {cells}"))

ash <- if(args$ash) "ash_" else ""
eff <- if(args$efficient) "_eff" else ""

pairs_dir <- file.path("/rds/project/rds-csoP2nj6Y6Y/biv22/data/pairs", cells)

# Read perturb-seq pairs
perturbation_pairs <- fread(file.path(pairs_dir, glue("{ash}perturbation_pairs{eff}.csv")))

# Read eQTL pairs
eQTL_pairs <- fread(file.path(pairs_dir, glue("{ash}eQTL_pairs{eff}.csv")))

# Identify the lead SNP with the smallest perturbation cis-eQTL p-value
# In case of SNPs with the same  cis-eQTL p-value pick the one with the biggest Beta in magnitude
lead_snp_pairs <- eQTL_pairs %>%
  group_by(perturbation) %>%
  filter(Pvalue.perturb == min(Pvalue.perturb)) %>%
  filter(abs(Beta.perturb) == max(abs(Beta.perturb))) %>%
  ungroup()

# Select unique perturb-seq pairs
perturb_effect <- perturbation_pairs %>%
  select(perturbation, effect) %>%
  filter(perturbation != effect) %>%
  unique

# Select only the pairs with significant trans-eQTLs
cis_trans <- lead_snp_pairs %>%
  filter(FDR.effect < 0.05) %>%
  select(perturbation, effect) %>%
  filter(perturbation != effect) %>%
  filter(effect != "") %>% # filter out genes mapped to empty gene symbol
  unique

# all possible pairs (cartesian product)
all_pairs <- expand.grid(perturbation = perturb_effect$perturbation %>% unique,
                         effect = perturb_effect$effect %>% unique)

pairs_df <- all_pairs %>%
  mutate(perturbation = as.character(perturbation),
         effect = as.character(effect)) %>%
  filter(perturbation != effect, # remove self-pairs
         paste(perturbation, effect) %in% paste(lead_snp_pairs$perturbation, lead_snp_pairs$effect)) %>%  # keep only pairs of genes that have been tested both for cis-eQTL on the perturbation and trans-eQTL on the effect
  mutate(perturb_pair = ifelse(paste(perturbation,
                                     effect) %in% paste(perturb_effect$perturbation,
                                     perturb_effect$effect), 1, 0),
         eQTL_pair = ifelse(paste(perturbation,
                                  effect) %in%
                            paste(cis_trans$perturbation,
                                  cis_trans$effect), 1, 0))

dat <- data.frame(x=pairs_df$perturb_pair,
                  y=pairs_df$eQTL_pair,
                  perturb=as.factor(pairs_df$perturbation),
                  effect=as.factor(pairs_df$effect))

# filter perturbations with no x=1 (no perturbation pair). Those appear because the perturbation pair was not tested for both cis-eQTL and trans-eQTL and was excluded when contructing pairs_df

dat <- dat %>%
  group_by(perturb) %>%
  filter(any(x == 1)) %>%
  filter(any(y == 1)) %>%
  ungroup()


fwrite(dat, file.path(pairs_dir, glue("{ash}logit_dat{eff}.csv")))

