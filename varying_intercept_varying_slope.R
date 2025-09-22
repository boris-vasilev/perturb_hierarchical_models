library(tidyverse)
library(here)
library(brms)
library(argparse)
library(glue)

supported_models <- c("varying_intercept_varying_slope",
		      "no_pooling_intercept_varying_slope",
		      "varying_intercept_fixed_slope")

parser <- ArgumentParser()
parser$add_argument("--cells", type = "character", help = "Cell type/line", required = TRUE)
parser$add_argument("--model", type = "character", help = "Model type", required = TRUE, choices = supported_models)
args <- parser$parse_args()
cells <- args$cells
model <- args$model

print(glue("Cells: {cells}            Model: {model}"))

# Use if-else for selecting the formula based on the model
if (model == "varying_intercept_varying_slope") {
  formula <- bf(y ~ x + (1 + x | perturb))
} else if (model == "no_pooling_intercept_varying_slope") {
  formula <- bf(y ~ 0 + factor(perturb) + x + (0 + x | perturb))
} else if (model == "varying_intercept_fixed_slope") {
  formula <- bf(y ~ x + (1 | perturb))
} else {
  stop("Invalid model type selected")  # Handle unexpected cases
}

cmdstanr::set_cmdstan_path("/home/biv22/rds/hpc-work/.cmdstan/cmdstan-2.36.0")

# Read perturb-seq pairs
perturbation_pairs <- read_csv(here(glue("data/perturb/pairs/{cells}/pertubration_pairs.csv")))

# Read eQTL pairs
eQTL_pairs <- read_csv(here(glue("data/perturb/pairs/{cells}/eQTL_pairs.csv")))

# Identify the lead SNP with the smallest perturbation cis-eQTL p-value
lead_snp_pairs <- eQTL_pairs %>%
  group_by(perturbation) %>%
  filter(Pvalue.perturb == min(Pvalue.perturb)) %>%
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
  unique

# all possible pairs (cartesian product)
all_pairs <- expand.grid(perturbation = perturbation_pairs$perturbation %>% unique,
			 effect = perturbation_pairs$effect %>% unique)

# Add indicators to the original datasets
perturb_effect <- perturb_effect %>%
  mutate(perturb_pair = 1)

cis_trans <- cis_trans %>%
  mutate(eQTL_pair = 1)

# Join with all_pairs, defaulting to 0 if not found
pairs_df <- all_pairs %>%
  left_join(perturb_effect, by = c("perturbation", "effect")) %>%
  left_join(cis_trans, by = c("perturbation", "effect")) %>%
  mutate(
    perturb_pair = ifelse(is.na(perturb_pair), 0, 1),
    eQTL_pair = ifelse(is.na(eQTL_pair), 0, 1)
  ) %>%
  select(perturbation, effect, perturb_pair, eQTL_pair) %>%  # Keep only relevant columns
  filter(perturbation != effect) %>%
  filter(paste(perturbation, effect) %in% paste(eQTL_pairs$perturbation, eQTL_pairs$effect))


dat <- data.frame(x=pairs_df$perturb_pair,
		  y=pairs_df$eQTL_pair,
		  perturb=as.factor(pairs_df$perturbation),
		  effect=as.factor(pairs_df$effect))

# Fit logistic regression
fit <- brm(
  formula = formula,
  data = dat,
  family = bernoulli(link = "logit"),
  threads = threading(2),
  chains = 4,
  cores = 4,
  iter = 2000,
  prior=c(prior(normal(0, 1), class = "b")),
  #algorithm="meanfield",
  backend = "cmdstanr"
)

saveRDS(fit, here(glue("models/{cells}/logit_{model}.rds")))

