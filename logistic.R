library(tidyverse)
library(brms)
library(argparse)
library(glue)
library(data.table)


supported_models <- c("vivs",
                      "no_pooling_intercept_varying_slope",
                      "varying_intercept_fixed_slope",
                      "vivs_student",
                      "vivs_horseshoe")

essential_screens <- c("K562_essential", "Jurkat", "RPE1", "HepG2")
gw_screens <- c("K562_GenomeWide")

parser <- ArgumentParser()
parser$add_argument("--screen", type = "character", help = "Cell type/line", required = TRUE, choices = c(essential_screens, gw_screens))
parser$add_argument("--model", type = "character", help = "Model type", required = TRUE, choices = supported_models)
# parser$add_argument("--efficient", action = "store_true", help="Efficient perturbations >70% only")

args <- parser$parse_args()
screen <- args$screen
model <- args$model
# efficient <- args$efficient

# eff <- if(args$efficient) "_eff" else ""

dat_file <- glue("/rds/project/rds-csoP2nj6Y6Y/biv22/data/pairs/dat_{screen}.csv")

print(glue("Screen: {screen}            Model: {model} "))

# Use if-else for selecting the formula based on the model
if (model %in% c("vivs", "vivs_horseshoe")) {
  formula <- bf(y ~ x + (1 + x | cis_gene))
} else if (model == "no_pooling_intercept_varying_slope") {
  formula <- bf(y ~ 0 + factor(cis_gene) + x + (0 + x | cis_gene))
} else if (model == "vivs_student"){
  formula <- bf(y ~ x + (1 + x | gr(cis_gene, dist="student")))
} else if (model == "varying_intercept_fixed_slope") {
  formula <- bf(y ~ x + (1 | cis_gene))
} else {
  stop("Invalid model type selected")  # Handle unexpected cases
}

cmdstanr::set_cmdstan_path("/home/biv22/rds/hpc-work/.cmdstan/cmdstan-2.36.0")

dat <- fread(dat_file) %>%
  mutate(perturb_eff_percent = 1 - 2^perturb_eff) %>%
  # Filter for efficient perturbations. Either efficiency >= 70% or NA (target gene not detected)
  filter(perturb_eff_percent >= 0.7 | is.na(perturb_eff_percent)) %>%
  select(cis_gene, trans_gene, x, y) %>%
  group_by(cis_gene) %>%
  filter(any(x == 1)) %>%
  filter(any(y == 1)) %>%
  ungroup

prior <- if(model == "vivs_horseshoe") {
  c(
    prior(horseshoe(df = 1, par_ratio = 0.1), class = "b"),
    prior(exponential(1), class = "sd", group = "perturb"),
    prior(lkj(2), class = "cor", group = "perturb")
  )
} else NULL

fit <- brm(
  formula = formula,
  data = dat,
  family = bernoulli(link = "logit"),
  prior = prior,
  chains = 6,
  cores = 6,
  threads = threading(4),  # 6 × 4 = 24 total threads
  iter = 2000,
  warmup = 1000,
  control = list(
    adapt_delta = 0.98,
    max_treedepth = 15
  ),
  backend = "cmdstanr"
)

saveRDS(fit, glue("/rds/project/rds-csoP2nj6Y6Y/biv22/models/{screen}/logit_{model}.rds"))
