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
parser$add_argument("--cells", type = "character", help = "Cell type/line", required = TRUE, choices = c(essential_screens, gw_screens))
parser$add_argument("--model", type = "character", help = "Model type", required = TRUE, choices = supported_models)
# parser$add_argument("--efficient", action = "store_true", help="Efficient perturbations >70% only")

args <- parser$parse_args()
cells <- args$cells
model <- args$model
# efficient <- args$efficient

# eff <- if(args$efficient) "_eff" else ""

dat_file <- if(args$cells %in% essential_screens) {
  "/rds/project/rds-csoP2nj6Y6Y/biv22/data/pairs/full_dat.csv"
} else {
  "/rds/project/rds-csoP2nj6Y6Y/biv22/data/pairs/full_dat_GW.csv"
}

print(glue("Cells: {cells}            Model: {model} "))

# Use if-else for selecting the formula based on the model
if (model %in% c("vivs", "vivs_horseshoe")) {
  formula <- bf(y ~ x + (1 + x | perturb))
} else if (model == "no_pooling_intercept_varying_slope") {
  formula <- bf(y ~ 0 + factor(perturb) + x + (0 + x | perturb))
} else if (model == "vivs_student"){
  formula <- bf(y ~ x + (1 + x | gr(perturb, dist="student")))
} else if (model == "varying_intercept_fixed_slope") {
  formula <- bf(y ~ x + (1 | perturb))
} else {
  stop("Invalid model type selected")  # Handle unexpected cases
}

cmdstanr::set_cmdstan_path("/home/biv22/rds/hpc-work/.cmdstan/cmdstan-2.36.0")

dat <- fread(dat_file) %>%
  mutate(perturb_eff_percent = 1 - 2^perturb_eff) %>%
  # Filter for efficient perturbations. Either efficiency >= 70% or NA (target gene not detected)
  filter(screen == cells, (perturb_eff_percent >= 0.7 | is.na(perturb_eff_percent))) %>%
  select(perturb, effect, x, y) %>%
  group_by(perturb) %>%
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

# saveRDS(fit, glue("/rds/project/rds-csoP2nj6Y6Y/biv22/models/{cells}/logit_{model}{eff}.rds"))
saveRDS(fit, glue("/rds/project/rds-csoP2nj6Y6Y/biv22/models/{cells}/logit_{model}.rds"))
