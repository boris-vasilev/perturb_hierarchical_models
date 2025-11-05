library(tidyverse)
library(brms)
library(argparse)
library(glue)
library(data.table)


supported_models <- c("varying_intercept_varying_slope",
                      "no_pooling_intercept_varying_slope",
                      "varying_intercept_fixed_slope",
                      "vivs_student")

parser <- ArgumentParser()
parser$add_argument("--cells", type = "character", help = "Cell type/line", required = TRUE)
parser$add_argument("--model", type = "character", help = "Model type", required = TRUE, choices = supported_models)
parser$add_argument("--ash", action = "store_true", help="Run adaptive shrinkage (ashr)")
parser$add_argument("--efficient", action = "store_true", help="Efficient perturbations >70% only")

args <- parser$parse_args()
cells <- args$cells
model <- args$model
efficient <- args$efficient

ash <- if(args$ash) "ash_" else ""
eff <- if(args$efficient) "_eff" else ""

print(glue("Cells: {cells}            Model: {model}           Efficient: {efficient} "))

# Use if-else for selecting the formula based on the model
if (model == "varying_intercept_varying_slope") {
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

dat <- fread(glue("/rds/project/rds-csoP2nj6Y6Y/biv22/data/pairs/{cells}/{ash}logit_dat{eff}.csv"))
# Fit logistic regression

fit <- brm(
  formula = formula,
  data = dat,
  family = bernoulli(link = "logit"),
  prior = c(
    prior(horseshoe(df = 1, par_ratio = 0.1), class = "b"),
    prior(exponential(1), class = "sd", group = "perturb"),
    prior(lkj(2), class = "cor", group = "perturb")
  ),
  chains = 6,
  cores = 6,
  threads = threading(4),  # 6 × 4 = 24 total threads
  iter = 4000,
  warmup = 2000,
  control = list(
    adapt_delta = 0.98,
    max_treedepth = 15
  ),
  backend = "cmdstanr"
)


saveRDS(fit, glue("/rds/project/rds-csoP2nj6Y6Y/biv22/models/{cells}/{ash}logit_{model}{eff}_HS.rds"))
