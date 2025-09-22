library(tidyverse)
library(brms)
library(here)
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
parser$add_argument("--bulk", action = "store_true", help="Pseudobulk")
parser$add_argument("--efficient", action = "store_true", help="Efficient perturbations >70% only")

args <- parser$parse_args()
cells <- args$cells
model <- args$model
efficient <- args$efficient

ash <- if(args$ash) "ash_" else ""
bulk <- if(args$bulk) "bulk_" else ""
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

dat <- read_csv(here(glue("data/perturb/pairs/{bulk}{cells}/{ash}{bulk}logit_dat{eff}.csv")))

# Fit logistic regression
fit <- brm(
  formula = formula,
  data = dat,
  family = bernoulli(link = "logit"),
#  threads = threading(2),
  chains = 8,
  cores = 8,
  iter = 2000,
  backend = "cmdstanr"
)

saveRDS(fit, here(glue("models/{bulk}{cells}/N_{ash}logit_{model}{eff}.rds")))
