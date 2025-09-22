library(tidyverse)
library(here)
library(brms)
library(argparse)
library(glue)
library(data.table)

supported_models <- c("varying_intercept_varying_slope",
                      "no_pooling_intercept_and_slope",
                      "varying_intercept_fixed_slope",
                      "vivs_student",
                      "vivs_student_concordance",
                      "vivs_concordance",
                      "vifs_concordance",
                      "vivs_student_pair_slope")

parser <- ArgumentParser()
parser$add_argument("--cells", type = "character", help = "Cell type/line", required = TRUE)
parser$add_argument("--response", type = "character", help = "Response type: trans-eQTL effect size or eQTL effect ratio", required = TRUE, choices = c("trans", "ratio"))
parser$add_argument("--model", type = "character", help = "Model type", required = TRUE, choices = supported_models)
parser$add_argument("--ash", action = "store_true", help="Run adaptive shrinkage (ashr)")
parser$add_argument("--bulk", action = "store_true", help="Pseudobulk")
parser$add_argument("--efficient", action = "store_true", help="Efficient")
parser$add_argument("--cpus", type= "numeric", help="Number of cores", default = 8)

args <- parser$parse_args()
cells <- args$cells
model <- args$model
efficient <- args$efficient
response <- args$response
cpus <- args$cpus

ash <- if(args$ash) "ash_" else ""
bulk <- if(args$bulk) "bulk_" else ""
eff <- if(args$efficient) "_eff" else ""

print(glue("Cells: {cells}             Response: {response}           Model: {model}           Efficient: {efficient}"))

dat <- fread(here(glue("data/perturb/pairs/{bulk}{cells}/{ash}{bulk}gaussian_{response}_dat{eff}.csv")))

# Make unsigned
#dat <- dat %>% mutate(x = abs(x), y = abs(y))

# Use if-else for selecting the formula based on the model
if (model == "varying_intercept_varying_slope") {
  formula <- bf(y ~ x + (1 + x | perturb))
} else if (model == "no_pooling_intercept_and_slope") {
  formula <- bf(y ~ 0 + perturb + x:perturb)
} else if (model == "vivs_student") {
  formula <- bf(y ~ 1 + x + (1 + x | gr(perturb, dist = "student")))
} else if (model == "vivs_student_concordance") {
  dat <- dat %>% mutate(c = sign(x) != sign(y))  # different signs means concordance. Check notes from 30/4 to see why!
  formula <- bf(c ~ x + (1 + x | gr(perturb, dist = "student")))
} else if (model == "vivs_concordance") {
  dat <- dat %>% mutate(c = sign(x) != sign(y))  # different signs means concordance. Check notes from 30/4 to see why!
  formula <- bf(c ~ x + (1 + x | perturb))
} else if (model == "vifs_concordance") {
  dat <- dat %>% mutate(c = sign(x) != sign(y))  # different signs means concordance. Check notes from 30/4 to see why!
  formula <- bf(c ~ x + (1 | perturb))
} else if (model == "vivs_student_pair_slope") {
  dat <- dat %>% mutate(pair = paste(perturb, effect))
  formula <- bf(y ~ 1 + x + (1 | gr(perturb, dist = "student")) + (x | gr(pair, dist = "student")))
} else if (model == "varying_intercept_fixed_slope") {
  formula <- bf(y ~ x + (1 | perturb))
} else {
  stop("Invalid model type selected")  # Handle unexpected cases
}

cmdstanr::set_cmdstan_path("/home/biv22/rds/hpc-work/.cmdstan/cmdstan-2.36.0")

if(!(model %in% c("vivs_concordance", "vivs_student_concordance", "vifs_concordance"))) {
# Fit logistic regression
fit <- brm(
  formula = formula,
  data = dat,
  family = gaussian(),
  #threads = threading(cpus/4 ),
  chains = 4,
  cores = 4,
  iter = 4000,
  backend = "cmdstanr"
)

saveRDS(fit, here(glue("models/{bulk}{cells}/N_{ash}gaussian_{response}_{model}{eff}.rds")))
} else {
fit <- brm(
  formula = formula,
  data = dat,
  family = bernoulli(link = "logit"),
 # threads = threading(cpus/4),
  chains = 4,
  cores = 4,
  iter = 8000,
  backend = "cmdstanr"
)
saveRDS(fit, here(glue("models/{bulk}{cells}/N_{ash}sign_{response}_{model}{eff}.rds")))
}
