#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(argparse)
  library(glue)
  library(data.table)
  library(cmdstanr)
})

# ----------------------------
# Enrichment models only
# ----------------------------
supported_models <- c(
  "binom_enrichment_fixed",                 # k ~ Binom(m, p)
  "binom_enrichment_varying_intercept",     # k ~ Binom(m, p_i), logit(p_i)=a+u_i
  "binom_enrichment_with_covariates"        # add per-perturb covariates + RE
)

essential_screens <- c("K562_essential", "Jurkat", "RPE1", "HepG2")
gw_screens <- c("K562_GenomeWide")

parser <- ArgumentParser()
parser$add_argument("--cells", type = "character", required = TRUE,
                    choices = c(essential_screens, gw_screens),
                    help = "Cell type/line")
parser$add_argument("--model", type = "character", required = TRUE,
                    choices = supported_models,
                    help = "Enrichment model type")

# optional knobs
parser$add_argument("--min_deg", type = "integer", default = 1,
                    help = "Minimum number of DEGs (x==1) per perturb (default 1).")
parser$add_argument("--chains", type = "integer", default = 4,
                    help = "Chains (default 4)")
parser$add_argument("--iter", type = "integer", default = 1500,
                    help = "Iterations per chain (default 1500)")
parser$add_argument("--warmup", type = "integer", default = 750,
                    help = "Warmup per chain (default 750)")
parser$add_argument("--threads", type = "integer", default = 1,
                    help = "Within-chain threads (default 1). Use 1–2 unless you know it helps.")

args <- parser$parse_args()
cells <- args$cells
model <- args$model

min_deg <- args$min_deg
chains_n <- args$chains
iter_n <- args$iter
warmup_n <- args$warmup
threads_n <- args$threads

print(glue("Cells: {cells} | Model: {model} | min_deg={min_deg} | chains={chains_n} iter={iter_n} warmup={warmup_n} threads={threads_n}"))

dat_file <- if (cells %in% essential_screens) {
  "/rds/project/rds-csoP2nj6Y6Y/biv22/data/pairs/full_dat.csv"
} else {
  "/rds/project/rds-csoP2nj6Y6Y/biv22/data/pairs/full_dat_GW.csv"
}

cmdstanr::set_cmdstan_path("/home/biv22/rds/hpc-work/.cmdstan/cmdstan-2.36.0")

# ----------------------------
# Load + filter
# ----------------------------
dat_raw <- fread(dat_file) %>%
  mutate(perturb_eff_percent = 1 - 2^perturb_eff) %>%
  # efficient perturbations only (or NA if target gene not detected)
  filter(screen == cells, (perturb_eff_percent >= 0.7 | is.na(perturb_eff_percent))) %>%
  select(perturb, effect, x, y, perturb_eff_percent)

# Filter to perturbs with at least min_deg DEGs
dat_edge <- dat_raw %>%
  group_by(perturb) %>%
  filter(sum(x == 1) >= min_deg) %>%
  ungroup()

# ----------------------------
# Construct per-perturb binomial data
# ----------------------------
# m = number of DEGs for perturb (x==1)
# k = number of overlaps among DEGs (x==1 & y==1)
dat_binom <- dat_edge %>%
  group_by(perturb) %>%
  summarise(
    m = sum(x == 1),
    k = sum(x == 1 & y == 1),

    # covariates computed among DEGs only
    mean_abs_effect_deg = ifelse(m > 0, mean(abs(effect[x == 1]), na.rm = TRUE), NA_real_),
    mean_effect_deg     = ifelse(m > 0, mean(effect[x == 1], na.rm = TRUE), NA_real_),

    # efficiency is ~constant within perturb; use mean for safety
    eff = mean(perturb_eff_percent, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(m >= min_deg)

if (nrow(dat_binom) == 0) {
  stop(glue("No perturbations left after filtering: min_deg={min_deg}."))
}

# ----------------------------
# Priors helper
# ----------------------------
logit <- function(p) log(p / (1 - p))
clamp <- function(x, lo, hi) pmin(pmax(x, lo), hi)

# overall overlap rate among DEGs
p0 <- sum(dat_binom$k) / sum(dat_binom$m)
p0 <- clamp(p0, 1e-6, 1 - 1e-6)
int_mu <- logit(p0)

# ----------------------------
# Choose model
# ----------------------------
formula <- NULL
priors <- NULL
dat_model <- dat_binom

if (model == "binom_enrichment_fixed") {
  # k ~ Binomial(m, p)
  formula <- bf(k | trials(m) ~ 1)
  priors <- c(
    prior(normal(int_mu, 1), class = "Intercept")
  )

} else if (model == "binom_enrichment_varying_intercept") {
  # k ~ Binomial(m, p_i), logit(p_i)=a+u_i
  formula <- bf(k | trials(m) ~ 1 + (1 | perturb))
  priors <- c(
    prior(normal(int_mu, 1), class = "Intercept"),
    prior(exponential(1), class = "sd", group = "perturb")
  )

} else if (model == "binom_enrichment_with_covariates") {
  # Add per-perturb predictors. Standardize for stability.
  dat_model <- dat_binom %>%
    mutate(
      log_m = log(pmax(m, 1)),
      z_log_m = as.numeric(scale(log_m)),
      z_mean_abs_effect_deg = as.numeric(scale(mean_abs_effect_deg)),
      z_eff = as.numeric(scale(eff))
    )

  formula <- bf(k | trials(m) ~ 1 + z_log_m + z_mean_abs_effect_deg + z_eff + (1 | perturb))
  priors <- c(
    prior(normal(int_mu, 1), class = "Intercept"),
    prior(normal(0, 0.5), class = "b"),
    prior(exponential(1), class = "sd", group = "perturb")
  )

} else {
  stop("Invalid enrichment model type selected.")
}

# ----------------------------
# Fit
# ----------------------------
threads_spec <- if (!is.null(threads_n) && threads_n > 1) threading(threads_n) else NULL
cores_n <- min(chains_n, parallel::detectCores())

fit <- brm(
  formula = formula,
  data = dat_model,
  family = binomial(link = "logit"),
  prior = priors,
  chains = chains_n,
  cores = cores_n,
  threads = threads_spec,
  iter = iter_n,
  warmup = warmup_n,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  backend = "cmdstanr",
  save_pars = save_pars(all = TRUE)
)

# ----------------------------
# Save
# ----------------------------
out_dir <- glue("/rds/project/rds-csoP2nj6Y6Y/biv22/models/{cells}")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_file <- glue("{out_dir}/enrichment_{model}_minDEG{min_deg}.rds")
saveRDS(fit, out_file)
print(glue("Saved: {out_file}"))
