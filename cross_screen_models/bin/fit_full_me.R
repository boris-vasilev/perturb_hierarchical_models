#!/usr/bin/Rscript

library(argparse)
library(glue)
library(tidyverse)
library(data.table)
library(cmdstanr)

parser <- ArgumentParser()
parser$add_argument("--perturb", type = "character", required = TRUE)

args <- parser$parse_args()
i <- args$perturb

print(glue("Fitting FULL ME model for {i}"))

set_cmdstan_path("/home/biv22/rds/hpc-work/.cmdstan/cmdstan-2.36.0")

models_dir <- "/rds/project/rds-csoP2nj6Y6Y/biv22/perturb_hierarchical_models/cross_screen_models/models"
fitted_dir <- "/rds/project/rds-csoP2nj6Y6Y/biv22/perturb_hierarchical_models/cross_screen_models/fitted"

dat <- fread("/rds/project/rds-csoP2nj6Y6Y/biv22/data/pairs/full_dat.csv")

dat_i <- dat %>%
  filter(perturb == i) %>%
  drop_na(logFC, lfcSE, perturb_eff, perturb_eff_se, screen, effect) %>%
  mutate(
    scaled_logFC = logFC / perturb_eff,
    scaled_SE    = lfcSE / abs(perturb_eff)
  )

effects <- sort(unique(dat_i$effect))
screens <- sort(unique(dat_i$screen))

J <- length(effects)
K <- length(screens)

dat_i <- dat_i %>%
  mutate(
    j_idx = match(effect, effects),
    k_idx = match(screen, screens)
  )

perturb_eff_df <- dat_i %>%
  select(screen, perturb_eff, perturb_eff_se) %>%
  distinct() %>%
  arrange(screen)

beta_i_obs <- perturb_eff_df$perturb_eff
se_beta_i  <- perturb_eff_df$perturb_eff_se

N <- nrow(dat_i)

stan_data.me_both <- list(
  N       = N,
  J       = J,
  K       = K,
  j_idx   = dat_i$j_idx,
  k_idx   = dat_i$k_idx,
  beta_j_obs = dat_i$logFC,
  se_beta_j  = dat_i$lfcSE,
  beta_i_obs    = beta_i_obs,
  se_beta_i     = se_beta_i
)
if(!file.exists("{fitted_dir}/{i}_full_me.RDS")) {
  mod.me_both <- cmdstan_model(stan_file=file.path(models_dir, "single_gene_full_me.stan"),
                               exe_file=file.path(models_dir, "single_gene_full_me"),
                               compile = FALSE)
  fit.me_both <- mod.me_both$sample(
    data = stan_data.me_both,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    seed = 123
  )

  fit.me_both$save_object(file = glue("{fitted_dir}/{i}_full_me.RDS"))
}