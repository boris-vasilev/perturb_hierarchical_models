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

print(glue("Fitting models for {i}"))

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

####################### MODEL 1: NO ME #######################
stan_data.no_me <- list(
  N    = N,
  J    = J,
  j_idx = dat_i$j_idx,
  y    = dat_i$scaled_logFC
)

if(!file.exists("{fitted_dir}/{i}_no_me.RDS")) {
  print(glue("NO ME model for {i}"))
  mod.no_me <- cmdstan_model(file.path(models_dir, "single_gene_no_me.stan"))
  fit.no_me <- mod.no_me$sample(
    data = stan_data.no_me,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    seed = 123
  )
  saveRDS(fit.no_me, glue("{fitted_dir}/{i}_no_me.RDS"))
}

####################### MODEL 2: ME ON J #######################

stan_data.me_j <- list(
  N    = N,
  J    = J,
  j_idx = dat_i$j_idx,
  y    = dat_i$scaled_logFC,
  se_y    = dat_i$scaled_SE
)


if(!file.exists("{fitted_dir}/{i}_me_j.RDS")) {
  print(glue("ME ON J model for {i}"))
  mod.me_j <- cmdstan_model(file.path(models_dir, "single_gene_me_j.stan"))
  fit.me_j <- mod.me_j$sample(
    data = stan_data.me_j,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    seed = 123
  )

  saveRDS(fit.me_j, glue("{fitted_dir}/{i}_me_j.RDS"))
}

####################### MODEL 3: FULL ME #######################

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
  print(glue("FULL ME model for {i}"))
  mod.me_both <- cmdstan_model(file.path(models_dir, "single_gene_full_me.stan"))
  fit.me_both <- mod.me_both$sample(
    data = stan_data.me_both,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    seed = 123
  )

  saveRDS(fit.me_both, glue("{fitted_dir}/{i}_full_me.RDS"))
}