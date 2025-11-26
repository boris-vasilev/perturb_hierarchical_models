#!/usr/bin/Rscript

library(cmdstanr)

set_cmdstan_path("/home/biv22/rds/hpc-work/.cmdstan/cmdstan-2.36.0")

models_dir <- "/rds/project/rds-csoP2nj6Y6Y/biv22/perturb_hierarchical_models/cross_screen_models/models"
mod.no_me <- cmdstan_model(file.path(models_dir, "single_gene_no_me.stan"))
mod.me_j <- cmdstan_model(file.path(models_dir, "single_gene_me_j.stan"))
mod.me_both <- cmdstan_model(file.path(models_dir, "single_gene_full_me.stan"))