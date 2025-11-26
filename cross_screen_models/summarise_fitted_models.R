library(tidyverse)
library(cmdstanr)
library(posterior)
library(loo)

fitted_dir <- "/rds/project/rds-csoP2nj6Y6Y/biv22/perturb_hierarchical_models/cross_screen_models/fitted"

files <- list.files(fitted_dir, full.names = TRUE, pattern = "\\.RDS$")

fits_table <- tibble(
  file = files,
  filename = basename(files)
) %>%
  mutate(
    perturb = sub("_(full_me|me_j|no_me)\\.RDS$", "", filename),
    model_type = case_when(
      grepl("full_me\\.RDS$", filename) ~ "full_me",
      grepl("me_j\\.RDS$",    filename) ~ "me_j",
      grepl("no_me\\.RDS$",   filename) ~ "no_me",
      TRUE ~ NA_character_
    ),
    fit = map(file, readRDS)
  )

####################### Convergence diagnostics #######################
extract_diags <- function(fit) {
  summ <- fit$summary()
  diag <- fit$diagnostic_summary()
  
  tibble(
    rhat_max   = max(summ$rhat, na.rm = TRUE),
    ess_min    = min(summ$ess_bulk, na.rm = TRUE),
    n_diverg   = sum(diag$divergent__),
    treedepth  = sum(diag$max_treedepth__),
    ebfmi_min  = min(diag$ebfmi__, na.rm = TRUE)
  )
}

diag_table <- fits_table %>%
  mutate(diag = map(fit, extract_diags)) %>%
  unnest(diag)

diag_summary <- diag_table %>%
  group_by(model_type) %>%
  summarise(
    bad_rhat = mean(rhat_max > 1.01),
    bad_div  = mean(n_diverg > 0),
    bad_ess  = mean(ess_min < 200),
    bad_ebfmi = mean(ebfmi_min < 0.2)
  )
####################### Posterior predictive check summaries #######################
ppc_summary <- function(fit, y_obs_varname = "y", yrep_varname = "y_rep") {
  
  # Extract observed y
  dat <- fit$metadata()$data
  y_obs <- dat[[y_obs_varname]]
  N <- length(y_obs)
  
  # Extract y_rep
  yrep <- fit$draws(yrep_varname, format = "matrix")
  ndraws <- min(30, nrow(yrep))
  yrep <- yrep[sample(1:nrow(yrep), ndraws), ]
  
  # 90% coverage
  lo <- apply(yrep, 2, quantile, 0.05)
  hi <- apply(yrep, 2, quantile, 0.95)
  coverage <- mean(y_obs > lo & y_obs < hi)
  
  # RMSE
  rmse <- mean((apply(yrep, 2, mean) - y_obs)^2)
  
  tibble(
    coverage_90 = coverage,
    rmse = rmse
  )
}

ppc_table <- fits_table %>%
  mutate(ppc = map(fit, ppc_summary)) %>%
  unnest(ppc)

ppc_summary_by_model <- ppc_table %>%
  group_by(model_type) %>%
  summarise(
    median_cov = median(coverage_90),
    frac_under = mean(coverage_90 < 0.8),
    median_rmse = median(rmse)
  )

####################### LOO and WAIC #######################
compute_loo <- function(fit) {
  ll <- fit$draws("log_lik")
  loo(ll)
}
loo_table <- fits_table %>%
  mutate(loo = map(fit, compute_loo))

loo_compare_table <- loo_table %>%
  select(perturb, model_type, loo) %>%
  unnest_wider(loo)

####################### Effect summaries #######################
extract_params <- function(fit) {
  draws <- fit$draws(c("mu_delta", "tau", "sigma", "sigma_pert"), format="df")
  summarise_all(draws, median)
}
params_table <- fits_table %>%
  mutate(params = map(fit, extract_params)) %>%
  unnest(params)

fwrite(fits_table, "/rds/project/rds-csoP2nj6Y6Y/biv22/perturb_hierarchical_models/cross_screen_models/fit_summary.csv")