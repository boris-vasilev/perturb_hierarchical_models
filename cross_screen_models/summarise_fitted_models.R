library(tidyverse)
library(cmdstanr)
library(posterior)
library(loo)
library(argparse)
library(parallel)

message("=== Starting model diagnostics pipeline ===")

parser <- ArgumentParser()
parser$add_argument("--cores", type = "numeric", required = TRUE)

args <- parser$parse_args()
cores <- args$cores

options(mc.cores = cores)

#-------------------------------------------------------------
# Load all fits
#-------------------------------------------------------------

fitted_dir <- "/rds/project/rds-csoP2nj6Y6Y/biv22/perturb_hierarchical_models/cross_screen_models/fitted"

message("Listing files in ", fitted_dir)

files <- list.files(fitted_dir, full.names = TRUE, pattern = "\\.RDS$")
message("Found ", length(files), " fitted model files.")
message("Loading fits into memory...")

fits <- mclapply(files, readRDS, mc.cores = cores)

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
    fit = fits
  )

message("All fits loaded successfully.\n")

#-------------------------------------------------------------
# Convergence diagnostics
#-------------------------------------------------------------

message("=== Running convergence diagnostics ===")

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

diag_list <- mclapply(fits_table$fit, extract_diags, mc.cores = cores)

diag_table <- fits_table %>%
  mutate(diag = diag_list) %>%
  unnest(diag)

message("Convergence diagnostics computed for all models.\n")

diag_summary <- diag_table %>%
  group_by(model_type) %>%
  summarise(
    bad_rhat = mean(rhat_max > 1.01),
    bad_div  = mean(n_diverg > 0),
    bad_ess  = mean(ess_min < 200),
    bad_ebfmi = mean(ebfmi_min < 0.2)
  )

message("=== Convergence Summary ===")
print(diag_summary)
message("\n")

# Warn if any serious issues:
if (any(diag_table$rhat_max > 1.05)) {
  message("WARNING: Some fits have Rhat > 1.05")
}
if (any(diag_table$n_diverg > 0)) {
  message("WARNING: Divergences detected in some fits")
}

#-------------------------------------------------------------
# Posterior predictive summaries
#-------------------------------------------------------------

message("=== Computing PPC summaries ===")

ppc_summary <- function(fit, y_obs_varname = "y", yrep_varname = "y_rep") {
  
  dat <- fit$metadata()$data
  y_obs <- dat[[y_obs_varname]]
  
  yrep <- fit$draws(yrep_varname, format = "matrix")
  ndraws <- min(30, nrow(yrep))
  yrep <- yrep[sample(1:nrow(yrep), ndraws), ]
  
  lo <- apply(yrep, 2, quantile, 0.05)
  hi <- apply(yrep, 2, quantile, 0.95)
  coverage <- mean(y_obs > lo & y_obs < hi)
  
  rmse <- mean((apply(yrep, 2, mean) - y_obs)^2)
  
  tibble(
    coverage_90 = coverage,
    rmse = rmse
  )
}

ppc_list <- mclapply(fits_table$fit, ppc_summary, mc.cores = cores)

ppc_table <- fits_table %>%
  mutate(ppc = ppc_list) %>%
  unnest(ppc)

ppc_summary_by_model <- ppc_table %>%
  group_by(model_type) %>%
  summarise(
    median_cov = median(coverage_90),
    frac_under = mean(coverage_90 < 0.8),
    median_rmse = median(rmse)
  )

message("=== PPC Summary ===")
print(ppc_summary_by_model)
message("\n")

#-------------------------------------------------------------
# LOO IC
#-------------------------------------------------------------

message("=== Computing LOO ===")

compute_loo <- function(fit) {
  ll <- fit$draws("log_lik")
  loo::loo(ll)
}

loo_list <- mclapply(fits_table$fit, compute_loo, mc.cores = cores)

# LOO can be slow — print progress
loo_table <- fits_table %>%
  mutate(loo = loo_list)

loo_compare_table <- loo_table %>%
  select(perturb, model_type, loo) %>%
  unnest_wider(loo)

message("LOO computed for all models.\n")

#-------------------------------------------------------------
# Parameter summaries
#-------------------------------------------------------------

message("=== Extracting parameter summaries ===")

extract_params <- function(fit) {
  draws <- fit$draws(c("mu_delta", "tau", "sigma", "sigma_pert"), format="df")
  summarise_all(draws, median)
}

params_list <- mclapply(fits_table$fit, extract_params, mc.cores = cores)

params_table <- fits_table %>%
  mutate(params = params_list) %>%
  unnest(params)

message("Parameter summaries extracted.\n")

#-------------------------------------------------------------
# Save main summary
#-------------------------------------------------------------

output_file <- "/rds/project/rds-csoP2nj6Y6Y/biv22/perturb_hierarchical_models/cross_screen_models/fit_summary.csv"

message("Saving model summary to: ", output_file)

summary_table <- fits_table %>%
  select(perturb, model_type) %>%
  bind_cols(
    diag_table %>% select(-fit, -file, -filename),
    ppc_table %>% select(-fit),
    params_table %>% select(-fit)
  )
write_csv(summary_table, output_file)

message("=== All diagnostics complete. ===")
