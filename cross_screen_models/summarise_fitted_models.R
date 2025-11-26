library(tidyverse)
library(cmdstanr)
library(posterior)
library(loo)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--cores", type = "numeric", required = TRUE)
args <- parser$parse_args()
cores <- args$cores

options(mc.cores = cores)
set_cmdstan_path("/home/biv22/rds/hpc-work/.cmdstan/cmdstan-2.36.0")

fitted_dir <- "/rds/project/rds-csoP2nj6Y6Y/biv22/perturb_hierarchical_models/cross_screen_models/fitted"
out_dir <- "/rds/project/rds-csoP2nj6Y6Y/biv22/perturb_hierarchical_models/cross_screen_models/fitted_summary"
dir.create(out_dir, showWarnings = FALSE)

files <- list.files(fitted_dir, full.names = TRUE, pattern = "\\.RDS$")

message("Found ", length(files), " files.")

for (file in files) {
  
  fname <- basename(file)
  out_csv <- file.path(out_dir, paste0(fname, ".summary.csv"))
  
  if (file.exists(out_csv)) {
    message("Skipping ", fname, " (already done)")
    next
  }
  
  message("Processing ", fname)
  
  fit <- readRDS(file)
  
  # extract perturb and model_type 
  perturb <- sub("_(full_me|me_j|no_me)\\.RDS$", "", fname)
  model_type <- case_when(
    grepl("full_me", fname) ~ "full_me",
    grepl("me_j", fname) ~ "me_j",
    TRUE ~ "no_me"
  )
  
  # -------- diagnostics ----------
  summ <- fit$summary(validate_csv = FALSE)
  diag <- fit$diagnostic_summary(validate_csv = FALSE)
  
  rhat_max   <- max(summ$rhat, na.rm = TRUE)
  ess_min    <- min(summ$ess_bulk, na.rm = TRUE)
  n_diverg   <- sum(diag$divergent__)
  treedepth  <- sum(diag$max_treedepth__)
  ebfmi_min  <- min(diag$ebfmi__, na.rm = TRUE)
  
  # -------- PPC ----------
  dat <- fit$metadata()$data
  y_obs <- dat$y
  
  yrep <- fit$draws("y_rep", format = "matrix")
  ndraw <- min(30, nrow(yrep))
  yrep <- yrep[sample(1:nrow(yrep), ndraw), ]
  
  lo <- apply(yrep, 2, quantile, 0.05)
  hi <- apply(yrep, 2, quantile, 0.95)
  coverage_90 <- mean(y_obs > lo & y_obs < hi)
  
  rmse <- mean((rowMeans(t(yrep)) - y_obs)^2)
  
  # -------- Params ----------
  avail <- fit$metadata()$model_params
  sigma_param <- intersect(c("sigma", "sigma_pert"), avail)
  
  draws <- fit$draws(c("mu_delta", "tau", sigma_param), format="df")
  med <- summarise_all(draws, median)
  names(med)[names(med)==sigma_param] <- "sigma"
  
  # -------- LOO ----------
  ll <- fit$draws("log_lik")
  loo_obj <- loo(ll, cores = cores)
  
  summary_row <- tibble(
    file = fname,
    perturb = perturb,
    model_type = model_type,
    rhat_max = rhat_max,
    ess_min = ess_min,
    n_diverg = n_diverg,
    treedepth = treedepth,
    ebfmi_min = ebfmi_min,
    coverage_90 = coverage_90,
    rmse = rmse,
    mu_delta = med$mu_delta,
    tau = med$tau,
    sigma = med$sigma,
    loo_ic = loo_obj$estimates["looic","Estimate"],
    loo_se = loo_obj$estimates["looic","SE"]
  )
  
  write_csv(summary_row, out_csv)
  
  rm(fit, yrep, draws)
  gc()
}
