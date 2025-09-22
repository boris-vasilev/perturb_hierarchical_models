library(tidyverse)
library(here)
library(data.table)
library(argparse)
library(meta)
library(glue)
library(broom)
library(logistf)
library(brms)

parser <- ArgumentParser()
parser$add_argument("--cells", type = "character", help = "Cell type/line", required = TRUE)
parser$add_argument("--ash", action = "store_true", help="Run adaptive shrinkage (ashr)")
parser$add_argument("--bulk", action = "store_true", help="Pseudobulk")
parser$add_argument("--efficient", action = "store_true", help="Efficient perturb")

args <- parser$parse_args()
cells <- args$cells

ash <- if(args$ash) "ash_" else ""
bulk <- if(args$bulk) "bulk_" else ""
eff <- if(args$efficient) "_eff" else ""

print(glue("Contingency tables for Cells: {cells}"))

dat <- fread(here(glue("data/perturb/pairs/{bulk}{cells}/{ash}{bulk}logit_dat{eff}.csv")))

# Build contingency table per perturb
# Continuity correction (add +0.5 to cells with 0 counts)

# Step 1: Build contingency tables
contingency_data <- dat %>%
  group_by(perturb) %>%
  summarise(
    N_11 = sum(x == 1 & y == 1),
    N_10 = sum(x == 1 & y == 0),
    N_01 = sum(x == 0 & y == 1),
    N_00 = sum(x == 0 & y == 0),
    .groups = "drop"
  )

# Step 2: Add 0.5 if necessary (continuity correction)
contingency_data <- contingency_data %>%
  rowwise() %>%
  mutate(
    needs_correction = (N_11 == 0 | N_10 == 0 | N_01 == 0),
    N_11_corr = ifelse(needs_correction, N_11 + 0.5, N_11),
    N_10_corr = ifelse(needs_correction, N_10 + 0.5, N_10),
    N_01_corr = ifelse(needs_correction, N_01 + 0.5, N_01),
    N_00_corr = ifelse(needs_correction, N_00 + 0.5, N_00)
  ) %>%
  ungroup()

# Step 3: Calculate ORs
or_results <- contingency_data %>%
  mutate(
    OR = (N_11_corr * N_00_corr) / (N_10_corr * N_01_corr),
    log_OR = log(OR),
    SE_log_OR = sqrt(1/N_11_corr + 1/N_10_corr + 1/N_01_corr + 1/N_00_corr)
  )

meta_res <- metafor::rma(yi = log_OR, sei = SE_log_OR, method = "REML", data = or_results)

png(here(glue("plots/{cells}/{ash}{bulk}contingency_forest_plot{eff}.png")), height=3000, width=1500)
forest(meta_res, transf=exp, refline=1, xlab="Odds Ratio")
dev.off()

png(here(glue("plots/{cells}/{ash}{bulk}contingency_log_forest_plot{eff}.png")), height=3000, width=1500)
forest(meta_res, transf=identity, refline=1, xlab="log(OR)")
dev.off()

# Compute the study-specific effect sizes (i.e. perturbation-specific)
# Firth's penalised regression (normal logistic regression fails because of the 0 counts in some N_01, N_10 cells)
# another option is continuity correction by adding 0.5 to those cells
perturb_effects <- dat %>%
  group_by(perturb) %>%
  group_split() %>%
  map_df(function(df) {
    fit <- tryCatch(
      logistf(y ~ x, data = df),
      error = function(e) return(NULL)  # skip problematic fits
    )
    
    if (is.null(fit)) return(NULL)
    
    # Extract the coefficient for 'x'
    coef_index <- which(names(fit$coefficients) == "x")
    
    tibble(
      perturb = unique(df$perturb),
      estimate = fit$coefficients[coef_index],
      std.error = sqrt(diag(vcov(fit)))[coef_index],  # SE not provided in the fit so have to extract it from the covariance matrix
      p.value = fit$prob[coef_index],
      conf.low = fit$ci.lower[coef_index],
      conf.high = fit$ci.upper[coef_index],
      odds_ratio = exp(fit$coefficients[coef_index])
    ) 
  })

or_results %>%
  ggplot(aes(x = SE_log_OR, y = log_OR)) +
  geom_point() +
  labs(x = expression(sigma[italic(j)]~("log-odds")),
       y = expression(italic(y[j])~("log-odds")))

perturb_effects %>%
  ggplot(aes(x = std.error, y = estimate)) +
  geom_point() +
  geom_text(
    data = subset(perturb_effects, std.error < 1),
    aes(label = perturb, color="red"),
    vjust = -0.5,
    size = 3,
    fontface = "bold"
  ) +
  labs(
    x = expression(sigma[italic(j)]~("log-odds")),
    y = expression(italic(y[j])~("log-odds"))
  )


#cmdstanr::set_cmdstan_path("/home/biv22/rds/hpc-work/.cmdstan/cmdstan-2.36.0")
## Bayesian meta-analysis
#print("Fitting meta-analysis")
#me0 <- 
#  brm(data = perturb_effects, 
#      family = gaussian,
#      estimate | se(std.error) ~ 1 + (1 | perturb),
#      prior = c(prior(normal(0, 1.5), class = Intercept),
#                prior(exponential(1), class = sd)),
#      iter = 2000, warmup = 1000, cores = 4, chains = 4,
#      seed = 15, backend = "cmdstanr")
#saveRDS(me0, here(glue("models/{cells}/OR_meta_model.rds")))
#
#
#print("Fitting hierarchical model")
## Hierarchical model
#me1 <- 
#  brm(data = dat, 
#      family = binomial,
#      y | trials(1) ~ 0 + Intercept + x + (1 + x | perturb),
#      prior = c(prior(normal(0, 1.5), class = b),
#                prior(exponential(1), class = sd),
#                prior(lkj(2), class = cor)),
#      iter = 2000, warmup = 1000, cores = 4, chains = 4,
#      seed = 15, backend = "cmdstanr")
#
#saveRDS(me1, here(glue("models/{cells}/OR_hierarchical_model.rds")))

# Frequentist meta-analysis with Firth correction
me_frq <- metabin(
  event.e = N_11,
  n.e = N_11 + N_10,
  event.c = N_01,
  n.c = N_01 + N_00,
  studlab = perturb,
  data = contingency_data,
  sm = "OR",
  #allstudies = T,
  method = "LRP"
)

png(here(glue("plots/{cells}/{ash}{bulk}frq_forest_plot{eff}.png")), height=7000, width=1000)
forest(me_frq, sortvar=TE)
dev.off()

png(here(glue("plots/{cells}/{ash}{bulk}frq_drapery_plot{eff}.png")), height=3000, width=1000)
drapery(me_frq)
dev.off()
