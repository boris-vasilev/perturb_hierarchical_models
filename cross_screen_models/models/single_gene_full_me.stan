data {
  int<lower=1> N;                       // number of (k,j) observations
  int<lower=1> J;                       // number of target genes j
  int<lower=1> K;                       // number of screens k

  array[N] int<lower=1, upper=J> j_idx; // gene index for each row
  array[N] int<lower=1, upper=K> k_idx; // screen index for each row

  vector[N] beta_j_obs;                   // observed logFC_ij
  vector[N] se_beta_j;                    // lfcSE_ij

  vector[K] beta_i_obs;                      // observed perturb_eff (logFC_ii)
  vector[K] se_beta_i;                       // lfcSE_ii
}

parameters {
  // Causal effect structure
  real mu_delta;                        // global mean effect across genes
  real<lower=1e-6> tau;                // SD of delta_j
  vector[J] z_delta;                  // non-centered random effects

  // Residual noise beyond DESeq2 SE
  real<lower=1e-6> sigma_pert;         // extra noise across screens/genes

  // True perturbation efficiencies per screen
  vector[K] beta_i_true;                     // latent true logFC_ii
}

transformed parameters {
  vector[J] delta;                      // actual delta_j
  delta = mu_delta + tau * z_delta;
}

model {
  // Priors
  mu_delta   ~ normal(0, 0.1);          // small effects
  tau        ~ exponential(1);          // most genes near mu_delta
  sigma_pert ~ exponential(1);          // extra noise small but >0

  z_delta  ~ normal(0, 1);            // non-centered

  // Measurement error for perturbation efficiency
  beta_i_true ~ normal(beta_i_obs, se_beta_i);         // beta_i_obs is noisy measure of beta_i_true

  // Likelihood: measurement error on beta, noise + beta_i_true * delta
  for (n in 1:N) {
    real mean_beta_j = beta_i_true[ k_idx[n] ] * delta[ j_idx[n] ];
    beta_j_obs[n] ~ normal(mean_beta_j, sqrt(se_beta_j[n]^2 + sigma_pert^2));
  }
}

generated quantities {
  vector[N] beta_rep;
  for (n in 1:N) {
    real mean_beta_j = beta_i_true[ k_idx[n] ] * delta[ j_idx[n] ];
    beta_rep[n] = normal_rng(mean_beta_j, sqrt(se_beta_j[n]^2 + sigma_pert^2));
  }
}
