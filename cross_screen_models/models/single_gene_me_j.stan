data {
  int<lower=1> N;                       // observations
  int<lower=1> J;                       // number of genes
  array[N] int<lower=1, upper=J> j_idx; // gene index for each row
  vector[N] y;                          // scaled logFC
  vector[N] se_y;                       // scaled SE(logFC_ij)
}

parameters {
  real mu_delta;                 // global mean delta
  real<lower=1e-6> tau;          // SD of delta_j
  vector[J] z_delta;           // raw random effects
  real<lower=1e-6> sigma_pert;        // extra noise beyond SE
}

transformed parameters {
  vector[J] delta;
  delta = mu_delta + tau * z_delta;
}

model {
  // Priors
  mu_delta ~ normal(0, 0.1);
  tau      ~ exponential(1);
  sigma_pert    ~ exponential(1);

  z_delta ~ normal(0, 1);

  // Likelihood: measurement error + extra noise
  for (n in 1:N) {
    y[n] ~ normal(delta[ j_idx[n] ], sqrt(se_y[n]^2 + sigma_pert^2));
  }
}

generated quantities {
  vector[N] y_rep;
  for (n in 1:N)
    y_rep[n] = normal_rng(delta[j_idx[n]], sqrt(se_y[n]^2 + sigma_pert^2));
}
