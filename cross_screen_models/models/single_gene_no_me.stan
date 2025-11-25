data {
  int<lower=1> N;                // number of observations (rows)
  int<lower=1> J;                // number of target genes (effects)
  array[N] int<lower=1, upper=J> j_idx; // gene index for each observation
  vector[N] y;                   // scaled_logFC for each observation
}

parameters {
  real mu_delta;                 // global mean effect
  real<lower=1e-6> tau;             // between-gene SD of delta_j
  vector[J] z_delta;           // non-centered random effects
  real<lower=1e-6> sigma;           // residual SD
}

transformed parameters {
  vector[J] delta;               // actual delta_j

  delta = mu_delta + tau * z_delta;
}

model {
  // Priors
  mu_delta ~ normal(0, 0.1);     // mean of the distribution of downstream effects
  // mu_delta small; most genes unaffected by perturb i
  tau      ~ exponential(1);     // SD of the distribution of downstream effects
  // small tau -> i doesn’t affect many genes; big tau -> i affects many genes
  sigma    ~ exponential(1);     // variation for the same gene across screens

  z_delta ~ normal(0, 1);      // normalised delta

  // Likelihood
  y ~ normal(delta[j_idx], sigma);
}

generated quantities {
  // For posterior predictive checks
  vector[N] y_rep;

  for (n in 1:N) {
    y_rep[n] = normal_rng(delta[j_idx[n]], sigma);
  }
}
