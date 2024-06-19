functions {
  /*
  * Alternative to poisson_log_rng() that 
  * avoids potential numerical problems during warmup
  */
  int poisson_log_safe_rng(real eta) {
    real pois_rate = exp(eta);
    if (pois_rate >= exp(20.79))
      return -9;
    return poisson_rng(pois_rate);
  }
}
data {
  int<lower=1> N;               // number of observations
  int<lower=1> K;               // number of nations
  int<lower=1> J;               // number of regions
  int<lower=1> k[N];            // nation index for each observation
  int<lower=1> j[N];            // region index for each observation
  int<lower=0> deaths[N];       // outcome variable
  vector<lower=0>[N] expected;  // exposure variable
  vector[N] UVBI;               // covariate UVBI
}
parameters {
  real beta0;                   // intercept
  real beta1;                   // fixed effect of UVBI
  real<lower=0> sigma_s;        // std deviation for nation random effect
  real<lower=0> sigma_u;        // std deviation for region random effect
  real<lower=0> sigma_e;        // std deviation for county random effect
  vector[K] s;                  // random intercepts for nations
  vector[J] u;                  // random intercepts for regions
  vector[N] e;                  // random intercepts for county
}
model {
  // Priors
  beta0 ~ normal(0, 1);
  beta1 ~ normal(0, 1);
  sigma_s ~ normal(0, 1);
  sigma_u ~ normal(0, 1);
  sigma_e ~ normal(0, 1);
  s ~ normal(0, sigma_s);   // random term with the intercept at level 3
  u ~ normal(0, sigma_u);   // random term with the intercept at level 2
  e ~ normal(0, sigma_e);   // random term with the intercept at level 1
  
  // Linear predictor
  vector[N] eta;
    for (n in 1:N) {
      eta[n] = log(expected[n]) + beta0 + beta1 * UVBI[n] + s[k[n]] + u[j[n]] + e[n];
    }  
    
  // Poisson regression
  // poisson_log(eta) is more efficient and stable alternative to poisson(exp(eta))
  deaths ~ poisson_log(eta);
} 
generated quantities {
  int y_rep[N];     // Generated quantities for replicated data
  real log_lik[N];  // Pointwise log-likelihood
  
  for (n in 1:N) {
    real eta_n = log(expected[n]) + beta0 + beta1 * UVBI[n] + s[k[n]] + u[j[n]] + e[n];
    y_rep[n] = poisson_log_safe_rng(eta_n);
    // Calculate Poisson log-likelihood for each observation
    log_lik[n] = poisson_log_lpmf(deaths[n] | eta_n);
  }
}
