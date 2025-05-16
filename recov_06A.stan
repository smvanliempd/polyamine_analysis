data {
  int<lower=0> N; // Number of observations
  int<lower=0> N_grp; // Number of groups
  int<lower=0> N_spl;
  int<lower=0> grp[N]; // Group index for each observation (1 to N_grp)
  int<lower=0> spl[N];
  vector[N_grp] mu_grp;
  vector[N] Y; // The duplicate measurements for each observation
}

transformed data {
  vector[N] y_sm;
  vector[N_grp] mu_grp_sm;

  y_sm = Y/mu_grp[1];
  mu_grp_sm = (mu_grp/mu_grp[1]);
  
}

parameters {
  ordered[N_grp] alpha_grp;
  vector[N_spl] z_spl;
  real<lower=0> sigma;
  real<lower=0> sigma_spl;
}

model {
  vector[N] mu;
  
  // Priors
  alpha_grp  ~ normal(mu_grp_sm, 0.1); // Prior for group means, 0.01 <- too tight?
  z_spl ~ normal(0,1);
  sigma_spl ~ exponential(1);
  sigma ~ exponential(1); // Prior for measurement variability
  

  // Likelihood for the duplicate measurements (measurement variability)
  for (n in 1:N) {
    mu[n] = alpha_grp[grp[n]] + z_spl[spl[n]] * sigma_spl;
    y_sm[n] ~ normal(mu[n], sigma) T[0,]; 
  }
}

generated quantities {
  real<lower=0> tau;
  real<lower=0> tau_spl;
  real<lower=0> R_100;
  real<lower=0> R_1000;

  tau = sigma * mu_grp[1];
  tau_spl = sigma_spl * mu_grp[1];

  R_100 = (alpha_grp[2] - alpha_grp[1])/(alpha_grp[3] - alpha_grp[1]);
  R_1000 = (alpha_grp[4] - alpha_grp[1])/(alpha_grp[5] - alpha_grp[1]);
}
