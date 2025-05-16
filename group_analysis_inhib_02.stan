data {
  int<lower=0> N;
  int<lower=0> N_grp;
  int<lower=0> N_spl;
  int<lower=0> grp[N];
  int<lower=0> spl[N];
  vector<lower=0>[N] Y;
}

transformed data {
  vector[N] Y_sm;
  real<lower=0> mu_Y;
  
  mu_Y = mean(Y);
  Y_sm = log(Y/mu_Y);
  
}

parameters {
  real a;
  vector[N_grp] a_grp;
  vector[N_spl] a_spl;
  real<lower=0> sigma;
  real<lower=0> sigma_spl;
}

model {
  vector[N] mu;
  
  a ~ normal(0,0.1);
  a_grp ~ normal(0,1);
  a_spl ~ normal(0,sigma_spl) ;
  sigma ~ exponential(5);
  sigma_spl ~ exponential(5);
  
  for (n in 1:N) {
    mu[n] =  a + a_grp[grp[n]] + a_spl[spl[n]];
    Y_sm[n] ~ normal(mu[n], sigma);
  }
}

generated quantities  {
  vector[N_grp] s_grp;
  vector[N] mu;
  vector[N] log_lik;      // Vector to store log-likelihoods
  
  for (n in 1:N) {
    mu[n] =  a + a_grp[grp[n]] + a_spl[spl[n]];
    log_lik[n] = normal_lpdf(Y_sm[n] | mu[n], sigma);  // Compute log-likelihood
  }
  
  s_grp = exp(a + a_grp) * mu_Y;
}
