// caluclation of signal attenuation
data {
  int<lower=0> N;
  int<lower=0> N_grp;
  int<lower=0> N_run;
  int<lower=0> grp[N];
  int<lower=0> run[N];
  vector[N] Y;
}

transformed data {
 vector[N] Y_sm;
 real mean_Y = mean(Y);

 Y_sm = Y/mean_Y;

 }

parameters {
  vector<lower=0>[N_grp] a_grp;
  vector[N_run] z_run;
  real<lower=0> sigma;
  real<lower=0> sigma_run;
}

model {
  vector[N] mu;
  
  a_grp ~ normal(0,0.5);
  z_run ~ normal(0,1);
  sigma ~ exponential(5);
  sigma_run ~ exponential(5);

  
  for (n in 1:N) {
    mu[n] = a_grp[grp[n]]+  z_run[run[n]] * sigma_run;
    Y_sm[n] ~ normal(mu[n], sigma) T[0,];
  }
}

generated quantities {
  real S_1000;
  vector[N] mu;
  vector[N] log_lik;      // Vector to store log-likelihoods

  S_1000  = (a_grp[2] - a_grp[1])/(a_grp[3]);
  
  for (n in 1:N) {
   mu[n] = a_grp[grp[n]]+  z_run[run[n]] * sigma_run;
    log_lik[n] = normal_lpdf(Y_sm[n] | mu[n], sigma);  // Compute log-likelihood
  }
}
