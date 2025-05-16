data {
  int<lower=0> N;      //total number of measurements
  int<lower=0> N_run;  // numer of runs
  int<lower=0> N_spl;  // numebr of samples
  int<lower=0> run[N]; 
  int<lower=0> spl[N];
  vector[N] Y;         // outcome variable
}

transformed data {
 vector[N] Y_sm;
 real mean_Y = mean(Y);
 
 // scale outcome var on global mean
 Y_sm = Y/mean_Y;

 }

parameters {
  real<lower=0> a;
  vector[N_run] z_run;
  vector[N_spl] z_spl;
  real<lower=0> sigma_run;
  real<lower=0> sigma_spl;
  real<lower=0> sigma;
}

model {
  vector[N] mu;
  
  a     ~ normal(1, 0.1);
  z_run ~ normal(0,1);
  z_spl ~ normal(0,1);
  sigma_run ~ exponential(5);  
  sigma_spl ~ exponential(5);  
  sigma ~ exponential(5);
  
  // likelihood. run and sample are random effects
  for (n in 1:N) {
    mu[n] = a + z_run[run[n]] * sigma_run + z_spl[spl[n]] * sigma_spl;
    Y_sm[n] ~ normal(mu[n], sigma) T[0,];
  }
}

generated quantities {
  real<lower = 0> sigma_pooled_run;
  real<lower = 0> sigma_pooled_spl;
  
  // calculate pooled standardevs
  sigma_pooled_run = sqrt(square(sigma_run) + square(sigma_spl) + square(sigma));
  sigma_pooled_spl = sqrt(square(sigma_spl) + square(sigma));
  
}
