data {
  int<lower=0> N;
  int<lower=0> N_grp;
  int<lower=0> N_spl;
  // int<lower=0> N_run;
  int<lower=0> grp[N];
  int<lower=0> spl[N];
  // int<lower=0> run[N];
  vector[N] Y_sm;
}

parameters {
  // real a;
  vector[N_grp] a_grp;
  // vector[N_run] z_run;
  vector[N_spl] z_spl;
  real<lower=0> sigma;
  // real<lower=0> sigma_run;
  real<lower=0> sigma_spl;
}

model {
  vector[N] mu;

  a_grp ~ normal(1,0.1);
  // z_run ~ normal(0,1) ;
  z_spl ~ normal(0,1) ;
  sigma ~ exponential(1);
  // sigma_run ~ exponential(1);
  sigma_spl ~ exponential(1);
  
  for (n in 1:N) {
    mu[n] =  a_grp[grp[n]] + z_spl[spl[n]] * sigma_spl;// + z_run[run[n]] * sigma_run;
    Y_sm[n] ~ normal(mu[n], sigma) T[0,];
  }
}
// 
// generated quantities  {
//  vector[N_grp] s;
//  
//  s = a + a_grp;
// }
