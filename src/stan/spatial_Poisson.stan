// Logistic growth model 
// Andrea Havron, 2021

data{
  int<lower=0> N;
  int<lower=0> NV;
  int y[N];
  int vi[N];
  matrix[NV,NV] M0;
  matrix[NV,NV] M1;
  matrix[NV,NV] M2;
  vector[2] lnkapPr;
  vector[2] lntauPr;
  int prior_type;
}

parameters{
  real b0;
  real ln_kappa;
  real ln_tau;
  vector[NV] omega;
}

transformed parameters {
  real kappa;
  real tau;
  matrix[NV,NV] Q;
  vector[N] eta;
  vector[NV] mu; 
  
  kappa = exp(ln_kappa);
  tau = exp(ln_tau);
  
  //precision matrix
  Q = pow(kappa,4) * M0 + 2 * pow(kappa,2) * M1 + M2;
  
  //linear mean
  for(i in 1:N){
    eta[i] = b0 + omega[vi[i]];
  }
  
  mu = rep_vector( 0.0, NV );
}

model {
  if(prior_type == 1){
     target += normal_lpdf(b0 | 0, 5);
     target += normal_lpdf(ln_kappa | lnkapPr[1], lnkapPr[2]);
     target += normal_lpdf(ln_tau | lntauPr[1], lntauPr[2]);
  }
  //spatial effects
  target += multi_normal_prec_lpdf(omega|mu,Q);
  
  //poisson obs
  for(i in 1:N){
    target += poisson_log_lpmf(y[i] | eta[i]);
  }
}


