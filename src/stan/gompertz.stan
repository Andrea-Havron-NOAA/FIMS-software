// Gompertz state-space model 
// Andrea Havron, 2021

data{
  int<lower=0> N;
  vector[N] y;
  vector[2] hyp_sig;
  vector[2] hyp_tau;
  real hyp_theta1;
  real hyp_theta2;
  int prior_type;
}

parameters{
  vector[2] theta;
  real ln_sig;
  real ln_tau;
  vector[N] u;
}

transformed parameters {
  real sigma;
  real tau;
  vector[N] umed;
  
  sigma = exp(ln_sig);
  tau = exp(ln_tau);
  
  
  // (conditional) prior distribution of Ps (from state equations):
  for (t in 2:N) { 
    umed[t] = theta[1] + theta[2]*u[t-1];
  }
  
}

model {
  //prior_type = 0:improper prior
  if(prior_type == 1){
    // prior distribution of theta: 
    for(i in 1:2){
      target += normal_lpdf(theta[i] | hyp_theta1,hyp_theta2);
    }
    
    // prior distribution of sigma2:
    target += inv_gamma_lpdf(sigma | hyp_sig[1],hyp_sig[2]);
  
    // prior distribution of tau2: 
    target += inv_gamma_lpdf(tau | hyp_tau[1],hyp_tau[2]);
  }
  //state likelihoods
  //target += normal_lpdf(u[1]| umed[1],sigma);
  for(t in 2:N){
     target += normal_lpdf(u[t] | umed[t],sigma); 
  }
  
  // sampling distribution:
  for (t in 1:N) { 
    target += normal_lpdf(y[t] | u[t], tau);
  }
}


