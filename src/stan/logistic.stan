// Logistic growth model 
// Andrea Havron, 2021

data{
  int<lower=0> N;
  vector[N] y;
  real hyp_sig;
  real hyp_tau;
  vector[2] hyp_theta1;
  vector[2] hyp_theta2;
  int prior_type;
}

parameters{
  vector<lower=log(0.01),upper=log(1000)>[2] theta;
  real ln_sig;
  real ln_tau;
  vector<lower=0.001>[N] u;
  //real initu;
}

transformed parameters {
  real r;
  real K;
  real sigma;
  real tau;
  vector[N] umed;
  
  r = exp(theta[1]);
  K = exp(theta[2]);
  sigma = exp(ln_sig);
  tau = exp(ln_tau);
  
  //umed[1] = initu;
  // (conditional) prior distribution of Ps (from state equations):
  for (t in 2:N) { 
    umed[t] =log(u[t-1] + r*u[t-1]*(1-u[t-1]/K));
  }
  
}

model {
  //prior_type = 0:
  if(prior_type == 1){
    // prior distribution of r:
    target += lognormal_lpdf(r | hyp_theta1[1],hyp_theta1[2]);
  
    //prior distribution of K:
    target += lognormal_lpdf(K | hyp_theta2[1],hyp_theta2[2]);
   
    // prior distribution of sigma:
    target += exponential_lpdf(sigma | hyp_sig);
  
    // prior distribution of tau: 
    target += exponential_lpdf(tau | hyp_tau);
  }
  //state likelihoods
  //target += lognormal_lpdf(u[1] | umed[1],sigma);
  
  for(t in 2:N){
     target += lognormal_lpdf(u[t] | umed[t],sigma); 
  }
  
  // sampling distribution:
  for (t in 1:N) { 
    target += lognormal_lpdf(y[t] | log(u[t]), tau);
  }
}


