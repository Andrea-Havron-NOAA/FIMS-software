// Shaefer surplus production model 
// adapted from Millar and Renate, 1999 
// Andrea Havron, 2021

data{
  int<lower=0> N;
  vector[N] C;
  vector[N] I;
}

parameters{
  real<lower=0.01, upper=1.2> r;
  real<lower=10, upper=1000> K;
  real<lower=0.5, upper=100> iq;
  real isigma2;
  real itau2;
  vector<lower=0.001, upper=2.0>[N] P;
  real initPmed;
}

transformed parameters {
  real q;
  real sigma;
  real tau;
  vector[N] Imed;
  vector<lower=0.001, upper=2.0>[N] Pmed;
  
  q = 1/iq;
  sigma = sqrt(1/isigma2);
  tau = sqrt(1/itau2);
  Pmed[1] = initPmed; 
  
  // (conditional) prior distribution of Ps (from state equations):
  for (t in 2:N) { 
    Pmed[t] = P[t-1] + r*P[t-1]*(1-P[t-1]) - C[t-1]/K;
  }
  
  // sampling distribution:
  for (t in 1:N) { 
    Imed[t] = log(q*K*P[t]);
  }

}

model {
  // prior distribution of r: lognormal with 10% and 90% quantile at 0.13 and 0.48
  target += lognormal_lpdf(r | -1.38,3.845);

  //prior distribution of K: lognormal with 10% and 90% quantile at 80 and 300
  target += lognormal_lpdf(K | 5.042905,3.7603664);
 
  // prior distribution of q: instead of improper (prop. to 1/q) use just proper IG
  target += gamma_lpdf(iq | 0.001,0.001);
  
  // prior distribution of sigma2: inv. gamma with 10% and 90% qu. at 0.04 and 0.08
  target += gamma_lpdf(isigma2 | 3.785518,0.010223);

  // prior distribution of tau2: inv. gamma with 10% and 90% qu. at 0.05 and 0.15
  target += gamma_lpdf(itau2 | 1.708603,0.008613854);

  //state likelihoods
  target += lognormal_lpdf(P[1] | Pmed[1],sigma);
  
  for(t in 2:N){
     target += lognormal_lpdf(P[t] | log(Pmed[t]),sigma); 
  }
  
  // sampling distribution:
  for (t in 1:N) { 
    target += lognormal_lpdf(I[t] | Imed[t], tau);
  }
}

generated quantities {
 real MSP;
 real EMSP;
 real P1990;
 real B1990;
 
  // further management parameters and predictions:
  MSP = r*K/4;
  EMSP = r/(2*q);
  P1990 = P[N]+r*P[N]*(1-P[N]) -C[N]/K;
  B1990 = P1990*K;
}

