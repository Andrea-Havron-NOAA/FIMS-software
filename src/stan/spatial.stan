//spatial model 
// Andrea Havron, 2021

data{
  int<lower=0> N;
  int<lower=0> NV;
  int y[N];
  int vi[N];
  matrix[NV,NV] M0;
  matrix[NV,NV] M1;
  matrix[NV,NV] M2;
//  vector[2] kap_tau_pr_mu;
//  vector[2] kap_tau_pr_var;
  int prior_type;
  real kappa;
  real tau;
}

parameters{
  real b0;
  //vector[2] ln_kap_tau;
  vector[NV] omega;
}

transformed parameters {
 // real kappa;
//  real tau;
 // real theta;
  matrix[NV,NV] Q;
 // matrix[2,2] Sigma;
  vector[N] eta;
  vector[NV] mu; 
  
 // kappa = sqrt(8)/50;//exp(ln_kap_tau[1]);
//  tau = 1/(sqrt(4*pi()*0.75)*kappa);//exp(ln_kap_tau[2]);
  //theta = -0.35;
  
  //precision matrix
  Q = pow(tau,2) * (pow(kappa,4) * M0 + 2 * pow(kappa,2) * M1 + M2);
  
  //linear mean
  for(i in 1:N){
    eta[i] = b0 + omega[vi[i]];
  }
  
  mu = rep_vector( 0.0, NV );
  
 // Sigma[1,1] = kap_tau_pr_var[1];
//  Sigma[2,2] = kap_tau_pr_var[2];
//  Sigma[1,2] = sqrt(kap_tau_pr_var[1])*sqrt(kap_tau_pr_var[2])*theta;
  //Sigma[2,1] = Sigma[1,2];
}

model {
  if(prior_type == 1){
     target += normal_lpdf(b0 | 0, 5);
     //target += cauchy_lpdf(theta | 0,2.5);
  //   target += multi_normal_lpdf(ln_kap_tau | kap_tau_pr_mu, Sigma);
     //target += normal_lpdf(ln_kappa | lnkapPr[1], lnkapPr[2]);
     //target += normal_lpdf(ln_tau | lntauPr[1], lntauPr[2]);
  }
  //spatial effects
  target += multi_normal_prec_lpdf(omega|mu,Q);
  
  //poisson obs
  for(i in 1:N){
    target += poisson_log_lpmf(y[i] | eta[i]);
  }
}


