//spatial Poisson model
#include <TMB.hpp>


template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density; //namespace with GMRF function
  using namespace Eigen; //namespace with SparseMatrix declaration

  DATA_VECTOR(y);
  DATA_IVECTOR(v_i); //location index
  DATA_SPARSE_MATRIX(M0); //sparse distance matrices
  DATA_SPARSE_MATRIX(M1); //sparse distance matrices
  DATA_SPARSE_MATRIX(M2); //sparse distance matrices
  //DATA_MATRIX(hyperpars); // r1:ln_kappa, r2:ln_tau, c1:mean, c2:sd
  DATA_VECTOR(kap_tau_pr_mu);
  DATA_VECTOR(kap_tau_pr_var);
  //DATA_VECTOR(kap_tau_pr_theta);
  DATA_INTEGER(prior_type);

  PARAMETER(b0); //intercept
  //PARAMETER(ln_kappa); //spatial correlation decay
  //PARAMETER(ln_tau); //spatial precision
  //PARAMETER_VECTOR(ln_kap_tau);
  //PARAMETER(theta);
 // PARAMETER(ln_phi);
 // PARAMETER(ln_spvar);
  PARAMETER_VECTOR(omega); //spatial random effect

  int i;
  int n = y.size();
  //int nv = omega.size()

  Type nll = 0.0;
  //Type theta = -0.35;
  vector<Type> resid(2);

  if(prior_type == 1){
    //nll -= dcauchy(theta, 0,2.5,true);
    //matrix<Type>Sigma(2,2);
    //Sigma(0,0) = kap_tau_pr_var(0);
    //Sigma(1,1) = kap_tau_pr_var(1);
    //Sigma(0,1) = theta * sqrt(kap_tau_pr_var(0)) * sqrt(kap_tau_pr_var(1));
    //Sigma(1,0) = Sigma(0,1);
    //resid = ln_kap_tau - kap_tau_pr_mu;
    //nll += MVNORM(Sigma)(resid);
    nll -= dnorm(b0, Type(0), Type(5), true);
   // nll -= dnorm(ln_phi, kap_tau_pr_mu(0), sqrt(kap_tau_pr_var(0)));
  //  nll -= dnorm(ln_spvar, kap_tau_pr_mu(1), sqrt(kap_tau_pr_var(1)));
    //nll -= dnorm(ln_kappa, hyperpars(0,0), hyperpars(0,1), true);
    //nll -= dnorm(ln_tau, hyperpars(1,0), hyperpars(1,1), true);
  }

  //Type kappa = exp(ln_kap_tau(0));
  //Type tau = exp(ln_kap_tau(1));
  //Type marg_sp_sd = 1/(2*sqrt(M_PI)*kappa*tau);
  //Type Range = sqrt(8)/kappa;
  Type marg_sp_sd = sqrt(0.75);//sqrt(exp(ln_spvar));
  Type Range = 50;//exp(ln_phi);
  Type kappa = sqrt(8)/Range;
  Type tau =  1/(2*sqrt(M_PI)*kappa*marg_sp_sd);

  //Spatial Likelihood
  //Define precision matrix
  SparseMatrix<Type> Q = pow(kappa,4) * M0 + 2 * pow(kappa,2) * M1 + M2;
  nll += SCALE(GMRF(Q), 1/tau)(omega);
  SIMULATE{
    omega = GMRF(Q).simulate()/tau;
    REPORT(omega);
  }

  //Data Likelihood

  vector<Type> eta(n);

  for(i=0; i<n; i++){
    eta(i) = b0 + omega(v_i(i));
    nll -= dpois(y(i), exp(eta(i)), true);
  }

  SIMULATE{
    y = rpois(exp(eta(i)));
    REPORT(y);
  }

  REPORT(Range);
  ADREPORT(Range);
  REPORT(marg_sp_sd);
  ADREPORT(marg_sp_sd);
  REPORT(omega);
  REPORT(nll);

  return nll;

}
