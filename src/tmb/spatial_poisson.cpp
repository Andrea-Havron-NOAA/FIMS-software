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

  PARAMETER(b0); //intercept
  PARAMETER(ln_kappa); //spatial correlation decay
  PARAMETER(ln_tau); //spatial precision 
  PARAMETER_VECTOR(omega); //spatial random effect
  
  int i; 
  int n = y.size();
  //int nv = omega.size()
  
  Type kappa = exp(ln_kappa);
  Type tau = exp(ln_tau);
  Type marg_sp_sd = 1/(2*sqrt(M_PI)*kappa*tau);
  Type Range = sqrt(8)/kappa;

  Type nll = 0.0;

  //Spatial Likelihood
  //Define precision matrix
  //if(sparse == 0) Matrix<Type> Q = pow(kappa,4) * M0 + 2 * pow(kappa,2) * M1 + M2;
  SparseMatrix<Type> Q = pow(kappa,4) * M0 + 2 * pow(kappa,2) * M1 + M2;
  nll += SCALE(GMRF(Q), 1/tau)(omega);  
  SIMULATE{
    omega = GMRF(Q).simulate()/tau;
    REPORT(omega);
  }
  
  //Data Likelihood
  
  vector<Type> eta(n);
  
  for(i=0; i<n; i++){ 
    eta(i) = b0 + omega(v_i(i),0);
    nll -= dpois(y(i), exp(eta(i)), true);
  } 

  SIMULATE{
    y = rpois(exp(eta(i)));
    REPORT(y);
  }

  REPORT(Range);
  REPORT(marg_sp_sd);
  REPORT(omega);
  REPORT(nll);

  return nll;

}
