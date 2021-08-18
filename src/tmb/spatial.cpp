//spatial delta-Gamma model
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

  PARAMETER_VECTOR(b0); //intercept
  PARAMETER_VECTOR(ln_kappa); //spatial correlation decay
  PARAMETER(ln_tau_z); //spatial precision for P/A
  PARAMETER(ln_tau_y); //spatial precision for Abundance
  PARAMETER_MATRIX(omega); //spatial random effect
  
  int i;j; 
  int ni = y.size();
  int nj = omega.rows(0).size()
  
  vector<Type> kappa = exp(ln_kappa);
  vector<Type> tau(nj);
  if(CppAD::Variable(ln_tau_z)) tau(0) = exp(ln_tau_z);
  if(CppAD::Variable(ln_tau_y)) tau(nj) = exp(ln_tau_y);
  vector<Type> marg_sp_sd = 1/(2*sqrt(M_PI)*kappa*tau);
  vector<Type> Range = sqrt(8)/kappa;

  Type nll = 0.0;

  //Spatial Likelihood
  //Define precision matrix
  for(j=0; j<n_j; j++){
    if(CppAD::Variable(ln_kappa(j))){ 
      SparseMatrix<Type> Q = pow(kappa(j),4) * M0 + 2 * pow(kappa(j),2) * M1 + M2;
    
      //Likelihood of spatial random effects
      nll += GMRF(Q)(omega.col(j));
      //scale omega using precision parameter
      omega.col(i) = omega.col(j)/tau(j);
      SIMULATE{
        omega.col(i) = GMRF(Q).simulate()/tau(i);
      }
    }
  }
  
  vector<Type> eta_z(n);
  vector<Type> eta_y(n);
  vector<Type> pi(n);
  vector<Type> mu(n);
  //Data Likelihood
  for(i=0; i<n; i++){ 
    
    eta_z(i) = b0(0);
    if(CppAD::Variable(tau_z)) eta_z(i) += b0(0) + omega(v_i(i),0);
    eta_y(i) = b0(nj);
    if(CppAD::Variable(tau_z)) eta_y(i) += b0(nj) + omega(v_i(i),nj);
    pi(i) = 1 / (1 + exp(-eta_z(i)));
    mu(i) = exp(eta_y(i));
    
    if(y(i)==0){
      nll -= 
    }
    nll -= dbinom(y(i), n, pi(i), true);
    SIMULATE{
      y = rbinom(n, pi(i));
    }
  } 

  SIMULATE{
    REPORT(y);
    REPORT(omega);
  }

  REPORT(Range);
  REPORT(marg_sp_sd);
  REPORT(pi);
  REPORT(omega);
  REPORT(nll);

  return nll;

}
