#include <TMB.hpp>

// dlnorm
template<class Type>
Type dlnorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
}

template<class Type>
Type objective_function<Type>::operator()()
{
  DATA_VECTOR(y);
  DATA_INTEGER(mod);
  
  PARAMETER_VECTOR(theta);
  PARAMETER(ln_sig);  Type sig = exp(ln_sig);
  PARAMETER(ln_tau);  Type tau = exp(ln_tau);
  
  PARAMETER_VECTOR(u);
  
  int n = y.size();
  int t;
  Type nll = 0;
  
  vector<Type> eta(n);
  if(mod == 0){
    for(t=1; t<n; t++){
      eta(t) = theta[0] + theta[1] * u[t-1];
      nll -= dnorm(u(t), eta(t), sig, true);
    }
    nll -= sum(dnorm(y,u,tau,true));
  }
  if(mod == 1){
    Type r = exp(theta[0]);
    Type K = exp(theta[1]);
    for(t=1; t<n; t++){
      eta(t) = log(u(t-1) + r*u(t-1)*(1-u(t-1)/K));
      nll -= dlnorm(u(t), eta(t), sig, true);
    }
    for(t=0; t<n; t++){
      nll -= dlnorm(y(t),log(u(t)),tau,true);
    }
  }
  return(nll);
}
