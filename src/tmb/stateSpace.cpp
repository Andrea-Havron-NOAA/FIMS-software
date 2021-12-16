#include <TMB.hpp>

// dlnorm
template<class Type>
Type dlnorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
}

//dinvgamma
template<class Type>
Type dinvgamma(Type x, Type shape, Type scale, int give_log=0){
  Type logres = -lgamma(shape) - shape*log(scale) - (shape+1)*x - 1/(scale*x);
  if(give_log) return logres; else return exp(logres);
}

template<class Type>
Type objective_function<Type>::operator()()
{
  DATA_VECTOR(y);
  DATA_INTEGER(mod);
  DATA_VECTOR(hyperpars); //theta1 mu, theta1 sd, theta2 mu, theta2 sd, sigma, tau

  PARAMETER_VECTOR(theta);
  PARAMETER(ln_sig);  Type sigma = exp(ln_sig);
  PARAMETER(ln_tau);  Type tau = exp(ln_tau);

  PARAMETER_VECTOR(u);

  int n = y.size();
  int t;
  Type nll = 0;


  vector<Type> eta(n);
  if(mod == 0){
    for(t=1; t<n; t++){
      eta(t) = theta[0] + theta[1] * u[t-1];
      nll -= dnorm(u(t), eta(t), sigma, true);
    }
    nll -= sum(dnorm(y,u,tau,true));
  }
  if(mod == 1){
    Type r = exp(theta[0]);
    Type K = exp(theta[1]);


    if(hyperpars.size() > 1){
      nll -= dlnrom(r, hyperpars_theta(0,0), hyperpars_theta(0,1), true);
      nll -= dlnrom(K, hyperpars_theta(1,0), hyperpars_theta(1,1), true);
      nll -= dexp(sigma, hyperpars_sd(0), true);
      nll -= dexp(tau, hyperpars_sd(1), true);
    }

    for(t=1; t<n; t++){
      eta(t) = log(u(t-1) + r*u(t-1)*(1-u(t-1)/K));
      nll -= dlnorm(u(t), eta(t), sigma, true);
    }
    for(t=0; t<n; t++){
      nll -= dlnorm(y(t),log(u(t)),tau,true);
    }
    REPORT(r);
    REPORT(K);
    ADREPORT(r);
    ADREPORT(K);
  }
  REPORT(sigma);
  REPORT(tau);
  ADREPORT(sigma);
  ADREPORT(tau);
  
  return(nll);
}
