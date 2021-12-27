//#include <TMB.hpp>
//#include "functions.hpp"
//#include "distributions.hpp"
#include "model.hpp"
#include "DerivedQuantities.hpp"


template<class Type>
Type objective_function<Type>::operator()(){
  logisticGrowth<Type>* inst = logisticGrowth<Type>::getinstance();
  
  DATA_VECTOR(y);
  
  PARAMETER_VECTOR(theta);
  PARAMETER(ln_sig);
  PARAMETER(ln_tau);
  PARAMETER_VECTOR(u);
  
  inst->y = y;
  inst->ln_sig = ln_sig;
  inst->ln_tau = ln_tau;
  inst->u = u;
  
  Type r = exp(theta[0]);
  Type K = exp(theta[1]);
  inst->r = exp(theta[0]);
  inst->K = exp(theta[1]);
  
  REPORT(r);
  REPORT(K);
  ADREPORT(r);
  ADREPORT(K);
  REPORT(u);
  
  Type nll = inst -> evaluate();
  
  Type sigma = exp(ln_sig);
  Type tau = exp(ln_tau);
  
  Type H = MSY(r,K);
  REPORT(H);
  ADREPORT(H);
  
  return nll;

}

