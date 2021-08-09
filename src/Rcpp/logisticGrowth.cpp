//#include <TMB.hpp>
//#include "functions.hpp"
//#include "distributions.hpp"
#include "model.hpp"


template<class Type>
Type objective_function<Type>::operator()(){
  logisticGrowth<Type>* inst = logisticGrowth<Type>::getinstance();
  
  DATA_VECTOR(y);
  
  PARAMETER_VECTOR(theta);
  PARAMETER(ln_sig);
  PARAMETER(ln_tau);
  PARAMETER_VECTOR(u);
  
  inst->y = y;
  inst->theta = theta;
  inst->ln_sig = ln_sig;
  inst->ln_tau = ln_tau;
  inst->u = u;
  
  return inst -> evaluate();
}

