// #include <TMB.hpp>
#include "model.hpp"


template<class Type>
Type objective_function<Type>::operator() ()
{
  ThetaLog<Type>* inst = ThetaLog<Type>::get_instance();
  
  /* Data section */
  DATA_VECTOR(Y);
  /* Parameter section */
  PARAMETER_VECTOR(X);
  PARAMETER(logr0);
  PARAMETER(logtheta);
  PARAMETER(logK);
  PARAMETER(logQ);
  PARAMETER(logR);
  
  inst->Y = Y;
  inst->X = X;
  inst->logr0 =logr0;
  inst->logtheta = logtheta;
  inst->logK = logK;
  inst->logQ = logQ;
  inst->logR = logR;
  
  
  return inst->evaluate();
}

