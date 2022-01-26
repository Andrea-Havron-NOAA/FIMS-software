// Separable covariance on lattice with AR1 structure in each direction.
// Modified from https://github.com/kaskr/adcomp/blob/master/TMB/inst/examples/ar1xar1.cpp
#include "model.hpp"

/* Parameter transform */
template <class Type>
Type f(Type x){return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);}

template<class Type>
Type objective_function<Type>::operator() ()
{

  ar1xar1<Type>* inst = ar1xar1<Type>::getinstance();
  
  DATA_VECTOR(y); 
  DATA_VECTOR_INDICATOR(keep,y);
  PARAMETER_ARRAY(eta);
  PARAMETER(transf_phi1); /* fastest running dim */
  PARAMETER(transf_phi2); /* slowest running dim */
  Type phi1 = f(transf_phi1);
  Type phi2 = f(transf_phi2);

  inst->y = y;
  inst->keep = keep;
  inst->phi1 = phi1;
  inst->phi2 = phi2;
  inst->eta = eta;
  
  Type nll = inst -> evaluate();

  SIMULATE{
    SEPARABLE(AR1(phi2), AR1(phi1)).simulate(eta);
    vector<Type> ln_lambda = eta;//convert array to vector for rpois simulation
    y = rpois(exp(ln_lambda));//rpois cannot accept array 
    REPORT(eta);
    REPORT(y);
  }
  
  REPORT(phi1);
  REPORT(phi2);
  REPORT(eta);
  ADREPORT(phi1);
  ADREPORT(phi2);

  return nll;
}
