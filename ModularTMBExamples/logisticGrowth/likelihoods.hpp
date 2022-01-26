//#include <Rcpp.h>
//#include <vector>
//#include <TMB.hpp>
//using namespace Rcpp;
//using namespace std;
#include <cmath>

template<class Type>
Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  Type resid = (log(x)-meanlog)/sdlog;
  Type logres = -log(sqrt(2*M_PI)) - log(sdlog) - Type(0.5)*resid*resid - log(x);
  if(give_log) return logres; else return exp(logres);
  return logres;
}

/* Not Implemented
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR of


template<class Type>
CppAD::vector<Type> sim_rlnorm(Type meanlog, Type sdlog, objective_function<Type>* of){
  Type y;
  SIMULATE {
    y = exp(rnorm(meanlog, sdlog));  // Simulate response
  }

  return y;
}
#undef  TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
*/
