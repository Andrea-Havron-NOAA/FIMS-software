//#include <Rcpp.h>
//#include <vector>
//#include <TMB.hpp>
//using namespace Rcpp;
//using namespace std;

template<class Type>
Type dlnorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  Type resid = (x-meanlog)/sdlog;
  Type logres = -log(sqrt(2*M_PI)) - log(sdlog) - Type(0.5)*resid*resid - log(x);
  //Type logres = -log(sdlog) - log(x) - 0.5*(pow((x-meanlog)/sdlog),2)
  if(give_log) return logres; else return exp(logres);
}

