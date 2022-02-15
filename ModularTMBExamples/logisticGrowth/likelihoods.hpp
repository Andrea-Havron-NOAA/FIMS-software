#ifndef LIKELIHOODS_HPP
#define LIKELIHOODS_HPP
#include <cmath>

/**
 * @brief lognormal density function
 * 
 * @tparam Type 
 * @param x observations
 * @param meanlog mean of the distribution on the log scale
 * @param sdlog standard deviate of the distribution on the log scale
 * @param give_log if true, return the log of the density function
 * @return lognormal density given observation and parameter values
 */
template<class Type>
Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  Type resid = (log(x)-meanlog)/sdlog;
  Type logres = -log(sqrt(2*M_PI)) - log(sdlog) - Type(0.5)*resid*resid - log(x);
  if(give_log) return logres; else return exp(logres);
  return logres;
}


#endif
