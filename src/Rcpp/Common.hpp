#ifndef COMMON_HPP
#define COMMON_HPP

#define TMB_MODEL

#ifdef TMB_MODEL

#include <TMB.hpp>
#include <Rcpp.h>
using namespace Rcpp;

template<typename Type>
struct model_traits{
 typedef typename CppAD::vector<Type> data_vector;
 typedef typename CppAD::vector<Type> parameter_vector;
};

template<typename T>
T exp(const T& x){
  return exp(x);
}


#endif


#endif