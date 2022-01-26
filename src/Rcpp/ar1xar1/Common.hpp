#ifndef COMMON_HPP
#define COMMON_HPP

#define TMB_MODEL

#ifdef TMB_MODEL

#include <TMB.hpp>
using namespace tmbutils;
using namespace density;
//#include "tmbutils_extra.hpp"
//#include <Rcpp.h>
//using namespace Rcpp;

template<typename Type>
struct model_traits{
 typedef typename CppAD::vector<Type> data_vector;
 typedef data_indicator<tmbutils::vector<Type> , Type>  data_indicator;
 typedef typename tmbutils::array<Type> array;
};

template<typename T>
T exp(const T& x){
  return exp(x);
}

template <class T>
const T log(const T& x){return std::log(x);}

#endif

#ifdef STD_LIB

#include <cmath>
#include <vector>



template<typename T>
T exp(const T& x){
  return std::exp(x);
}

template <class T>
const T log(const T& x){return std::log(x);}

// not implemented yet
// 
// template<class distribution1, class distribution2>
// inline density::SEPARABLE_t<distribution1,distribution2> SEPARABLE_(distribution1 a, distribution2 b){
//   return density::SEPARABLE(a,b);
// }
// 
// template<class scalartype, class distribution>
// // see https://github.com/kaskr/adcomp/blob/cf781e381b17d911e9cfbbfcdbb8748e966530fc/TMB/inst/include/tmbutils/density.hpp#L453
// inline density::AR1_t<distribution> AR1_(scalartype phi_, distribution f_){
//   return density::AR1(phi_, f_);
// }
// template <class scalartype>
// inline density::AR1_t<N01<scalartype> > AR1_(scalartype phi_){
//   return density::AR1(phi_);
// }

 template<class Type>
 Type squeeze(Type u){
   Type eps = std::numeric_limits<double>::epsilon();
   u = (1.0 - eps) * (u - .5) + .5;
   return u;
 }

#endif

#endif
