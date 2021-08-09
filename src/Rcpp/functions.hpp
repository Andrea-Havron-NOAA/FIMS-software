#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

template <typename T>
  std::vector<T> logisticGrowth(const std::vector<T>& x, T r, T K)
{
  int n = x.size();
  std::vector<T> res(n);
  
  for(int j=0; j<n; ++j){
    res[j] =log(x[j] + r*x[j]*(1-x[j]/K));
  }
  
  return res;
}



// [[Rcpp::export]]
SEXP logisticGrowth(SEXP x, SEXP r, SEXP K) {
  
  switch (TYPEOF(x)){
    case REALSXP:{
      return wrap(logisticGrowth<double>(as<std::vector<double> >(x),as<double>(r),as<double>(K)));
    }
    case INTSXP:{
      return wrap(logisticGrowth<double>(as<std::vector<double> >(x),as<double>(r), as<double>(K)));
    }
    default: {
      warning(
        "Invalid SEXPTYPE %d (%s).\n",
        TYPEOF(x), type2name(x)
      );
      return R_NilValue;
    }
  }
  
}
