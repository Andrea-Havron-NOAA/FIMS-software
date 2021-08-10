#include <Rcpp.h>
#include "likelihoods.hpp"
using namespace Rcpp;

// [[Rcpp::export(rng=false)]]
SEXP dlognorm(SEXP x, SEXP meanlog, SEXP sdlog, SEXP give_log){
  return wrap(dlognorm<double>(as<double>(x), as<double>(meanlog), as<double>(sdlog), as<int>(give_log)));
}

/***R
y <- rlnorm(100)
r.res <- sum(dlnorm(y, log = TRUE))
rcpp.res <- rep(0,length(y))
for(i in 1:length(y)){
  rcpp.res[i] <- dlognorm(y[i], 0, 1, 1)
}
print(c(r.res, sum(rcpp.res)))
*/
