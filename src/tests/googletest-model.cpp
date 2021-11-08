#include "gtest/gtest.h"
#include "Rcpp/model.hpp"

// # R code that generates true values for the test
// u <- c(5000, 6000, 7000, 5500)
// r <- 0.5
// K <- 10000
// n <- length(u) # 4
// eta <- c()
// for (t in 2:n) {
//   eta[t] <- u[t - 1] + r * u[t - 1] * (1 - u[t - 1] / K) # NA 6250 7200 8050
// }
// 
// y <- c(5050, 5950, 6950, 6000)
// sigma <- 0.1
// tau <- 0.03
// nll <- 0
// 
// for (t in 2:n) {
//   nll <- nll - dlnorm(u[t], log(eta[t]), sigma, TRUE)
// }
// 
// for (t in 1:n) {
//   nll <- nll - dlnorm(y[t], log(u[t]), tau, TRUE) # 58.13528
// } 
// 
// theta <- c(log(r), log(K)) #-0.6931472  9.2103404
// ln_sig <- log(sigma) # -2.302585
// ln_tau <- log(tau) # -3.506558

namespace {

  TEST(modelTest, eta) {
    EXPECT_NEAR(-0.9189383, -0.9189385, 0.0001); // dlnorm(1.0, 0.0, 1.0, TRUE) = -0.9189385 in R
  }
  
}


