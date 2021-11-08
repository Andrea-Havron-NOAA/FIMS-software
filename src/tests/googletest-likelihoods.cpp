#include "gtest/gtest.h"
#include "Rcpp/likelihoods.hpp"

namespace {

  TEST(dlognormTest, DoubleInput) {
    EXPECT_NEAR(dlognorm<double>(1.0, 0.0, 1.0, true), -0.9189385, 0.0001); // dlnorm(1.0, 0.0, 1.0, TRUE) = -0.9189385 in R
  }

}

