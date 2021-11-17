#include "gtest/gtest.h"
#include "../Rcpp/likelihoods.hpp"

// # R code that generates true values for the test
// # TestSuiteName: dlognormTest; TestName: DoubleInput, IntInput, givelog
// dlnorm(1.0, 0.0, 1.0, TRUE) = -0.9189385
// dlnorm(1.0, 0.0, 1.0, FALSE) = 0.3989423
// dlnorm(5.0, 10.0, 2.5, TRUE) = -9.07679
// dlnorm(5.0, 10.0, 2.5, FALSE) = 0.0001142879
// dlnorm(5.0, 10.0, 3.0, TRUE) = -7.538185


namespace {

  // dlognorm(x, meanlog, sdlog, give_log=0)
  // Test dlognorm with double input values
  TEST(dlognormTest, DoubleInput) {
    
    EXPECT_NEAR( dlognorm<double>(1.0, 0.0, 1.0, true) , -0.9189385 , 0.0001 ); 
    EXPECT_NEAR( dlognorm<double>(1.0, 0.0, 1.0, false) , 0.3989423 , 0.0001 ); 
    EXPECT_NEAR( dlognorm<double>(5.0, 10.0, 2.5, true) , -9.07679 , 0.0001 ); 
    EXPECT_NEAR( dlognorm<double>(5.0, 10.0, 2.5, false) , 0.0001142879 , 0.0001 ); 
    EXPECT_NE( dlognorm<double>(5.0, 10.0, 2.5, false) , 1.0 ); 
  
  }
  
  // Test dlognorm with integer input values
  TEST(dlognormTest, IntInput) {
    
    EXPECT_NE( dlognorm(1, 0, 1, true) , -0.9189385 );
    EXPECT_NE( dlognorm<double>(5, 10, 3, true) , -7.538185 );
    
  }
  
  // Test dlognorm with different values of give_log
  TEST(dlognormTest, give_log) {
    
    EXPECT_NEAR( dlognorm<double>(1.0, 0.0, 1.0, 1) , -0.9189385 , 0.0001 );
    EXPECT_NEAR( dlognorm<double>(1.0, 0.0, 1.0, 1.0) , -0.9189385 , 0.0001 );
    EXPECT_NEAR( dlognorm<double>(1.0, 0.0, 1.0, 0) , 0.3989423 , 0.0001 ); 
    EXPECT_NEAR( dlognorm<double>(1.0, 0.0, 1.0) , 0.3989423 , 0.0001 ); 
    
  }
  
}

