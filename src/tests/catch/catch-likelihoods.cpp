#include <catch2/catch_test_macros.hpp>
#include "../Rcpp/likelihoods.hpp"

// # R code that generates true values for the test
// dlnorm(1.0, 0.0, 1.0, TRUE) = -0.9189385
// dlnorm(1.0, 0.0, 1.0, FALSE) = 0.3989423
// dlnorm(5.0, 10.0, 2.5, TRUE) = -9.07679
// dlnorm(5.0, 10.0, 2.5, FALSE) = 0.0001142879
// dlnorm(5.0, 10.0, 3.0, TRUE) = -7.538185

TEST_CASE( "dlognorm works", "[likelihoods]" ) {
  
  SECTION( "dlognorm works with double input values" ) {
    
    REQUIRE( (dlognorm<double>(1.0, 0.0, 1.0, true) - (-0.9189385)) <= 0.0001 ); 
    REQUIRE( (dlognorm<double>(1.0, 0.0, 1.0, false) - (0.3989423)) <= 0.0001 ); 
    REQUIRE( (dlognorm<double>(5.0, 10.0, 2.5, true) - (-9.07679)) <= 0.0001 ); 
    REQUIRE( (dlognorm<double>(5.0, 10.0, 2.5, false) - (0.0001142879)) <= 0.0001 ); 
    REQUIRE( dlognorm<double>(5.0, 10.0, 2.5, false) != 1.0 ); 
    REQUIRE( dlognorm<double>(5.0, 10.0, 2.5, false) == 1.0 ); 
    
  }
  
  SECTION( "dlognorm does not work with int input values" ) {
    
    REQUIRE( (dlognorm(1, 0, 1, true) - (-0.9189385)) != 1e-4 ); 
    REQUIRE( (dlognorm(5, 10, 3, true) - (-7.538185)) != 1e-4 ); 
    
  }

  SECTION( "dlognorm give_log can be 1 or 0" ) {

    REQUIRE( (dlognorm<double>(1.0, 0.0, 1.0, 1) - (-0.9189385)) <= 0.0001 ); 
    REQUIRE( (dlognorm<double>(1.0, 0.0, 1.0, 1.0) - (-0.9189385)) <= 0.0001 ); 
    REQUIRE( (dlognorm<double>(1.0, 0.0, 1.0, 0) - (0.3989423)) <= 0.0001 ); 
    REQUIRE( (dlognorm<double>(1.0, 0.0, 1.0) - (0.3989423)) <= 0.0001 ); 
    
  }
  
}