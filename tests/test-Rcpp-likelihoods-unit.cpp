#include "include/catch_amalgamated.hpp"
#include "../src/Rcpp/likelihoods.hpp"

TEST_CASE( "dlognorm works", "[likelihoods]" ) {
  
  SECTION( "dlognorm works with different input values" ) {
    
    REQUIRE( (dlognorm(1.0, 0.0, 1.0, true) - (-0.9189385)) <= 1e-4 ); // dlnorm(1.0, 0.0, 1.0, TRUE) = -0.9189385 in R
    REQUIRE( (dlognorm(1.0, 0.0, 1.0, false) - (0.3989423)) <= 1e-4 ); // dlnorm(1.0, 0.0, 1.0, FALSE) = 0.3989423 in R
    REQUIRE( (dlognorm(5.0, 10.0, 2.5) - (0.0001142879)) <= 1e-4 ); // dlnorm(5.0, 10.0, 2.5, FALSE) = 0.0001142879 in R
    REQUIRE( dlognorm(5.0, 10.0, 2.5, false) != 1.0 ); // dlnorm(5.0, 10.0, 2.5, FALSE) = 0.0001142879 in R
  }

  SECTION( "dlognorm give_log can be 1 or 0" ) {

    REQUIRE( (dlognorm(1.0, 0.0, 1.0, 1) - (-0.9189385)) <= 1e-4 ); // dlnorm(1.0, 0.0, 1.0, TRUE) = -0.9189385 in R
    REQUIRE( (dlognorm(1.0, 0.0, 1.0, 0) - (0.3989423)) <= 1e-4 ); // dlnorm(1.0, 0.0, 1.0, FALSE) = 0.3989423 in R

  }

}

