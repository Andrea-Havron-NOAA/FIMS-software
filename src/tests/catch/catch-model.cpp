#include <catch2/catch_test_macros.hpp>
#include "../Rcpp/model.hpp"

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

TEST_CASE( "Logistic growth model eta test", "[logistic-growth]" ) {
  
  //logisticGrowth<Type>* inst = logisticGrowth<Type>::getinstance();
  logisticGrowth<double>* inst = new logisticGrowth<double>(); //non-singleton usage
  
  model_traits<double>::parameter_vector u={5000, 6000, 7000, 5500};
  double r = 0.5;
  double K = 10000;
  
  auto eta_val = inst->calculateEta(u, r, K);
  
  model_traits<double>::parameter_vector true_eta={0, 6250, 7200, 8050};
  
  REQUIRE((eta_val == true_eta));

  SECTION( "Empty u test" ) {
  
    u = {};
    r = 0.5;
    K = 1000;
    
    eta_val = inst->calculateEta(u, r, K);
    true_eta = {};
    
    REQUIRE( eta_val == true_eta );
  }
  
  SECTION( "Different r test" ) {
    
    u = {5000, 6000, 7000, 5500};
    r = 1.5;
    eta_val = inst->calculateEta(u, r, K);
    
    true_eta = {0, 6250, 7200, 8050};
    
    REQUIRE( eta_val != true_eta );
  }
}

TEST_CASE( "Logistic growth model nll test", "[logistic-growth]" ) {
  
  //logisticGrowth<Type>* inst = logisticGrowth<Type>::getinstance();
  logisticGrowth<double>* inst = new logisticGrowth<double>(); //non-singleton usage
  
  model_traits<double>::parameter_vector u={5000, 6000, 7000, 5500};
  model_traits<double>::parameter_vector eta={0, 6250, 7200, 8050};
  double sigma = 0.1;
  double tau = 0.03;
  model_traits<double>::data_vector y={5050, 5950, 6950, 6000};
  
  auto nll_val = inst->calculateNll(u, eta, sigma, tau, y);
  double true_nll = 58.13528;
  
  REQUIRE( (nll_val - true_nll) <= 0.0001 );
  
  SECTION( "Different sigma test" ) {
    
    sigma = 0.3;
    nll_val = inst->calculateNll(u, eta, sigma, tau, y);
    REQUIRE( nll_val != true_nll );
  }
  
}

TEST_CASE( "Logistic growth model evaluate test", "[logistic-growth]" ) {
  
  //logisticGrowth<Type>* inst = logisticGrowth<Type>::getinstance();
  logisticGrowth<double>* inst = new logisticGrowth<double>(); //non-singleton usage
  
  model_traits<double>::parameter_vector theta = {-0.6931472, 9.2103404};
  double ln_sig = -2.302585;
  double ln_tau = -3.506558;
  
  model_traits<double>::parameter_vector u={5000, 6000, 7000, 5500};
  model_traits<double>::data_vector y={5050, 5950, 6950, 6000};
  
  double true_nll = 58.13528;
  
  inst->y = y;
  inst->theta = theta;
  inst->ln_sig = ln_sig;
  inst->ln_tau = ln_tau;
  inst->u = u;

  double nll_val = inst->evaluate();

  REQUIRE( (nll_val - true_nll) <= 0.0001 );
  

  SECTION( "Different ln_sig test" ) {
    
    inst->ln_sig = -5;
    nll_val = inst->evaluate();
    REQUIRE( nll_val != true_nll );
    
  }
  
}


