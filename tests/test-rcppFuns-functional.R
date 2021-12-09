# Functional tests that are related to src/Rcpp/rcppFuns.cpp
library(testthat)
library(Rcpp)
library(here)

dir_path <- here::here("src", "Rcpp")

Rcpp::sourceCpp(file.path(dir_path, "rcppFuns.cpp"))

test_that(
  "rcppFuns dlognorm works with different input values", {
  set.seed(1)
  x_vec <- rlnorm(100)
  expected_val <- dlnorm(x_vec, 0.0, 1.0, log = TRUE)

  # Expect object values match expected values
  object_val <- vector(mode = "numeric", length = length(expected_val))
  for (i in seq_along(expected_val)) {
    object_val[i] <- dlognorm(x_vec[i], 0.0, 1.0, TRUE)
  }

  expect_equal(object_val, expected_val)
  
  # Expect object value matches expected value with a different meanlog value
  expected_val <- dlnorm(x=5.0, meanlog = 10.0, sdlog = 1.0, log=TRUE)
  object_val <- dlognorm(5.0, 10.0, 1.0, TRUE)
  expect_equal(object_val, expected_val)
  
  # Expect object value matches expected value with a different sdlog value
  expected_val <- dlnorm(x=5.0, meanlog = 10.0, sdlog = 4.0, log=TRUE)
  object_val <- dlognorm(5.0, 10.0, 4.0, TRUE)
  expect_equal(object_val, expected_val)
  
  # Expect object value matches expected value with a different give_log value
  expected_val <- dlnorm(x=5.0, meanlog = 10.0, sdlog = 4.0, log=FALSE)
  object_val <- dlognorm(5.0, 10.0, 4.0, FALSE)
  expect_equal(object_val, expected_val)
  
  # Expect object value does not match expected value with a different x value 
  expected_val <- dlnorm(x=5.0, meanlog = 10.0, sdlog = 4.0, log=TRUE)
  object_val <- dlognorm(4.0, 10.0, 4.0, TRUE)
  expect_false(isTRUE(object_val==expected_val))
  
  # Expect object value matches expected value when using default give_log value (FALSE)
  expected_val <- dlnorm(x=5.0, meanlog = 10.0, sdlog = 4.0, log=FALSE)
  object_val <- dlognorm(5.0, 10.0, 4.0)
  expect_equal(object_val, expected_val)
  
  
  }
)

test_that(
  "rcppFuns dlognorm works with different types of inputs", {
    
    # Expect object value matches expected value when using numeric values with digit=0
    expected_val <- dlnorm(x=5, meanlog = 10, sdlog = 4, log=TRUE)
    object_val <- dlognorm(5, 10, 4, TRUE)
    expect_equal(object_val, expected_val)
    
    # Expect object value matches expected value when using numeric value for give_log parameter
    expected_val <- dlnorm(x=5, meanlog = 10, sdlog = 4, log=0)
    object_val <- dlognorm(5, 10, 4, 0)
    expect_equal(object_val, expected_val)
    
    # Expect dlognorm does not work with a vector of x values
    expect_error(dlognorm(rlnorm(100), 0.0, 1.0, 0))
    
    # Expect dlognorm does not work with a non 0/1 double value
    expect_error(dlognorm(5.0, 0.0, 1.0, 0.5))
    expect_error(dlognorm(5.0, 0.0, 1.0, 3.5))
    expect_error(dlognorm(5.0, 0.0, 1.0, -3.5))
    
  }
)
