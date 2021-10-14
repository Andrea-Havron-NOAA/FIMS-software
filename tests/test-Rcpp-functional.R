# Functional tests that are related to src/Rcpp/rcppFuns.cpp
library(testthat)
library(Rcpp)
library(here)

dir_path <- here::here("src", "Rcpp")

Rcpp::sourceCpp(file.path(dir_path, "rcppFuns.cpp"))

test_that("rcppFuns dlognorm works correctly", {
  set.seed(1)
  expected_x <- rlnorm(100)
  expected_res <- dlnorm(expected_x, log = TRUE)
  expected_res_sum <- sum(expected_res)

  object_res <- vector(mode = "numeric", length = length(expected_x))
  for (i in seq_along(expected_x)) {
    object_res[i] <- dlognorm(expected_x[i], 0, 1, 1)
  }
  object_res_sum <- sum(object_res)

  expect_equal(object_res_sum, expected_res_sum)
})
