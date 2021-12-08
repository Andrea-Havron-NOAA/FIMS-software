# Functional tests that are related to src/Rcpp/logisticGrowth.cpp
library(testthat)
library(TMB)
library(here)

dir_path <- here::here("src", "Rcpp")

compile(file.path(dir_path, "logisticGrowth.cpp"), flags= "-DTMB_MODEL -w")

dyn.load(dynlib(file.path(dir_path, "logisticGrowth")))
# dyn.unload(dynlib(file.path(dir_path, "logisticGrowth")))

test_that("Logistic growth model returns correct outputs", {

  # Prepare expect values
  seed <- 1234
  N <- 2^12
  r <- 0.2
  K <- 100
  u1 <- 4
  var <- list(proc = 0.01, obs = 0.001)
  mod <- 1 # 0: Gompertz; 1: Logistic

  expected_val <- list(
    r = r,
    K = K,
    sigma_proc = sqrt(var$proc),
    tau_obs = sqrt(var$obs)
  )

  # Prepare input data for state-space Gompertz TMB
  set.seed(seed)
  eta <- u <- rep(0, N)
  u[1] <- u1
  for (i in 2:N) {
    eta[i] <- log(u[i - 1] + r * u[i - 1] * (1 - u[i - 1] / K))
    u[i] <- rlnorm(1, eta[i], sqrt(var$proc))
  }
  obs_y <- rlnorm(N, log(u), sqrt(var$obs))

  # Run TMB
  data <- list(
    y = obs_y,
    mod = mod
  )

  parameters <- list(
    theta = c(log(0.5), log(80)),
    ln_sig = -1,
    ln_tau = -1,
    u = rep(1, N)
  )

  random <- "u"

  obj <- MakeADFun(
    data = data,
    parameters = parameters,
    random = random
  )
  
  tmb.mod <- nlminb(obj$par, obj$fn, obj$gr)

  sdr <- sdreport(obj)
  
  # Get TMB results
  object_val <- list(
    r = exp(summary(sdr)[1,1]),
    K = exp(summary(sdr)[2,1]),
    sigma_proc = exp(summary(sdr)[3,1]),
    tau_obs = exp(summary(sdr)[4,1])
  )

  # Expect success: relative errors in r, K, sigma_proc, and tau_obs are less than 10% or 50% when there is no model mismatch

  for (i in 1:length(expected_val)) {
    re <- abs((object_val[[i]] - expected_val[[i]]) / expected_val[[i]])
    tolerance_val <- ifelse(expected_val[[i]] < 0.1, 0.5, 0.1)
    expect_success(expect_lt(object = re, expected = tolerance_val))
  }

  # Expect failure: relative errors in r, K, sigma_proc, and tau_obs are greater than 10% or 50% after modifying input data y

  data <- list(
    y = obs_y + runif(length(obs_y), 1, 100),
    mod = mod
  )

  obj <- MakeADFun(
    data = data,
    parameters = parameters,
    random = random
  )

  tmb.mod <- nlminb(obj$par, obj$fn, obj$gr)

  sdr <- sdreport(obj)

  # Get TMB results
  object_val <- list(
    r = exp(summary(sdr)[1,1]),
    K = exp(summary(sdr)[2,1]),
    sigma_proc = exp(summary(sdr)[3,1]),
    tau_obs = exp(summary(sdr)[4,1])
  )

  for (i in 1:length(expected_val)) {
    re <- abs((object_val[[i]] - expected_val[[i]]) / expected_val[[i]])
    tolerance_val <- ifelse(expected_val[[i]] < 0.1, 0.5, 0.1)
    expect_failure(expect_lt(object = re, expected = tolerance_val))
  }
})
