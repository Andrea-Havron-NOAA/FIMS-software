# Functional tests that are related to src/tmb/stateSpace.cpp
library(testthat)
library(TMB)
library(here)

dir_path <- here::here("src", "tmb")

compile(file.path(dir_path, "stateSpace.cpp"))
dyn.load(dynlib(file.path(dir_path, "stateSpace")))
# dyn.unload(dynlib(file.path(dir_path, "stateSpace")))

test_that("State-space Gompertz growth model works correctly", {

  # Prepare expected values
  seed <- 1234
  N <- 1000
  theta <- c(2, 0.8)
  u1 <- 4
  var <- list(proc = 0.1, obs = 0.5)
  mod <- 0 # 0: Gompertz; 1: Logistic

  expected_val <- list(
    theta1 = theta[1],
    theta2 = theta[2],
    sigma_proc = sqrt(var$proc),
    tau_obs = sqrt(var$obs)
  )

  # Prepare input data for state-space Gompertz TMB
  set.seed(seed)
  eta <- u <- rep(0, N)
  u[1] <- u1
  for (i in 2:N) {
    eta[i] <- theta[1] + theta[2] * u[i - 1]
    u[i] <- rnorm(1, eta[i], sqrt(var$proc))
  }
  obs_y <- rnorm(N, u, sqrt(var$obs))

  # Run TMB
  data <- list(
    y = obs_y,
    mod = mod
  )

  parameters <- list(
    theta = c(0, 0),
    ln_sig = 0,
    ln_tau = 0,
    u = rep(0, N)
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
    theta1 = summary(sdr)[1, 1],
    theta2 = summary(sdr)[2, 1],
    sigma_proc = summary(sdr, "report")[1, 1],
    tau_obs = summary(sdr, "report")[2, 1]
  )

  # Expect success: relative errors in theta1, theta2, sigma_proc, and tau_obs are less than 10% or 50% (when expected value is really small) when there is no model mismatch

  for (i in 1:length(expected_val)) {
    re <- abs((object_val[[i]] - expected_val[[i]]) / expected_val[[i]])
    tolerance_val <- ifelse(expected_val[[i]] < 0.1, 0.5, 0.1)
    expect_success(expect_lt(object = re, expected = tolerance_val))
  }

  # Expect failure: relative errors in theta1, theta2, sigma_proc, and tau_obs are greater than 10% or 50% after modifying input data y

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

  object_val <- list(
    theta1 = summary(sdr)[1, 1],
    theta2 = summary(sdr)[2, 1],
    sigma_proc = summary(sdr, "report")[1, 1],
    tau_obs = summary(sdr, "report")[2, 1]
  )

  for (i in 1:length(expected_val)) {
    re <- abs((object_val[[i]] - expected_val[[i]]) / expected_val[[i]])
    tolerance_val <- ifelse(expected_val[[i]] < 0.1, 0.5, 0.1)
    expect_failure(expect_lt(object = re, expected = tolerance_val))
  }
})

test_that("State-space Logistic growth model works correctly", {

  # Prepare expect values
  seed <- 1234
  N <- 1000
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
    r = summary(sdr, "report")[1, 1],
    K = summary(sdr, "report")[2, 1],
    sigma_proc = summary(sdr, "report")[3, 1],
    tau_obs = summary(sdr, "report")[4, 1]
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
    r = summary(sdr, "report")[1, 1],
    K = summary(sdr, "report")[2, 1],
    sigma_proc = summary(sdr, "report")[3, 1],
    tau_obs = summary(sdr, "report")[4, 1]
  )

  for (i in 1:length(expected_val)) {
    re <- abs((object_val[[i]] - expected_val[[i]]) / expected_val[[i]])
    tolerance_val <- ifelse(expected_val[[i]] < 0.1, 0.5, 0.1)
    expect_failure(expect_lt(object = re, expected = tolerance_val))
  }
})
