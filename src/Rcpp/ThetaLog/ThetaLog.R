setwd('src/Rcpp/ThetaLog')
library(TMB)
compile("ThetaLog.cpp", flags= "-w")
dyn.load(dynlib("ThetaLog"))

## Read data
Y <- scan("ThetaLog.dat", skip=3, quiet=TRUE)
data <- list(Y=Y)

## Parameter initial guess
parameters <- list(
  X = data$Y*0,
  logr0 = 0,
  logtheta = 0,
  logK = 6,
  logQ = 0,
  logR = 0
)

## Fit model
obj <- MakeADFun(data, parameters, random="X", DLL="ThetaLog")
newtonOption(obj, smartsearch=FALSE)
obj$fn()
obj$gr()
system.time(opt <- nlminb(obj$par, obj$fn, obj$gr))
rep <- sdreport(obj)
rep
