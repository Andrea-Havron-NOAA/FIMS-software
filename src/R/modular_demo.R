library(TMB)
compile("src/Rcpp/logisticGrowth.cpp", flags= "-w")
dyn.load(dynlib("src/Rcpp/logisticGrowth"))
y <- read.csv('data/logistic/logistic_n128.csv')$y

dat <- list(y=y)
par <- list(theta = 0,
            ln_sig = 0,
            ln_tau = 0,
            u = dat$y*0)
obj <- MakeADFun(dat, par, random="u", DLL="logisticGrowth")
newtonOption(obj, smartsearch=FALSE)
obj$fn()
obj$gr()
system.time(opt <- nlminb(obj$par, obj$fn, obj$gr))
rep <- sdreport(obj)
rep
