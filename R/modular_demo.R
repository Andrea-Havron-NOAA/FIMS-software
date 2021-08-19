library(TMB)
compile("src/Rcpp/logisticGrowth.cpp", flags= "-w")
dyn.load(dynlib("src/Rcpp/logisticGrowth"))
y <- read.csv('data/logistic/logistic_n128.csv')$y

dat <- list(y=y)
par <- list(theta = c(log(0.5), log(80)), ln_sig=-1,ln_tau=-1,
            u = rep(1,length(y)))
obj <- MakeADFun(dat, par, random="u", DLL="logisticGrowth")
newtonOption(obj, smartsearch=FALSE)
obj$fn()
obj$gr()
system.time(opt <- nlminb(obj$par, obj$fn, obj$gr))
rep <- sdreport(obj)
exp(summary(rep, 'fixed')[,1])

