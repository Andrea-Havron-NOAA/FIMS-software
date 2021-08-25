library(TMB)
compile("src/Rcpp/logisticGrowth.cpp", flags= "-w")
dyn.load(dynlib("src/Rcpp/logisticGrowth"))
simdata <- gendat(seed=123,
                  N=100,
                  theta = c(0.2,100),
                  u1 = 4,
                  var = list(proc=0.01,obs=0.001),
                  mod.name = 'logistic')

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

