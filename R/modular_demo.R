library(TMB)
source('data/simdata.R')
#compile("src/Rcpp/logisticGrowth.cpp", flags= "-DTMB_MODEL -w")
TMB::compile("src/Rcpp/logisticGrowth.cpp")
dyn.load(dynlib("src/Rcpp/logisticGrowth"))
y <- gendat(seed=123,
                  N=100,
                  theta = c(0.2,100),
                  u1 = 4,
                  var = list(proc=0.01,obs=0.001),
                  mod.name = 'logistic')
dat <- list(y=y)
par <- list(theta = c(log(0.5), log(80)), ln_sig=-1,ln_tau=-1,
            u = rep(1,length(y)))
a <- Sys.time()
obj <- MakeADFun(dat, par, random="u", DLL="logisticGrowth")
newtonOption(obj, smartsearch=FALSE)
obj$fn()
obj$gr()
system.time(opt <- nlminb(obj$par, obj$fn, obj$gr))
sdr <- sdreport(obj)
b <- Sys.time()
difftime(b,a)
summary(sdr, 'fixed')
summary(sdr, 'random')[,1]
summary(sdr, 'report')

report <- obj$report()
report$u
report$K
report$r

dyn.unload(dynlib("src/Rcpp/logisticGrowth"))

