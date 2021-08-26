## install CRAN packages
# install.packages("TMB")
library(TMB)

source('data/simdata.R')
source('R/utils.R')

Mod <- 'logistic'
n <- 100

logistic.init <- function(){
  list(theta = c(log(0.5), log(80)), ln_sig=-1,ln_tau=-1,
       u = rep(1,n))
}

simdata <- gendat(seed=123,
                  N=n,
                  theta = c(0.2,100),
                  u1 = 4,
                  var = list(proc=0.01,obs=0.001),
                  mod.name = Mod)

setupTMB(dll.name = 'stateSpace')
Dat <- mkTMBdat(simdata, Mod, prType=0)
Random <- 'u'
set.seed(123)
obj <- MakeADFun(Dat, logistic.init(), random = Random)
tmb.mod <- nlminb( obj$par, obj$fn, obj$gr )
sdr <- sdreport(obj)

