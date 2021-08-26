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
Dat <- mkTMBdat(simdata, Mod, prType = 0)

set.seed(123)
Map0 <- list(ln_sig = factor(NA), ln_tau = factor(NA))
obj0 <- MakeADFun(Dat, logistic.init(), map = Map0)
tmb.mod0 <- nlminb( obj0$par, obj0$fn, obj0$gr )
Map1 <- list(theta = rep(factor(NA),2), u = rep(factor(NA),n))
obj1 <- MakeADFun(Dat, obj0$env$parList(), map = Map1)
tmb.mod1 <- nlminb( obj1$par, obj1$fn, obj1$gr )
sdr <- sdreport(obj1)
