## install CRAN packages
# install.packages("tmbstan")
library(tmbstan)

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
Dat <- mkTMBdat(simdata, Mod)
Random <- 'u'
set.seed(123)
obj <- MakeADFun(Dat, logistic.init(), random = Random)
tmb.mod <- tmbstan.mod <- try(tmbstan(obj, seed = 123, init = 'last.par.best', iter = 4000))
