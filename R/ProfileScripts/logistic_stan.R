## install packages not on CRAN
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

library(cmdstanr)

source('data/simdata.R')
source('R/model_setup.R')
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

runSTAN(simdata, Mod,0)