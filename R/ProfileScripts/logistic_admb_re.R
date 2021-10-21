## install CRAN packages
# install.packages("R2admb")
library(R2admb)


source('data/simdata.R')
source('R/utils.R')

wd <- getwd()
setwd('src/admb')

Mod <- 'logistic'
n <- 100

simdata <- gendat(seed=123,
                  N=n,
                  theta = c(0.2,100),
                  u1 = 4,
                  var = list(proc=0.01,obs=0.001),
                  mod.name = Mod)

 
Dat <- list(n=n, y = simdata)
Par <- list(ln_sig=-1,ln_tau=-1, ln_r = log(0.5), ln_K = log(80), u = rep(1,n) )
write_dat('logisticGrowth', Dat)
write_pin('logisticGrowth', Par)
compile_admb("logisticGrowth", re = TRUE, verbose = TRUE)
admb.mod <- run_admb("logisticGrowth", verbose = TRUE)
admb.rep <-  readLines('logisticGrowth.rep')
parm.rep <- function(f){
  rep.out <- strsplit(f,"=")
  parm.nm <- rep.out[[1]][1] 
  parm.val <- rep.out[[1]][2] 
  out <- list()
  out[[parm.nm]] = parm.val
  return(out)
}
sapply(1:length(admb.rep), function(x) parm.rep(admb.rep[x]))



setwd(wd)
