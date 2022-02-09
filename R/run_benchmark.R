## install CRAN packages
# install.packages("TMB", "tmbstan")

## install packages not on CRAN
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

library(TMB)
library(tmbstan)
library(cmdstanr)
library(R2admb)

source('data/simdata.R')
source('R/utils.R')
source('R/model_setup.R')

setupADMB('logisticGrowth', doRE=TRUE)
setupADMB('logisticGrowth_nonsep', doRE=TRUE)

#Implement simulation for n = 2^seq(5,11,1)
n.seq <- seq(5,12,1)
for(i in 1:length(n.seq)){
  n <- 2^n.seq[i]
  
  #logistic model
  Mod <- 'logistic'
  
  logistic.init <- function(){
    list(theta = c(log(0.5), log(80)), ln_sig=-1,ln_tau=-1,
         u = rep(1,n))
  }
  
  simdata <- gendat(seed=i,
                    N=n,
                    theta = c(0.2,100),
                    u1 = 4,
                    var = list(proc=0.01,obs=0.001),
                    mod.name = Mod)
  write.csv(data.frame(y=simdata), file = paste0('data/logistic/logistic', '_n', n, '.csv'))
  results <- runTMB(simdata,Mod) 
  #modify and save for Julia
  inits <- c(r = unname(exp(results$inits[1])), K = unname(exp(results$inits[2])), sigma = unname(exp(results$init[3])), 
             tau = unname(exp(results$inits[4])), u_init = unname(results$inits[5]), 
             results$inits[6:length(results$inits)])
  save(inits, file = paste0('data/logistic/logisticInits', '_n', n, '.RData'))
  logistic.results <-  list(tmb = results$tmb, tmbstan = results$tmbstan)
  
  # Init functions
  logistic.init <- function(){
    list(theta = results$inits[1:2],
         ln_sig = results$inits[3],
         ln_tau = results$inits[4],
         u = results$inits[5:length(results$inits)])
  }
  #stan
  #use vague priors
  logistic.results$stan <- runSTAN(simdata, Mod,1)
  
  save(logistic.results, file = paste0('results/logistic/logistic', '_n', n, '.RData'))
  
  
  
}


## ADMB
## Issues with admb so moving runs into its own loop

for(i in 1:length(n.seq)){
  n <- 2^n.seq[i]
  
  #logistic model
  Mod <- 'logistic'
  
  simdata <- gendat(seed=i,
                    N=n,
                    theta = c(0.2,100),
                    u1 = 4,
                    var = list(proc=0.01,obs=0.001),
                    mod.name = Mod)
  load(paste0('results/logistic/logistic', '_n', n, '.RData'))
  logistic.results$admb <- runADMB(simdata, Mod)
  save(logistic.results, file = paste0('results/logistic/logistic', '_n', n, '.RData'))
}

