library(rstan)
library(TMB)
library(tmbstan)
library(cmdstanr)
library(profvis)
library(INLA)

source('data/simdata.R')
source('src/R/utils.R')
source('src/R/model_setup.R')

#Implement simulation for n = 2^seq(5,11,1)
n.seq <- seq(5,11,1)
for(i in 1:length(n.seq)){
  n <- 2^n.seq[i]
  gompertz.results <- list()
  logistic.results <- list()
  
  gompertz.init <- function(){
    list(theta = c(0,0), ln_sig = 0,
         ln_tau = 0,
         u = rep(0,n))
  }
  
  
  #run simulation models
  #gompertz
  Mod <- 'gompertz'
  simdata <- gendat(seed=123,
                    N=n,
                    theta = c(2,0.8),
                    u1 = 4,
                    var = list(proc=0.1,obs=0.5),
                    mod.name = Mod)
  write.csv(data.frame(y=simdata), file = paste0('data/gompertz/gompertz', '_n', n, '.csv'))
  results <- runTMB(simdata,Mod)
  
  #modify and save for Julia
  inits <- c(alpha = unname(results$inits[1]), beta = unname(results$inits[2]), sigma = unname(exp(results$init[3])), 
             tau = unname(exp(results$inits[4])), u_init = unname(results$inits[5]), 
             results$inits[6:length(results$inits)])
  save(inits, file = paste0('data/gompertz/gompertzInits', '_n', n, '.RData'))
  gompertz.results <- list(tmb = results$tmb, tmbstan = results$tmbstan)
  
  # Init functions
  gompertz.init <- function(){
    list(theta = results$inits[1:2],
         ln_sig = results$inits[3],
         ln_tau = results$inits[4],
         u = results$inits[5:length(results$inits)])
  }
  
  #stan
  #use improper priors to compare with tmbstan
  gompertz.results$stanP0 <- runSTAN(simdata, Mod,0)
  #use vague priors
  gompertz.results$stanP1 <- runSTAN(simdata, Mod,1)
  
  save(gompertz.results, file = paste0('results/gompertz/gompertz', '_n', n, '.RData'))
  
  #Compare stan, tmbstan, tmb
  cbind(true=c(2,0.8,0.1,0.5),sapply(gompertz.results, function(x) x$par.est))
  sapply(gompertz.results, function(x) x$se.est)
  sapply(gompertz.results, function(x) x$time)
  sapply(gompertz.results, function(x) x$meanESS)
  sapply(gompertz.results, function(x) x$minESS)
  
  
  #logistic model
  Mod <- 'logistic'
  
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
  write.csv(data.frame(y=simdata), file = paste0('data/logistic/logistic', '_n', n, '.csv'))
  results <- runTMB(simdata,Mod) 
  #modify and save for Julia
  inits <- c(r = unname(exp(results$inits[1])), K = unname(exp(results$inits[2])), sigma = unname(exp(results$init[3])), 
             tau = unname(exp(results$inits[4])), u_init = unname(results$inits[5]), 
             results$inits[6:length(results$inits)])
  save(inits, file = paste0('data/logistic/logisticInits', '_n', n, '.RData'))
  logistic.results <-  list(tmb = results$tmb, tmbstan = results$tmbstan)
  
  
  #stan
  #use improper priors to compare with tmbstan
  logistic.results$stanP0 <- runSTAN(simdata, Mod,0)
  #use vague priors
  logistic.results$stanP1 <- runSTAN(simdata, Mod,1)
  
  save(logistic.results, file = paste0('results/logistic/logistic', '_n', n, '.RData'))
  #Compare rstan, tmbstan, tmb
  cbind(true=c(0.2,100,0.01,0.001),round(sapply(logistic.results, function(x) x$par.est),3))
  sapply(logistic.results, function(x) x$se.est)
  sapply(logistic.results, function(x) x$time)
  sapply(logistic.results, function(x) x$meanESS)
  sapply(logistic.results, function(x) x$minESS)
  
  #spatial
  Mod <- 'spatial'
  simdata <- gendat(seed=123,
                    N=n,
                    theta = c(2,50,0.75), #c(b0,Range,sp.var)
                    u1 = NA,
                    var = NA,
                    mod.name = Mod)
  spatial.results <- runTMB(simdata, Mod)
}
