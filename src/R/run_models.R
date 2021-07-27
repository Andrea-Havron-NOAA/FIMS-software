library(rstan)
library(TMB)
library(tmbstan)
library(cmdstanr)

source('data/simdata_logistic.R')
source('src/R/utils.R')
source('src/R/model_setup.R')

gompertz.results <- list()
logistic.results <- list()

# Init functions
gompertz.init <- function(){
  list(theta = c(0,0), ln_sig = 0, ln_tau = 0,
       u = rep(0,n))
}
logistic.init <- function(){
  list(theta = c(log(0.5), log(80)), ln_sig=-1,ln_tau=-1,
       u = rep(1,n))
}

n <- 100
#run simulation models
#gompertz
Mod <- 'gompertz'
simdata <- gendat(seed=123,
                  N=n,
                  theta = c(2,0.8),
                  u1 = 4,
                  var = list(proc=0.1,obs=0.5),
                  mod.name = Mod)
gompertz.results <- runTMB(simdata,Mod)

#stan
#use improper priors to compare with tmbstan
a <- Sys.time()
hyperParameters <- list(
  hyperSig = c(0,0), #dgamma
  hyperTau = c(0,0), #dgamma
  hyperTheta1 = 0, #normal mean
  hyperTheta2 = 0 #normal sd
)
Dat <- mkSTANdat(simdata, hyperParameters) 
file <- file.path('src', 'stan', 'gompertz.stan')
mod <- cmdstan_model(file) #compile stan model into C++
fit.p0 <- mod$sample(
  data = Dat,
  init = gompertz.init,
  iter_warmup = 4000, iter_sampling = 5000
)

b <- Sys.time()
gompertz.results$stanP0 <- list(par.est = as.vector(as.matrix(fit.p0$summary(c('theta', 'ln_sig', 'ln_tau'), 'mean')[,2])),
                                se.est = as.vector(as.matrix(fit.p0$summary(c('theta', 'ln_sig', 'ln_tau'), 'sd')[,2])),
                                time = difftime(b,a, units = 'mins'))
gompertz.results$stanP0$par.est[3:4] <- exp(gompertz.results$stanP0$par.est[3:4])
#use proper vague priors
a <- Sys.time()
hyperParameters <- list(
  hyperSig = c(0.001,0.001), #dgamma
  hyperTau = c(0.001,0.001), #dgamma
  hyperTheta1 = 0, #normal mean
  hyperTheta2 = 100 #normal sd
)
Dat <- mkSTANdat(simdata, hyperParameters)
Dat$prior_type <- 1
fit.p1 <- mod$sample(
  data = Dat,
  init = gompertz.init,
  iter_warmup = 4000, iter_sampling = 5000
)
b <- Sys.time()
gompertz.results$stanP1 <-  list(par.est = as.vector(as.matrix(fit.p1$summary(c('theta', 'ln_sig', 'ln_tau'), 'mean')[,2])),
                                 se.est = as.vector(as.matrix(fit.p1$summary(c('theta', 'ln_sig', 'ln_tau'), 'sd')[,2])),
                                 time = difftime(b,a, units = 'mins'))
gompertz.results$stanP1$par.est[3:4] <- exp(gompertz.results$stanP1$par.est[3:4])
save(gompertz.results, file = 'results/gompertz.RData')
#Compare stan, tmbstan, tmb
cbind(true=c(2,0.8,0.1,0.5),sapply(gompertz.results, function(x) x$par.est))
sapply(gompertz.results, function(x) x$se.est)
sapply(gompertz.results, function(x) x$time)


#logistic model
Mod <- 'logistic'

simdata <- gendat(seed=123,
                  N=100,
                  theta = c(0.2,100),
                  u1 = 4,
                  var = list(proc=0.01,obs=0.001),
                  mod.name = Mod)

logistic.results <- runTMB(simdata,Mod) 


#stan
#use improper priors to compare with tmbstan
a <- Sys.time()
hyperParameters <- list(
  hyperSig = c(0,0), #dgamma
  hyperTau = c(0,0), #dgamma
  hyperTheta1 = c(0,0), #normal
  hyperTheta2 = c(0,0) #normal 
)
Dat <- mkSTANdat(simdata, hyperParameters) 
file <- file.path('src', 'stan', 'logistic.stan')
mod <- cmdstan_model(file) #compile stan model into C++
fit.p0 <- mod$sample(
  data = Dat,
  init = logistic.init,
  iter_warmup = 4000, iter_sampling = 5000
)
b <- Sys.time()
logistic.results$stanP0 <-  list(par.est = as.vector(as.matrix(fit.p0$summary(c('r', 'K', 'sigma', 'tau'), 'mean')[,2])),
                                 se.est = as.vector(as.matrix(fit.p0$summary(c('r', 'K', 'sigma', 'tau'), 'sd')[,2])),
                                 time = difftime(b,a, units = 'mins'))

#use proper vague priors
a <- Sys.time()
hyperParameters <- list(
  hyperSig = c(0.001,0.001), #dgamma
  hyperTau = c(0.001,0.001), #dgamma
  hyperTheta1 = c(-1,4), #lognormal
  hyperTheta2 = c(5,4) #lognormal
)
Dat <- mkSTANdat(simdata, hyperParameters)
Dat$prior_type <- 1
fit.p1 <- mod$sample(
  data = Dat,
  init = logistic.init,
  iter_warmup = 4000, iter_sampling = 5000
)
b <- Sys.time()
logistic.results$stanP1 <- list(par.est = as.vector(as.matrix(fit.p1$summary(c('r','K', 'sigma', 'tau'), 'mean')[,2])),
                                se.est = as.vector(as.matrix(fit.p1$summary(c('r', 'K', 'sigma', 'tau'), 'sd')[,2])),
                                time = difftime(b,a, units = 'mins'))

save(logistic.results, file = 'results/logistic.RData')
#Compare rstan, tmbstan, tmb
cbind(true=c(0.2,100,0.01,0.001),round(sapply(logistic.results, function(x) x$par.est),3))
sapply(logistic.results, function(x) x$se.est)
sapply(logistic.results, function(x) x$time)
