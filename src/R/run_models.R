library(rstan)
library(TMB)
library(tmbstan)
source('data/simdata_logistic.R')
source('src/R/utils.R')
source('src/R/model_setup.R')

gompertz.results <- list()
logistic.results <- list()

# Init functions
gompertz.init <- function(N=Dat$N){
  list(theta = c(0,0), ln_sig = 0, ln_tau = 0,
       u = rep(0,N))
}
logistic.init <- function(N=Dat$N){
  list(theta = c(log(0.5), log(80)), ln_sig=-1,ln_tau=-1,
       u = rep(1,N))
}

#run simulation models
#gompertz
Mod <- 'gompertz'
simdata <- gendat(seed=123,
                  N=100,
                  theta = c(2,0.8),
                  u1 = 4,
                  var = list(proc=0.1,obs=0.5),
                  mod.name = Mod)
gompertz.results <- runTMB(simdata,Mod)

#rstan
#use improper priors to compare with tmbstan
a <- Sys.time()
hyperParameters <- list(
  hyperSig = c(0,0), #dgamma
  hyperTau = c(0,0), #dgamma
  hyperTheta1 = 0, #normal mean
  hyperTheta2 = 0 #normal sd
)
Dat <- mkSTANdat(simdata, hyperParameters) 
stan.p0 <- stan(file = 'src/stan/gompertz.stan',
                 data = Dat, iter = 9000,
                 init = gompertz.init)
b <- Sys.time()
gompertz.results$stanP0 <- list(par.est = summary(stan.p0)[[1]][1:4,1],
                                se.est = summary(stan.p0)[[1]][1:4,2],
                                time = difftime(b,a, units = 'mins'))
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
stan.p1 <- stan(file = 'src/stan/gompertz.stan',
                    data = Dat, iter = 9000,
                    init = gompertz.init)
b <- Sys.time()
gompertz.results$stanP1 <- list(par.est = summary(stan.p1)[[1]][1:4,1],
                                se.est = summary(stan.p1)[[1]][1:4,2],
                                time = difftime(b,a, units = 'mins'))

save(gompertz.results, file = 'results/gompertz.RData')
#Compare rstan, tmbstan, tmb
cbind(true=c(2,0.8,log(0.1),log(0.5)),sapply(gompertz.results, function(x) x$par.est))
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

logistic.results <- runTMB(simdata,Mod) #not working


#rstan
#use improper priors to compare with tmbstan
a <- Sys.time()
hyperParameters <- list(
  hyperSig = c(0,0), #dgamma
  hyperTau = c(0,0), #dgamma
  hyperTheta1 = c(0,0), #normal
  hyperTheta2 = c(0,0) #normal 
)
Dat <- mkSTANdat(simdata, hyperParameters) 
stan.p0 <- stan(file = 'src/stan/logistic.stan',
                    data = Dat, iter = 9000,
                    init = logistic.init)
b <- Sys.time()
logistic.results$stanP0 <- list(par.est = summary(stan.p0)[[1]][1:4,1],
                                se.est = summary(stan.p0)[[1]][1:4,2],
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
stan.p1 <- stan(file = 'src/stan/logistic.stan',
                    data = Dat, iter = 9000,
                    init = logistic.init)
b <- Sys.time()
logistic.results$stanP1 <- list(par.est = summary(stan.p1)[[1]][1:4,1],
                                se.est = summary(stan.p1)[[1]][1:4,2],
                                time = difftime(b,a, units = 'mins'))

save(logistic.results, file = 'results/logistic.RData')
#Compare rstan, tmbstan, tmb
cbind(true=c(log(0.2),log(100),log(0.01),log(0.001)),sapply(logistic.results, function(x) x$par.est))
sapply(logistic.results, function(x) x$se.est)
sapply(logistic.results, function(x) x$time)

