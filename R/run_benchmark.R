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

#Implement simulation for n = 2^seq(5,11,1)
n.seq <- seq(5,12,1)
for(i in 1:length(n.seq)){
  n <- 2^n.seq[i]
  # 
  # gompertz.init <- function(){
  #   list(theta = c(0,0), ln_sig = 0,
  #        ln_tau = 0,
  #        u = rep(0,n))
  # }
  # 
  # 
  # #run simulation models
  # #gompertz
  # Mod <- 'gompertz'
  # simdata <- gendat(seed=i,
  #                   N=n,
  #                   theta = c(2,0.8),
  #                   u1 = 4,
  #                   var = list(proc=0.1,obs=0.5),
  #                   mod.name = Mod)
  # write.csv(data.frame(y=simdata), file = paste0('data/gompertz/gompertz', '_n', n, '.csv'))
  # results <- runTMB(simdata,Mod)
  # 
  # #modify and save for Julia
  # inits <- c(alpha = unname(results$inits[1]), beta = unname(results$inits[2]), sigma = unname(exp(results$init[3])),
  #            tau = unname(exp(results$inits[4])), u_init = unname(results$inits[5]),
  #            results$inits[6:length(results$inits)])
  # save(inits, file = paste0('data/gompertz/gompertzInits', '_n', n, '.RData'))
  # gompertz.results <- list(tmb = results$tmb, tmbstan = results$tmbstan)
  # 
  # # Init functions
  # gompertz.init <- function(){
  #   list(theta = results$inits[1:2],
  #        ln_sig = results$inits[3],
  #        ln_tau = results$inits[4],
  #        u = results$inits[5:length(results$inits)])
  # }
  # 
  # #stan
  # # #use improper priors to compare with tmbstan
  # # gompertz.results$stanP0 <- runSTAN(simdata, Mod,0)
  # #use vague priors
  # gompertz.results$stan <- runSTAN(simdata, Mod, 1)
  # 
  # save(gompertz.results, file = paste0('results/gompertz/gompertz', '_n', n, '.RData'))
  # 
  # #Compare stan, tmbstan, tmb
  # cbind(true=c(2,0.8,0.1,0.5),sapply(gompertz.results, function(x) x$par.est))
  # sapply(gompertz.results, function(x) x$se.est)
  # sapply(gompertz.results, function(x) x$time)
  # sapply(gompertz.results, function(x) x$meanESS)
  # sapply(gompertz.results, function(x) x$minESS)


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
  
  
  #admb
  logistic.results$admb <- runADMB(simdata, Mod)
  save(logistic.results, file = paste0('results/logistic/logistic', '_n', n, '.RData'))
  
  #Compare rstan, tmbstan, tmb
  cbind(true=c(0.2,100,sqrt(0.01),sqrt(0.001)),round(sapply(logistic.results, function(x) x$par.est),3))
  sapply(logistic.results, function(x) x$se.est)
  sapply(logistic.results, function(x) x$time)
  sapply(logistic.results, function(x) x$meanESS)
  sapply(logistic.results, function(x) x$minESS)
  
  
  
}

library(ggplot2)
library(magrittr)
library(tidyr)
plot.res <- c()
tmb.est <- tmbstan.est <- stan.est <- c()
tmb.se <- tmbstan.se <- stan.se <- c()
for(i in 1:length(n.seq)){
  n <- 2^n.seq[i]
  load( paste0('results/logistic/logistic', '_n', n, '.RData'))
  plot.res <- rbind(plot.res, sapply(logistic.results, function(x) x$meanESS)[2:3]/
    sapply(logistic.results, function(x) x$time)[2:3])
  plot.res <- rbind(plot.res, sapply(logistic.results, function(x) x$minESS)[2:3])
  plot.res <- rbind(plot.res, sapply(logistic.results, function(x) x$time)[2:3])
  tmb.est <- rbind(tmb.est, logistic.results$tmb$par.est)
  tmbstan.est <- rbind(tmbstan.est, logistic.results$tmbstan$par.est)
  stan.est <- rbind(stan.est, logistic.results$stan$par.est)
  tmb.se <- rbind(tmb.se, logistic.results$tmb$se.est)
  tmbstan.se <- rbind(tmbstan.se, logistic.results$tmbstan$se.est)
  stan.se <- rbind(stan.se, logistic.results$stan$se.est)
}
colnames(plot.res) <- names(logistic.results)[2:3]
plot.res %<>% as.data.frame()
plot.res$metric <- rep(c('MCMC efficiency', 'min ESS', 'time'), length(n.seq))
plot.res$nsamp <- rep(2^n.seq, each = 3)
png(filename = 'results/plots/logistic.png', width = 1200, height = 500, res=200)
plot.res %>% 
  pivot_longer(., 1:2, names_to = 'model', values_to = 'value') %>%
  ggplot(., aes(x=nsamp, y=value,col=model)) + geom_line() + 
  theme_classic() + facet_wrap(~metric, scales = 'free', nrow=1)
dev.off()

tmb.est <- as.data.frame(tmb.est)
tmbstan.est <- as.data.frame(tmbstan.est)
stan.est <- as.data.frame(stan.est)
tmb.est$model <- rep('tmb-re',8)
tmb.est$nsamp <- 2^n.seq
tmbstan.est$model <- rep('tmbstan',8)
tmbstan.est$nsamp <- 2^n.seq
stan.est$model <- rep('stan',8)
stan.est$nsamp <- 2^n.seq
colnames(stan.est) <- colnames(tmb.est)
colnames(tmbstan.est) <- colnames(tmb.est)

est.res <- rbind(tmb.est, tmbstan.est, stan.est)
est.res %<>% pivot_longer(.,1:4,names_to = 'parm', values_to = 'est') %>% as.data.frame()
levels(est.res$model)
png(filename = 'results/plots/logistic_est.png', width = 1000, height = 1000, res = 200)
ggplot(est.res, aes(x=nsamp, y = est)) + 
  facet_grid(parm~model, scales = "free_y") + 
  geom_point() + theme_bw() +
  geom_hline(data=dplyr::filter(est.res, parm=="K"), aes(yintercept=100), color = 'blue')+
  geom_hline(data=dplyr::filter(est.res, parm=="r"), aes(yintercept=0.2), color = 'blue')+
  geom_hline(data=dplyr::filter(est.res, parm=="sigma"), aes(yintercept=sqrt(0.01)), color = 'blue')+
  geom_hline(data=dplyr::filter(est.res, parm=="tau"), aes(yintercept=sqrt(0.001)), color = 'blue')
dev.off()

