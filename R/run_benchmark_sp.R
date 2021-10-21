## install CRAN packages
# install.packages("TMB", "tmbstan")

## install packages not on CRAN
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

library(TMB)
library(tmbstan)
library(cmdstanr)
library(INLA)

source('data/simdata.R')
source('R/utils.R')
source('R/model_setup.R')



simdata <- gendat(seed=5,
                  N=NA,
                  theta = c(2,50,0.75), #c(b0,Range,sp.var)
                  u1 = NA,
                  var = NA,
                  mod.name = 'spatial')

#Implement sampling for n = 2^seq(5,11,1)
#n.seq <- seq(5,11,1)
n.seq = 6
for(i in 1:length(n.seq)){
  n <- 2^n.seq[i]

  set.seed(i)
  samp.idx <- sample(1:nrow(simdata), n)
  sampdata <- simdata[samp.idx,]
 
  Mod <- 'spatial'

  spatial.results <- runTMB(sampdata, Mod, 1)
  save(spatial.results, file = paste0('results/spatial/spatial_n',n,'.RData'))
  
  
  #stan
  # Init functions
  spatial.inits <- function(){
    list(b0 = spatial.results$inits[1], 
         omega = spatial.results$inits[2:length(spatial.results$inits)])
  }
  spatial.results$stan <- runSTAN(simdata, Mod,1)
  
  save(spatial.results, file = paste0('results/spatial/spatial', '_n', n, '.RData'))
  #Compare rstan, tmbstan, tmb
  cbind(true=c(0.2,100,0.01,0.001),round(sapply(logistic.results, function(x) x$par.est),3))
  sapply(logistic.results, function(x) x$se.est)
  sapply(logistic.results, function(x) x$time)
  sapply(logistic.results, function(x) x$meanESS)
  sapply(logistic.results, function(x) x$minESS)
  
}
library(magrittr)
library(tidyr)
plot.res <- c()
for(i in 1:length(n.seq)){
  n <- 2^n.seq[i]
  load( paste0('results/logistic/logistic', '_n', n, '.RData'))
  plot.res <- rbind(plot.res, sapply(logistic.results, function(x) x$meanESS)[2:4]/
    sapply(logistic.results, function(x) x$time)[2:4])
  plot.res <- rbind(plot.res, sapply(logistic.results, function(x) x$minESS)[2:4])
  plot.res <- rbind(plot.res, sapply(logistic.results, function(x) x$time)[2:4])
}
colnames(plot.res) <- names(logistic.results)[2:4]
plot.res %<>% as.data.frame()
plot.res$metric <- rep(c('MCMC efficiency', 'min ESS', 'time'), length(n.seq))
plot.res$nsamp <- rep(2^n.seq, each = 3)
plot.res %>% 
  pivot_longer(., 1:3, names_to = 'model', values_to = 'value') %>%
  ggplot(., aes(x=nsamp, y=value,col=model)) + geom_line() + 
  theme_classic() + facet_wrap(~metric, scales = 'free', ncol=1)

