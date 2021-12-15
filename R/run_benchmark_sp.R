## install CRAN packages
# install.packages("TMB", "tmbstan")

## install packages not on CRAN
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

library(TMB)
library(tmbstan)
library(cmdstanr)
library(INLA)
library(R2admb)

source('data/simdata.R')
source('R/utils.R')
source('R/model_setup.R')



simdata <- gendat(seed=5,
                  N=NA,
                  theta = c(2,50,0.75), #c(b0,Range,sp.var)
                  u1 = NA,
                  var = NA,
                  mod.name = 'spatial')
png(filename = 'results/plots/spatialsim.png', height = 1200, width = 1500, res = 200)
ggplot(mapping = aes(x=simdata[,1], y=simdata[,2],fill=simdata[,3])) + 
  geom_raster() + xlab('') + ylab('') + theme_classic()
dev.off()
#Implement sampling for n = 2^seq(5,11,1)
n.seq <- seq(5,11,1)

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
  inits <- spatial.inits()
  save(inits,  file = paste0('data/spatial/spatialInits', '_n', n, '.RData'))
  spatial.results$stan <- runSTAN(simdata, Mod,1)
  
  save(spatial.results, file = paste0('results/spatial/spatial', '_n', n, '.RData'))

}
#run admb spatial
wd <- getwd()
setwd('src/admb')
results <- list()
parm.rep <- function(f){
  rep.out <- strsplit(f,"=")
  parm.nm <- rep.out[[1]][1] 
  parm.val <- as.numeric(strsplit(rep.out[[1]][2], " ")[[1]]) 
  if(is.na(parm.val[1])) parm.val <- parm.val[-1]
  out <- list()
  out[[parm.nm]] = parm.val
  return(out)
}
for(i in 1:length(n.seq)){
  n <- 2^n.seq[i]
  Mod <- 'spatial'
  set.seed(i)
  samp.idx <- sample(1:nrow(simdata), n)
  sampdata <- simdata[samp.idx,]
  mesh <- INLA::inla.mesh.create(sampdata[,1:2])
  spde <- INLA::inla.spde2.matern(mesh)
  m <- mesh$n
  
  Dat <- list(n=n, m = m,  y = sampdata$z,
              v_i = mesh$idx$loc, M0 = as.matrix(spde$param.inla$M0),
              M1 = as.matrix(spde$param.inla$M1), 
              M2 = as.matrix(spde$param.inla$M2),
              log_kappa = log(sqrt(8)/50), 
              log_tau = log( 1/( sqrt(4*pi)*sqrt(8)/50*sqrt(.75) )))
  Par <- list(b = 1,  u = rep(0,m) )
  write_dat('spatial_gmrf', Dat)
  write_pin('spatial_gmrf', Par)
  file.remove('spatial_gmrf.rep')
  a <- Sys.time()
  admb.mod <- run_admb("spatial_gmrf", verbose = TRUE)
  b <- Sys.time()
  admb.rep <-  readLines('spatial_gmrf.rep')
  parm.est <- sapply(1:length(admb.rep), function(x) parm.rep(admb.rep[x]))
  results <- list(parms = parm.est, time = difftime(b,a, units = 'min'))
  save(results, file = paste0('../../results/spatial/admb_re_n', n, '.RData'))
}
 setwd(wd)
  


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

