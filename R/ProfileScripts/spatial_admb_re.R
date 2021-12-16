## install CRAN packages
# install.packages("R2admb")
library(R2admb)
library(INLA)


source('data/simdata.R')
source('R/utils.R')

wd <- getwd()
setwd('src/admb')
compile_admb("spatial_gmrf", re = TRUE, verbose = TRUE)

Mod <- 'spatial'

simdata <- gendat(seed=5,
                  N=NA,
                  theta = c(2,50,0.75), #c(b0,Range,sp.var)
                  u1 = NA,
                  var = NA,
                  mod.name = 'spatial')
set.seed(1)

n <- 100
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
difftime(b,a)
admb.rep <-  readLines('spatial_gmrf.rep')
admb.rep
parm.rep <- function(f){
  rep.out <- strsplit(f,"=")
  parm.nm <- rep.out[[1]][1] 
  parm.val <- as.numeric(strsplit(rep.out[[1]][2], " ")[[1]]) 
  if(is.na(parm.val[1])) parm.val <- parm.val[-1]
  out <- list()
  out[[parm.nm]] = parm.val
  return(out)
}
#Kappa = 0.05656854
#log(Kappa) = -2.8723
#Tau = 8.143375
#log(Tau) = 2.097205
parm.est <- sapply(1:length(admb.rep), function(x) parm.rep(admb.rep[x]))
parm.est$`b `
u.scale <- parm.est$`u `[mesh$idx$loc]*sqrt(4*pi)*sqrt(8)/50*sqrt(.75) #u/tau
plot(u.scale, sampdata$u);abline(0,1)


setwd(wd)
