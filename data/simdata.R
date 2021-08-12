## Simulate data FIMS architechture comparison
## Currently three models:
## 1. State Space Gompertz
## 2. State Space Logistic growth model
## 3. Spatial Poisson

library(mvtnorm)

gendat <- function(seed = NULL, N=100,
                   theta = c(2,.8), #gompertz: c(alpha=2,beta=.8); logistic:c(r=.2,K=100); spatial: c(b0 = 2,Range=50,SpVar=.75)
                   u1 = 4,
                   var = list(proc=0.1, obs=0.5), #gompertz: proc=0.1,obs=0.5; logistic: proc=0.01, obs=0.001
                   mod.name,
                   prior.type = 0
                   ){
 
  if(mod.name == 'gompertz'){
    set.seed(seed)
    eta <- u <-  rep(0,N)
    u[1] <- u1
    for(i in 2:N){
      eta[i] <- theta[1] + theta[2]*u[i-1]
      u[i] <- rnorm(1, eta[i], sqrt(var$proc))
    }
    y <- rnorm(N,u,sqrt(var$obs))
  }
  
  if(mod.name == 'logistic'){

    r <- theta[1]
    K <- theta[2]
   
    set.seed(seed)
    eta <- u <-  rep(0,N)
    u[1] <- u1
    for(i in 2:N){
      eta[i] <- log(u[i-1] + r*u[i-1]*(1-u[i-1]/K))
      u[i] <- rlnorm(1, eta[i], sqrt(var$proc))
    }
    y <- rlnorm(N,log(u),sqrt(var$obs))
  }
  
  if(mod.name == 'spatial'){
    set.seed(seed)
    grid.xy <- expand.grid(x=seq(0,100,1), y=seq(0,100,1))
    d <- as.matrix(dist(grid.xy, upper = TRUE))
    Omega <- sim.omega(theta[2], theta[3], d)
    y.sim <- rpois(nrow(grid.xy), exp(theta[1]+as.vector(Omega)))
    
    samp.idx <- sample(1:nrow(grid.xy), N)
    y <- data.frame(grid.xy[samp.idx,], z=y.sim[samp.idx])
  }
  
  
  return(y)
}

## function to simulate spatial MVN using a Matern covariance function
sim.omega <- function(Range, sig2, Dmat, Nu = 1){
  Kappa <- sqrt(8)/Range
  N <- dim(Dmat)[1]
  Sigma <- matrix(NA,N,N)
  #matern covariance
  Sigma <- sig2 * 2^(1-Nu)/gamma(Nu) * (Kappa * Dmat) * besselK(Kappa * Dmat, Nu)
  diag(Sigma) <- sig2
  omega <- rmvnorm(1, rep(0,N), Sigma, method = "chol")
  
  return(omega)
}
