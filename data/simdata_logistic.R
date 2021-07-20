## Simulate data for simple and more complex state-space models:
## 1. Gompertz
## 2. Logistic growth model

gendat <- function(seed = NULL, N=100,
                   theta = c(2,.8), #gompertz: c(2,.8); logistic:c(.2,100)
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
  return(y)
}
