n.seq <- seq(5,11,1)
library(tidyr)
library(magrittr)
library(ggplot2)
gompertz.theta1 <- gompertz.theta2 <- gompertz.sigma <- 
  gompertz.tau <- gompertz.meanESS <- gompertz.minESS <- gompertz.time <- 
  logistic.r <- logistic.K <- logistic.sigma <- logistic.tau <- 
  logistic.meanESS <- logistic.minESS <- logistic.time <-
    as.data.frame(matrix(0,4,length(n.seq), 
                       dimnames = list(c('tmb', 'tmbstan', 'stanP0', 'stanP1'),
                                       paste0('n_', 2^n.seq))))
for(i in 1:length(n.seq)){
  n <- 2^n.seq[i]
  #gompertz
  load(paste0('results/gompertz/gompertz', '_n', n, '.RData'))
  parEst <- sapply(gompertz.results, function(x) x$par.est)
  gompertz.theta1[,i] <- parEst[1,]
  gompertz.theta2[,i] <- parEst[2,]
  gompertz.sigma[,i] <- parEst[3,]
  gompertz.tau[,i] <- parEst[4,]
  gompertz.meanESS[,i] <- sapply(gompertz.results, function(x) x$meanESS)
  gompertz.minESS[,i] <- sapply(gompertz.results, function(x) x$minESS)
  gompertz.time[,i] <- sapply(gompertz.results, function(x) x$time)
  #logistic
  load(paste0('results/logistic/logistic', '_n', n, '.RData'))
  parEst <- sapply(logistic.results, function(x) x$par.est)
  logistic.r[,i] <- parEst[1,]
  logistic.K[,i] <- parEst[2,]
  logistic.sigma[,i] <- parEst[3,]
  logistic.tau[,i] <- parEst[4,]
  logistic.meanESS[,i] <- sapply(logistic.results, function(x) x$meanESS)
  logistic.meanESS$n <-
  logistic.minESS[,i] <- sapply(logistic.results, function(x) x$minESS)
  logistic.time[,i] <- sapply(logistic.results, function(x) x$time)
}

g.minESS <- gompertz.minESS %>% t() %>% as.data.frame() %>% 
            pivot_longer(., 1:4, names_to = 'model', values_to = 'minESS')
g.minESS$n <- rep(2^n.seq,each = 4)
ggplot(g.minESS, aes(x=n,y=minESS, col = model)) + geom_line() 
g.meanESS <- gompertz.meanESS %>% t() %>% as.data.frame() %>% 
  pivot_longer(., 1:4, names_to = 'model', values_to = 'meanESS')
g.meanESS$n <- rep(2^n.seq,each = 4)
g.time <- gompertz.time %>% t() %>% as.data.frame() %>% 
  pivot_longer(., 1:4, names_to = 'model', values_to = 'time')

ggplot(mapping = aes(x=g.meanESS$n,y=g.meanESS$meanESS/g.time$time, col = g.meanESS$model)) + geom_line() 

l.minESS <- logistic.minESS %>% t() %>% as.data.frame() %>% 
  pivot_longer(., 1:4, names_to = 'model', values_to = 'minESS')
l.minESS$n <- rep(2^n.seq,each = 4)
ggplot(l.minESS, aes(x=n,y=minESS, col = model)) + geom_line() 
l.meanESS <- logistic.meanESS %>% t() %>% as.data.frame() %>% 
  pivot_longer(., 1:4, names_to = 'model', values_to = 'meanESS')
l.meanESS$n <- rep(2^n.seq,each = 4)
l.time <- logistic.time %>% t() %>% as.data.frame() %>% 
  pivot_longer(., 1:4, names_to = 'model', values_to = 'time')
ggplot(mapping = aes(x=l.meanESS$n,y=l.meanESS$meanESS/l.time$time, col = l.meanESS$model)) + geom_line() 
