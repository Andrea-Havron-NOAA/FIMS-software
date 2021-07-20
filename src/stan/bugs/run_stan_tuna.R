library(rstan)
source('data/tuna.R')

#run surplus production model - tuna example
dat <- list(N=23,
            C = tuna.dat$C,
            I = tuna.dat$I)

s <- stan(file = 'src/stan/bugs/SP.stan',
          data = dat, iter = 12000,
          init = surplus.init)

