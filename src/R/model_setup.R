#runs TMB and tmbstan
runTMB <- function(dat,mod){
  results <- list()
  setupTMB()
  Dat <- mkTMBdat(dat, mod)
  Par <- do.call(paste0(mod,'.init'), args = list())
  obj <- MakeADFun(Dat, Par, random = 'u')
  a<- Sys.time()
  tmb.mod <- nlminb( obj$par, obj$fn, obj$gr )
  sdr <- sdreport(obj)
  b <- Sys.time()
  if(mod == 'gompertz'){
    results$tmb <- list(par.est = c(summary(sdr)[1:2,1],summary(sdr,'report')[,1]),
                        se.est = c(summary(sdr)[1:2,1],summary(sdr,'report')[,2]),
                        time = difftime(b,a, units = 'mins'))
  } else {
    results$tmb <- list(par.est = summary(sdr,'report')[,1],
                        se.est = summary(sdr,'report')[,2],
                        time = difftime(b,a, units = 'mins'))
  }
  print('TMB model complete')
  
  a<- Sys.time()
  tmbstan.mod <- try(tmbstan(obj, init = 'last.par.best', iter = 4000))
  b<- Sys.time()
  if(is.null(summary(tmbstan.mod))){
    results$tmbstan <- list(par.est = c(theta = NA, theta = NA, ln_sig = NA, ln_tau = NA),
                            se.est = c(theta = NA, theta = NA, ln_sig = NA, ln_tau = NA),
                            time = difftime(b,a, units = 'mins'))
  } else {
    if(mod == 'logistic'){
      results$tmbstan <- list(par.est = exp(summary(tmbstan.mod)[[1]][1:4,1]),
                              se.est = summary(tmbstan.mod)[[1]][1:4,2], 
                              time = difftime(b,a, units = 'mins'))
      warning('se in log space, not comparable to other model runs')
    } else {
      results$tmbstan <- list(par.est = c(summary(tmbstan.mod)[[1]][1:2,1],exp(summary(tmbstan.mod)[[1]][3:4,1]) ),
                              se.est = summary(tmbstan.mod)[[1]][1:4,2],
                              time = difftime(b,a, units = 'mins'))
      warning('sigma and tau SE in log space, not comparable to other model runs')
    }
  }
  opt.par <- obj$env$last.par.best
  results$inits <- opt.par
  print('tmbstan model complete')
  return(results)
}

#runs stan - not working yet
runSTAN <- function(dat,mod,hyp,
                    rstanArgs = list(
                      chains = 4, iter = 2000, warmup = floor(iter/2), thin = 1
                    )){
  results <- list()
  Dat <- mkSTANdat(dat,ptype=0,hyp)
  mod0 <- stan(file = paste0('src/stan/',mod,'.stan'),
               data = Dat, 
               chains = rstanArgs$chains,
               iter = rstanArgs$iter,
               warmup)
}