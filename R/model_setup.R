#runs TMB and tmbstan
runTMB <- function(dat,mod,seed=123,prType=0){
  results <- list()
  if(mod != 'spatial'){
    setupTMB(dll.name = 'stateSpace')
    Dat <- mkTMBdat(dat, mod, prType)
    Par <- do.call(paste0(mod,'.init'), args = list())
    Random <- 'u'
    if(mod == 'gompertz'){
      upr.limit <- rep(Inf, length(unlist(Par)))
      lwr.limit <- rep(-Inf, length(unlist(Par)))
    } else {
      lwr.limit <- c(log(0.01), log(0.01), -Inf, -Inf, rep(0.001, length(Par$u)))
      upr.limit <- c(log(1000), log(1000), rep(Inf, (length(Par$u)+2))) #+2 for ln_sig and ln_tau
    }
  }
  if(mod == 'spatial'){
    setupTMB(dll.name = 'spatial_poisson')
    inits <- mkSpatialInits(dat)
    Dat <- inits$Dat
    Par <- inits$Par()
    Random = 'omega'
    upr.limit <- rep(Inf, length(unlist(Par)))
    lwr.limit <- rep(-Inf, length(unlist(Par)))
  }
  set.seed(seed)
  obj <- MakeADFun(Dat, Par, random = Random)
  a<- Sys.time()
  tmb.mod <- nlminb( obj$par, obj$fn, obj$gr, 
                     lower = lwr.limit, upper = upr.limit )
  sdr <- sdreport(obj)
  b <- Sys.time()
  opt.par <- obj$env$last.par.best
  if(mod == 'gompertz'){
    results$tmb <- list(par.est = c(summary(sdr)[1:2,1],summary(sdr,'report')[,1]),
                        se.est = c(summary(sdr)[1:2,2],summary(sdr,'report')[,2]),
                        time = as.numeric(difftime(b,a, units = 'mins')), 
                        stan.time = NA, meanESS = NA, minESS = NA)
  }
  if(mod == 'spatial'){
    results$tmb <- list(par.est = c(summary(sdr,'fixed')[1,1],summary(sdr,'report')[,1]),
                        se.est = c(summary(sdr,'fixed')[1,2],summary(sdr,'report')[,2]),
                        time = as.numeric(difftime(b,a, units = 'mins')), 
                        stan.time = NA, meanESS = NA, minESS = NA)
  }
  if(mod == 'logistic'){
    results$tmb <- list(par.est = summary(sdr,'report')[,1],
                        se.est = summary(sdr,'report')[,2],
                        time = as.numeric(difftime(b,a, units = 'mins')), 
                        stan.time = NA, meanESS = NA, minESS = NA)
  }
  print('TMB model complete')
  
  a<- Sys.time()
  obj$env$data$hyperpars <- mkTMBdat(dat, mod, prType=1)$hyperpars
  tmbstan.mod <- try(tmbstan(obj, seed = seed, init = 'last.par.best', iter = 4000))
  b<- Sys.time()
  mon <- monitor(tmbstan.mod)
  
  if(is.null(summary(tmbstan.mod))){
    results$tmbstan <- list(par.est = c(theta = NA, theta = NA, ln_sig = NA, ln_tau = NA),
                            se.est = c(theta = NA, theta = NA, ln_sig = NA, ln_tau = NA),
                            time = as.numeric(difftime(b,a, units = 'mins')),
                            stan.time = NA, meanESS = NA, minESS = NA)
  } else {
    if(mod == 'logistic'){
      results$tmbstan <- list(par.est = exp(summary(tmbstan.mod)[[1]][1:4,1]),
                              se.est = summary(tmbstan.mod)[[1]][1:4,2], 
                              time = difftime(b,a, units = 'mins'), 
                              stan.time = as.numeric(difftime(b,a, units = 'mins')), 
                              meanESS = mean(mon$Bulk_ESS), minESS = min(mon$Bulk_ESS))
      warning('se in log space, not comparable to other model runs')
    } else {
      results$tmbstan <- list(par.est = c(summary(tmbstan.mod)[[1]][1:2,1],exp(summary(tmbstan.mod)[[1]][3:4,1]) ),
                              se.est = summary(tmbstan.mod)[[1]][1:4,2],
                              time = as.numeric(difftime(b,a, units = 'mins')), 
                              meanESS = mean(mon$Bulk_ESS), minESS = min(mon$Bulk_ESS))
      warning('sigma and tau SE in log space, not comparable to other model runs')
    }
  }

  results$inits <- opt.par
  print('tmbstan model complete')
  return(results)
}

#runs stan 
runSTAN <- function(dat,mod,prType,seed=123){
  if(prType == 0){
    hyperParameters <- list(
      hyperSig = 0, #dexp
      hyperTau = 0, #dexp
      hyperTheta1 = 0, #normal mean
      hyperTheta2 = 0 #normal sd
    )
    if(mod == 'logistic'){
      hyperParameters$hyperTheta1 <- c(0,0)
      hyperParameters$hyperTheta2 <- c(0,0)
    }
  }
  if(prType == 1){
    hyperParameters <- list(
      hyperSig = 0.1, #vague exp prior
      hyperTau = 0.1, #vague exp prior
      hyperTheta1 = 0, #normal mean
      hyperTheta2 = 100 #normal sd
    )
    if(mod == 'logistic'){
      hyperParameters$hyperTheta1 = c(-1,4) #lognormal
      hyperParameters$hyperTheta2 = c(5,4) #lognormal
    }
  }
  Dat <- mkSTANdat(dat,hyperParameters)
  Dat$prior_type <- prType
  file <- file.path('src', 'stan', paste0(mod, '.stan'))
  model <- cmdstan_model(file) #compile stan model into C++
  Par <- get(paste0(mod,'.init'))
  a <- Sys.time()
  fit <- model$sample(
    data = Dat,
    seed = seed,
    init = Par,
    iter_warmup = 2000, iter_sampling = 2000
  )
  b <- Sys.time()
  if(mod == 'gompertz') par.names <- c('theta', 'sigma', 'tau')
  if(mod == 'logistic') par.names <- c('r', 'K', 'sigma', 'tau')
  results <- list(par.est = as.vector(as.matrix(fit$summary(par.names, 'mean')[,2])),
                  se.est = as.vector(as.matrix(fit$summary(par.names, 'sd')[,2])),
                  time = as.numeric(difftime(b,a, units = 'mins')),
                  stan.time = as.numeric(difftime(b,a, units = 'mins')),
                  meanESS =  mean(fit$summary()$ess_bulk, na.rm=TRUE),
                  minESS =  min(fit$summary()$ess_bulk, na.rm=TRUE))
  return(results)
}