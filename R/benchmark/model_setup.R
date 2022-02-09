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
    inits <- mkSpatialInits(dat,0,'tmb')
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
    results$tmb <- list(par.est = c(summary(sdr,'fixed')[1,1]),#,summary(sdr,'report')[,1]),
                        se.est = c(summary(sdr,'fixed')[1,2]),#,summary(sdr,'report')[,2]),
                        omega.est = summary(sdr, 'random')[,1],
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
  
  if(mod == 'spatial'){
    if(prType==1){
      obj$env$data$prior_type <- 1
      new.dat <- mkSpatialInits(dat,1, 'tmb')
      obj$env$data$kap_tau_pr_mu <- new.dat$Dat$kap_tau_pr_mu
      obj$env$data$kap_tau_pr_var <- new.dat$Dat$kap_tau_pr_var
    }
  } else {
    obj$env$data$hyperpars <- mkTMBdat(dat, mod, prType=1)$hyperpars
  }
  a<- Sys.time()
  tmbstan.mod <- try(tmbstan(obj, seed = seed, init = 'last.par.best', iter = 4000))
  b<- Sys.time()
  mon <- monitor(tmbstan.mod)
  
  if(is.null(summary(tmbstan.mod))){
    if(mod != 'spatial'){
      results$tmbstan <- list(par.est = c(theta = NA, theta = NA, ln_sig = NA, ln_tau = NA),
                              se.est = c(theta = NA, theta = NA, ln_sig = NA, ln_tau = NA),
                              time = as.numeric(difftime(b,a, units = 'mins')),
                              stan.time = NA, meanESS = NA, minESS = NA)
    } else {
      results$tmbstan <- list(par.est = c(b0 = NA),#, ln_phi = NA, ln_spvar = NA),
                              se.est = c(b0 = NA),#, ln_phi = NA, ln_spvar = NA),
                              time = as.numeric(difftime(b,a, units = 'mins')),
                              stan.time = NA, meanESS = NA, minESS = NA)
    }
  } else {
    if(mod == 'logistic'){
      results$tmbstan <- list(par.est = exp(summary(tmbstan.mod)[[1]][1:4,1]),
                              se.est = summary(tmbstan.mod)[[1]][1:4,2], 
                              time = difftime(b,a, units = 'mins'), 
                              stan.time = as.numeric(difftime(b,a, units = 'mins')), 
                              meanESS = mean(mon$Bulk_ESS), minESS = min(mon$Bulk_ESS))
      warning('se in log space, not comparable to other model runs')
    }
    if(mod == 'gompertz'){
      results$tmbstan <- list(par.est = c(summary(tmbstan.mod)[[1]][1:2,1],exp(summary(tmbstan.mod)[[1]][3:4,1]) ),
                              se.est = summary(tmbstan.mod)[[1]][1:4,2],
                              time = as.numeric(difftime(b,a, units = 'mins')), 
                              meanESS = mean(mon$Bulk_ESS), minESS = min(mon$Bulk_ESS))
      warning('sigma and tau SE in log space, not comparable to other model runs')
    }
    if(mod == 'spatial'){
      results$tmbstan <- list(par.est = c(summary(tmbstan.mod)[[1]][1,1] ),
                              se.est = summary(tmbstan.mod)[[1]][1,2],
                              omega.est = c(summary(tmbstan.mod)[[1]][2:(nrow(mon)-1)]),
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
  if(mod == 'logistic' | mod == 'gompertz'){
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
    Par <- get(paste0(mod,'.init'))
  }
  if(mod == 'spatial'){
    inits <- mkSpatialInits(dat,1,'stan')
    
    Dat <- inits$Dat
    Par <- get('spatial.inits')
  }
  file <- file.path('src', 'stan', paste0(mod, '.stan'))
  model <- cmdstan_model(file) #compile stan model into C++
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
  if(mod == 'spatial') par.names <- c('b0', 'omega')#, 'kappa', 'tau')
  results <- list(par.est = as.vector(as.matrix(fit$summary(par.names, 'mean')[,2])),
                  se.est = as.vector(as.matrix(fit$summary(par.names, 'sd')[,2])),
                  time = as.numeric(difftime(b,a, units = 'mins')),
                  stan.time = as.numeric(difftime(b,a, units = 'mins')),
                  meanESS =  mean(fit$summary()$ess_bulk, na.rm=TRUE),
                  minESS =  min(fit$summary()$ess_bulk, na.rm=TRUE))
  
  return(results)
}
setupADMB <- function(modName, doRE){ 
  wd <- getwd()
  setwd('src/admb')
  compile_admb(modName, re = doRE, verbose = TRUE)
  setwd(wd)

}
runADMB <- function(dat,mod){
  parm.rep <- function(f){
    rep.out <- strsplit(f,"=")
    parm.nm <- rep.out[[1]][1] 
    parm.val <- as.numeric(strsplit(rep.out[[1]][2], " ")[[1]]) 
    if(is.na(parm.val[1])) parm.val <- parm.val[-1]
    out <- list()
    out[[parm.nm]] = parm.val
    return(out)
  }
  wd <- getwd()
  setwd('src/admb')
  if(mod == 'gompertz'){
    results <- list()
  }
  if(mod == 'logistic'){
    models <- c('logisticGrowth', 'logisticGrowth_nonsep')
    results <- list()
    for(m in 1:length(models)){
      mod.name <- models[m]
      Dat <- list(n=length(dat), y = dat)
      Par <- list(ln_r = log(0.5), ln_K = log(80), ln_sig=-1,ln_tau=-1,  u = rep(1,n) )
      write_dat(mod.name, Dat)
      write_pin(mod.name, Par)
      a <- Sys.time()
      admb.mod <- run_admb(mod.name, verbose = TRUE, extra.args = '-noinit')
      b <- Sys.time()
      rep.name <- paste0(mod.name, '.rep')
      if(rep.name %in% list.files()){
        admb.rep <-  readLines(rep.name)
        parm.out <- sapply(1:length(admb.rep), function(x) parm.rep(admb.rep[x]))
        results[[mod.name]] <- list(par.est = unlist(parm.out),
                        se.est = rep(NA,4),
                        time = as.numeric(difftime(b,a, units = 'mins')),
                        stan.time = as.numeric(difftime(b,a, units = 'mins')),
                        meanESS =  NA,
                        minESS =  NA)
        file.remove(rep.name)
      } else {
        results[[mod.name]] <- list(par.est = rep(NA,4),
                        se.est = rep(NA,4),
                        time = as.numeric(difftime(b,a, units = 'mins')),
                        stan.time = as.numeric(difftime(b,a, units = 'mins')),
                        meanESS =  NA,
                        minESS =  NA)
      }
    }
   
  }
  setwd(wd)
  return(results)
  
}
