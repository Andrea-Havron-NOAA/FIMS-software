source('../../data/simdata.R')
source('../../R/utils.R')
source('../../R/model_setup.R')



run_models <- function(n.seq){
#for(i in 1:length(n.seq)){
  n <- 2^n.seq
 
  #logistic model
  Mod <- 'logistic'
  
  logistic.init <- function(){
    list(theta = c(log(0.5), log(80)), ln_sig=-1,ln_tau=-1,
         u = rep(1,n))
  }
  
  simdata <- gendat(seed=1,
                    N=n,
                    theta = c(0.2,100),
                    u1 = 4,
                    var = list(proc=0.01,obs=0.001),
                    mod.name = Mod)
  #write.csv(data.frame(y=simdata), file = paste0('./data/logistic/logistic', '_n', n, '.csv'))
  results <- runTMB(simdata,Mod) 
  print(n.seq)
  #modify and save for Julia
  inits <- c(r = unname(exp(results$inits[1])), K = unname(exp(results$inits[2])), sigma = unname(exp(results$init[3])), 
             tau = unname(exp(results$inits[4])), u_init = unname(results$inits[5]), 
             results$inits[6:length(results$inits)])
  #save(inits, file = paste0('data/logistic/logisticInits', '_n', n, '.RData'))
  logistic.results <-  list(tmb = results$tmb, tmbstan = results$tmbstan)
  
  # Init functions
  logistic.init <- function(){
    list(theta = results$inits[1:2],
         ln_sig = results$inits[3],
         ln_tau = results$inits[4],
         u = results$inits[5:length(results$inits)])
  }
  #stan
  #use vague priors
  logistic.results$stan <- runSTAN(simdata, Mod,1)
 

}


