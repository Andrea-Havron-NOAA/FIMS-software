## TMB specific functions


#' Setup C++ TMB code 
#'
#' @param comp boolean that determines whether or not to compile the C++ file
#'
#' @return DLL is compiled and loaded
#' @export
setupTMB <- function(comp=FALSE){
  if(comp==TRUE){
    try(dynlib(dyn.unload('src/tmb/stateSpace')))
    TMB::compile('src/tmb/stateSpace.cpp')
  }
  suppressWarnings(dynlib(dyn.load('src/tmb/stateSpace')))
}

#' Make TMB data list 
#'
#' @param obs vector of observation data
#' @param modname string specifying 'gompertz' or 'logistic'
#'
#' @return data list
#' @export
#'
mkTMBdat <- function(obs, modname){
  dat <- list(y=obs)
  if(modname=='gompertz'){
    dat$mod <- 0
  }
  if(modname=='logistic'){
    dat$mod <- 1
  }
  return(dat)
}

## Stan specific functions


#' Make STAN data list
#'
#' @param obs vector of observations
#' @param ptype integer specifying if using improper priors (ptype=0) 
#'              or vague proper priors (ptype=1)
#'
#' @return data list
#' @export
mkSTANdat <- function(obs,
                      hyperParms = list(
                        hyperSig = c(0.001,0.001),
                        hyperTau = c(0.001,0.001),
                        hyperTheta1 = c(-1,4),
                        hyperTheta2 = c(5,4)
                      )){
 return(list(N = length(obs), 
                y = obs,
                hyp_sig = hyperParms$hyperSig,
                hyp_tau = hyperParms$hyperTau,
                hyp_theta1 = hyperParms$hyperTheta1,
                hyp_theta2 = hyperParms$hyperTheta2,
                prior_type = 0))
  
}

