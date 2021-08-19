## TMB specific functions


#' Setup C++ TMB code 
#'
#' @param comp boolean that determines whether or not to compile the C++ file
#'
#' @return DLL is compiled and loaded
#' @export
setupTMB <- function(dll.name, comp=FALSE){
  if(comp==TRUE){
    try(dynlib(dyn.unload(paste0('src/tmb/', dll.name))))
    TMB::compile(paste0('src/tmb/', dll.name))
  }
  suppressWarnings(dynlib(dyn.load(paste0('src/tmb/', dll.name))))
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


mkSpatialInits <- function(df){
  Loc <- df[,1:2]
  y <- df[,3]
  
  #build INLA mesh
  mesh <- inla.mesh.2d(Loc, max.edge = c(10,20), offset = c(5,25))
  #calculate sparse distance matrix components
  spde <- inla.spde2.matern(mesh)
  
  Dat <- list(y = y,
              v_i = mesh$idx$loc-1,
              sparse = 1,
              M0 = spde$param.inla$M0,
              M1 = spde$param.inla$M1,
              M2 = spde$param.inla$M2)
  Par <- list(b0 = 0,
              ln_kappa = 0,
              ln_tau = 0,
              omega = rep(0,mesh$n))
  init.list <- list(Dat = Dat, Par = Par)
  return(init.list)
}