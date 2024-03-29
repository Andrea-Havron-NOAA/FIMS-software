## TMB specific functions


#' Setup C++ TMB code 
#'
#' @param comp boolean that determines whether or not to compile the C++ file
#'
#' @return DLL is compiled and loaded
#' @export
setupTMB <- function(dll.name, comp=FALSE){
  if(!(paste0(dll.name, '.dll') %in% list.files('src/tmb')) |
      comp==TRUE){
    try(dyn.unload(dynlib(paste0('src/tmb/', dll.name))))
    TMB::compile(paste0('src/tmb/', dll.name,'.cpp'))
  }
  suppressWarnings(dyn.load(dynlib(paste0('src/tmb/', dll.name))))
}

#' Make TMB data list 
#'
#' @param obs vector of observation data
#' @param modname string specifying 'gompertz' or 'logistic'
#'
#' @return data list
#' @export
#'
mkTMBdat <- function(obs, modname, prType){
  dat <- list(y=obs)
  if(modname=='gompertz'){
    dat$mod <- 0
  }
  if(modname=='logistic'){
    dat$mod <- 1
  }
  if(prType == 0){
    dat$hyperpars <- 0
  } else {
    dat$hyperpars <- c(-1,4,5,4,0.1,0.1) 
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
                        hyperSig = c(0.1),
                        hyperTau = c(0.1),
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


mkSpatialInits <- function(df,pr,method){
  Loc <- df[,1:2]
  y <- df[,3]
  
  #build INLA mesh
  mesh <- inla.mesh.2d(Loc, max.edge = c(10,20), offset = c(5,25))
  #calculate sparse distance matrix components
  spde <- inla.spde2.matern(mesh)
  
  if(method == 'tmb'){
    Dat <- list(y = y,
                v_i = mesh$idx$loc-1,
                M0 = spde$param.inla$M0,
                M1 = spde$param.inla$M1,
                M2 = spde$param.inla$M2,
                #kap_tau_pr_mu = c(0,0),
                #kap_tau_pr_var = c(0,0),
                prior_type = 0,
                kappa = sqrt(8)/50,
                tau = 1/(sqrt(4*pi*0.75)*sqrt(8)/50))
    if(pr == 1){
      #Dat$lnkapPr <- c(max(dist(Loc))*.2, 10)
      #Dat$lntauPr <- c(0,2)
      # Dat$kap_tau_pr_mu <- c(log(max(dist(Loc))*.2), 0)
      # Dat$kap_tau_pr_var <- c(10, 2)
      Dat$prior_type = 1
    }
    Par.fn <- function(){
      list(b0 = 0,#ln_phi = 0, ln_spvar = 0,
           omega = rep(0,mesh$n))
    }
  } 
  if(method == 'stan'){
    Dat <- list(N = length(y),
                NV = mesh$n,
                y = y,
                vi = mesh$idx$loc,
                M0 = as.matrix(spde$param.inla$M0),
                M1 = as.matrix(spde$param.inla$M1),
                M2 = as.matrix(spde$param.inla$M2),
                #lnkapPr = c(0,0),
                #lntauPr = c(0,0),
                # kap_tau_pr_mu = c(0,0),
                # kap_tau_pr_var = c(0,0),
                prior_type = 0,
                kappa = sqrt(8)/50,
                tau = 1/(sqrt(4*pi*0.75)*sqrt(8)/50))
    if(pr == 1){
      #Dat$lnkapPr <- c(max(dist(Loc))*.2, 10)
      #Dat$lntauPr <- c(0,2)
      #Dat$kap_tau_pr_mu <- c(log(max(dist(Loc))*.2), 0)
      #Dat$kap_tau_pr_var <- c(10, 2)
      Dat$prior_type = 1
    }
  }
  
  init.list <- list(Dat = Dat, Par = Par.fn)
  return(init.list)
}
