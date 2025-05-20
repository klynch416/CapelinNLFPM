#' Create parameter list for for general functional response model
#'
#'
#' @param tmb_data TMB data object from the model_data() function
#' @param lbeta A logged starting value for the beta parameter
#' @param lchi A logged starting value for the chi parameter
#' @param lk A placeholder parameter for the linear model that allows it to run, since otherwise there wouldn't be any fixed effects parameters. Don't change from default
#' @param lscl  A logged starting value for the scl (temperature scaling) parameter
#' @param lmink A logged starting value for the mink (minimum consumption rate) parameter
#' @param lmidp A logged starting value for the midp (mid point) parameter
#'
#' @return A list of parameters for the functional response model, including a vector the length of the number of years (iye), 
#' and a starting value for the beta, chi, and random walk variability parameters
#' @export
#'
#' @examples model_param(tmb_data=tmb_data, lbeta=log(2), lchi=log(1), type="nonlinear")
model_param <- function(tmb_data, lbeta, lchi, lscl = NULL, lmink = NULL, lmidp = NULL){
  
  if(tmb_data$modType == "linear"){
    tmb_params <- list(
      iye = seq(from=log(10), to=log(20), length.out=tmb_data$nyrs),
      lk = 0
    )
  }
  
  if(tmb_data$modType == "nonlinear"){
    
    
    if(tmb_data$tempType == "noTemp"){
      tmb_params <- list(
        iye = seq(from=log(10), to=log(20), length.out=tmb_data$nyrs),
        logrw_var = log(1),
        lbeta = c(rep(lbeta, tmb_data$ndex-1)),
        lchi = c(rep(lchi, tmb_data$ndex-1))
      )
    }
    
    
    
    
    if(tmb_data$tempType == "scaledTemp"){
      if(is.null(lscl) == T){
        warning("Automatic value lscl = log(1). Please assign lscl.")
        
        tmb_params <- list(
          iye = seq(from=log(10), to=log(20), length.out=tmb_data$nyrs),
          logrw_var = log(1),
          lbeta = c(rep(lbeta, tmb_data$ndex-1)),
          lchi = c(rep(lchi, tmb_data$ndex-1)),
          lscl = c(rep(lscl, tmb_data$ndex-1))
        )
      } else{
        tmb_params <- list(
          iye = seq(from=log(10), to=log(20), length.out=tmb_data$nyrs),
          logrw_var = log(1),
          lbeta = c(rep(lbeta, tmb_data$ndex-1)),
          lchi = c(rep(lchi, tmb_data$ndex-1)),
          lscl = c(rep(log(1), tmb_data$ndex-1))
        )
      }
      
    }
    
    
    
    
    if(tmb_data$tempType == "monoTemp"){
      if(is.null(lmink) == T){
        Warning("Automatic value lmink = log(1). Please assign lscl.")
        
        tmb_params <- list(
          iye = seq(from=log(10), to=log(20), length.out=tmb_data$nyrs),
          logrw_var = log(1),
          lbeta = c(rep(lbeta, tmb_data$ndex-1)),
          lchi = c(rep(lchi, tmb_data$ndex-1)),
          # parameters for theta temperature curve estimation
          lmink = c(rep(log(1), length(unique(tmb_data$sppex))-1)),
          slo = c(rep(1, (length(tmb_data$spp_n)-1))),
          lmidp = c(rep(lmidp, length(unique(tmb_data$sppex))-1))
        )
        
      }
      
      
      if(is.null(lmidp) == T){
        Warning("Automatic value lmidp = log(15). Please assign lscl value.")
        
        tmb_params <- list(
          iye = seq(from=log(10), to=log(20), length.out=tmb_data$nyrs),
          logrw_var = log(1),
          lbeta = c(rep(lbeta, tmb_data$ndex-1)),
          lchi = c(rep(lchi, tmb_data$ndex-1)),
          # parameters for theta temperature curve estimation
          lmink = c(rep(lmink, length(unique(tmb_data$sppex))-1)),
          slo = c(rep(1, (length(tmb_data$spp_n)-1))),
          lmidp = c(rep(log(15), length(unique(tmb_data$sppex))-1))
        )
      }
      
      if(is.null(lmink) == T & is.null(lmidp) == T){
        Warning("Automatic value lmink = log(1) and lmidp = log(15). Please assign lmink and lminp value.")
        
        tmb_params <- list(
          iye = seq(from=log(10), to=log(20), length.out=tmb_data$nyrs),
          logrw_var = log(1),
          lbeta = c(rep(lbeta, tmb_data$ndex-1)),
          lchi = c(rep(lchi, tmb_data$ndex-1)),
          # parameters for theta temperature curve estimation
          lmink = c(rep(log(1), length(unique(tmb_data$sppex))-1)),
          slo = c(rep(1, length(unique(tmb_data$sppex))-1)),
          lmidp = c(rep(log(15), length(unique(tmb_data$sppex))-1))
        )
        
      } else{
        tmb_params <- list(
          iye = seq(from=log(10), to=log(20), length.out=tmb_data$nyrs),
          logrw_var = log(1),
          lbeta = c(rep(lbeta, tmb_data$ndex-1)),
          lchi = c(rep(lchi, tmb_data$ndex-1)),
          # parameters for theta temperature curve estimation
          lmink = c(rep(lmink, length(unique(tmb_data$sppex))-1)),
          slo = c(rep(1, length(unique(tmb_data$sppex))-1)),
          lmidp = c(rep(lmidp, length(unique(tmb_data$sppex))-1))
        )
      }
      
    }
    
    
    

  }
  
  return(tmb_params)
  
}
