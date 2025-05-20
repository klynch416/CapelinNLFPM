#' Create data list for general functional response model
#'
#'
#' @param fishy_dat Object output from the make_data() function
#' @param n A vector the length of the number of datasets used where each number indicates the number of observations for that dataset
#' @param modType A character identifying whether you are creating data for a linear or nonlinear functional response model. Input "linear" for a linear model and "nonlinear" for a nonlinear model
#' @param tempType Argument for the implementation of temperature: "noTemp", "monoTemp", "scaledTemp" (character) 
#'
#'
#'
#' @return A list of data to be used when generating the TMB model. The list includes the total number of observations (n), the number of years (nyrs),
#' the number of datasets (ndex), an index vector for the year of each observation (iyear), an index vector for the dataset associated with each observation (idex), and
#' a vector for all presence/absence data (pa)
#' @export
#'
#' @examples model_data(fishy_dat=fishy_dat, n=c(length(df$pa),length(df1$pa),length(df2$pa)), modType="nonlinear", tempType="noTemp")
model_data <- function(fishy_dat, n, modType, tempType){
  
  if(modType == "linear"){
    tmb_data <- list(n = length(fishy_dat$pa),
                     nyrs = length(unique(fishy_dat$year)),
                     ndex = length(unique(fishy_dat$idex)),
                     iyear = fishy_dat$year-min(fishy_dat$year),
                     idex = fishy_dat$idex,
                     pa = fishy_dat$pa,
                     modType = "linear"
    )
    
  }
  
  
  
  if(modType == "nonlinear"){
    
    if(tempType == "noTemp"){
      tmb_data <- list(n = length(fishy_dat$pa),
                       nyrs = length(unique(fishy_dat$year)),
                       ndex = length(unique(fishy_dat$idex)),
                       iyear = fishy_dat$year-min(fishy_dat$year),
                       idex = fishy_dat$idex,
                       pa = fishy_dat$pa,
                       k = 1, # maximum consumtion rate
                       idmod = 0,
                       tempType = "noTemp",
                       modType = "nonlinear"
      )
      
    }
    
    if(tempType == "scaledTemp"){
      
      # Generate theta values
      theta <- data.frame(species = fishy_dat$names, value = rep(NA, length(fishy_dat$pa)))
      for(i in 1:length(fishy_dat$pa)){
        theta[i,2] = (1/(fishy_dat$isd_temp[i]*sqrt(2*pi)))*exp((-0.5)*((fishy_dat$temp[i]-fishy_dat$itopt[i])/fishy_dat$isd_temp[i])^2)
      }
      theta <- theta %>% group_by(species) %>% mutate(min = min(value), max = max(value))
      
      
      tmb_data <- list(n = length(fishy_dat$pa),
                       nyrs = length(unique(fishy_dat$year)),
                       iyear = fishy_dat$year - min(fishy_dat$year),
                       k = 1, # maximum consumtion rate
                       pa = fishy_dat$pa,
                       idmod = 5,
                       ndex = length(unique(fishy_dat$idex)),
                       idex = fishy_dat$idex,
                       theta = theta$value,
                       Tmin = unique(theta$min)[2:length(n)], # start indexing at 2 as to not include trawl information
                       Tmax = unique(theta$max)[2:length(n)],
                       tempType = "scaledTemp",
                       modType = "nonlinear"
      )
      
    }
    
    
    if(tempType == "monoTemp"){
      tmb_data <- list(n = length(fishy_dat$pa),
                       nyrs = length(unique(fishy_dat$year)),
                       iyear = fishy_dat$year - min(fishy_dat$year),
                       k = 1, # maximum consumtion rate
                       pa = fishy_dat$pa,
                       idmod = 5,
                       ndex = length(unique(fishy_dat$idex)),
                       idex = fishy_dat$idex,
                       sppex = fishy_dat$sppex,
                       temp = (fishy_dat$temp+2), # 2:constant added to ensure temperature is always a positive value
                       tempType = "monoTemp",
                       modType = "nonlinear"
      )
      
    }
    
    
  }
  
  return(tmb_data)
  
}