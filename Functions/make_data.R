#' Create basic data.frame for the model
#'
#'
#' Combine presence/absence data from bottom trawl surveys, stomach contents, and an index for the different indices and different species/ ontogeny
#'
#' @param pa Vector of presence/absence data from all surveys (binary)
#' @param year Vector of years associated with each observation (numeric)
#' @param n Vector of the number of observations for each survey (integer)
#' @param names Vector of names for each type of data, used for plotting (character)
#' @param tempType Argument for the implementation of temperature: "noTemp", "monoTemp", "scaledTemp" (character) 
#' @param onto Optional, Bollean whether to account for ontogenetic diet shifts (True/False), defaults to False
#' @param temp Optional, Vector of individual temperature observations (vector)
#' @param spp_id Optional, Vector of species identifiers, only needed for \code{tempType = "monoTemp"} when \code{onto = T} (integer)
#' @param spp_n Optional, Vector of the number of observations for each survey and species, only needed for \code{tempType = "monoTemp"} when \code{onto = T} (integer)
#'
#' @return data.frame
#' @export
#'
#' @examples
#' make_data(pa=rbinom(n=20, size=1, prob=0.5), year=rep(seq(from=1, to=10, by=1),2), n=c(10,10), id =c(0,1,1), names=c("Trawl", "Call", "Full"))
make_data <- function(n, pa, year, names, tempType, onto = F, temp = NULL, spp_id = NULL, spp_n = NULL){
  
  if(tempType == "noTemp" | tempType == "scaledTemp" | tempType == "monoTemp"){
    
    id <- c(seq(0,(length(names)-1)))
    
    #############################
    # For no temperature effect #
    #############################
    
    if(onto == T | onto == F){  
      if(tempType == "noTemp"){
        
        fishy_dat <- data.frame(pa = pa, 
                                year = year, 
                                names = rep(names, n),
                                idex = rep(id, n)
        )
      }
    }
    
    
    #############################
    # For no ontogenetic effect #
    #############################
    
    if(onto == F){
      
      # For monotonic temperature shape
      if(tempType == "monoTemp"){
        
        if(is.null(spp_id) == T | is.null(spp_n) == T){
          stop("Please ensure spp_id and spp_n are defined.")
        } else{
          if(is.null(temp == T)){
            stop("Please ensure temp is defined.")
          } else{
            
            fishy_dat <- data.frame(pa = pa, 
                                    year = year, 
                                    names = rep(names, n), 
                                    idex = rep(id, n), 
                                    sppex =  rep(id, n),
                                    temp = temp
            )
          }
        }
      }
      
      # For bell curve temperature shape
      if(tempType == "scaledTemp"){
        
        if(is.null(temp == T)){
          stop("Please ensure temp is defined.")
        }else{
          fishy_dat <- data.frame(pa = pa, 
                                  year = year, 
                                  names = rep(names, n)
          )
          
          temp_eqn <- data.frame(idex = rep(id, n), temp = temp) %>% group_by(idex) %>% mutate(isd_temp = sd(temp), itopt = median(temp))
          
          fishy_dat <- cbind(fishy_dat, temp_eqn)
        }
        
      }
    }
    
    
    ##########################
    # For ontogenetic effect #
    ##########################
    
    if(onto == T){
      
      # For monotonic temperature shape
      if(tempType == "monoTemp"){
        
        if(is.null(spp_id) == T | is.null(spp_n) == T){
          stop("Please ensure spp_id and spp_n are defined.")
        } else{
          if(is.null(temp == T)){
            stop("Please ensure temp is defined.")
          } else{
            
            fishy_dat <- data.frame(pa = pa, 
                                    year = year, 
                                    names = rep(names, n),
                                    idex = rep(id, n),
                                    sppex = rep(spp_id, spp_n),
                                    temp = temp
            )
          }
        }
      }
      
      # For bell curve temperature shape
      if(tempType == "scaledTemp"){
        if(is.null(temp == T)){
          stop("Please ensure temp is defined.")
        }else{
          fishy_dat <- data.frame(pa = pa, 
                                  year = year, 
                                  names = rep(names, n)
          )
          
          temp_eqn <- data.frame(idex = rep(id, n), temp = temp) %>% group_by(idex) %>% mutate(isd_temp = sd(temp), itopt = median(temp))
          
          fishy_dat <- cbind(fishy_dat, temp_eqn)
        }
        
      }
    }
    
    return(fishy_dat)}
  
  else{
    stop("Please select temperature effect: noTemp, scaledTemp, or monoTemp.")
  }
  
}
