# Set up environment
library(tidyverse)
library(TMB)

# Load functions
source("./Functions/make_data.R")
source("./Functions/model_data.R")
source("./Functions/model_param.R")

# Load TMB models
dyn.load("./TMBcode/SclNLFPM")
dyn.load("./TMBcode/SppmonoTempNLFPM")


# Set up data for model
load(file = "./ModelSaves/Set_Up/FC-2020.RData")

## This structure goes by species, so need to aggregate species data together
fishy_dat <- make_data(n = c(length(trawlf$year), 
                             length(c(ac1$year, ac2$year)), 
                             length(c(gh1$year, gh2$year)), 
                             length(c(ap1$year))), 
                       pa = c(trawlf$pa, 
                              c(ac1$pa, ac2$pa), 
                              c(gh1$pa, gh2$pa), 
                              c(ap1$pa)), 
                       year = c(trawlf$year, 
                                c(ac1$year, ac2$year), 
                                c(gh1$year, gh2$year), 
                                c(ap1$year)), 
                       names = c("Trawl", "Atlantic cod", "Greenland halibut", "American plaice"), 
                       onto = F,
                       tempType = "monoTemp",
                       spp_id = spp_id,
                       spp_n = spp_n,
                       temp = c(trawlf$bot.temp, 
                                c(ac1$bot.temp, ac2$bot.temp),
                                c(gh1$bot.temp, gh2$bot.temp),
                                ap1$bot.temp)
)



# Set up data to be passed into TMB
tmb_data <- model_data(fishy_dat = fishy_dat, modType = "nonlinear", tempType = "monoTemp",
                       n = c(length(trawlf$year), 
                             length(c(ac1$year, ac2$year)), 
                             length(c(gh1$year, gh2$year)), 
                             length(c(ap1$year)))
)



# Set up parameters estimated by TMB
tmb_params <- model_param(tmb_data = tmb_data, lbeta = log(0.5), lchi = log(0.5), lmink = log(1), lmidp = log(15))


# Run model
obj <- MakeADFun(data = tmb_data,
                 parameters = tmb_params,
                 DLL = "SppmonoTempNLFPM", 
                 random = c("iye"))

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(trace = 10, eval.max = 2000, iter.max = 1000), silent = TRUE)

# Save model outputs
rep <- obj$report()
sdrep <- sdreport(obj)
sdrep
2*opt$objective + 2*length(opt$par)
