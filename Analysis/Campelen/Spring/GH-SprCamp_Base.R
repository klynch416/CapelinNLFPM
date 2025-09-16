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
load(file = "./ModelSaves/Set_Up/SC-2019.RData")

gh_fishy_dat <- make_data(n = c(length(trawls$year), 
                                length(c(gh1$year, gh2$year))),
                          pa = c(trawls$pa, 
                                 c(gh1$pa, gh2$pa)), 
                          year = c(trawls$year, 
                                   c(gh1$year, gh2$year)), 
                          names = c("Trawl", "Greenland halibut"),
                          tempType = "noTemp",
                          onto = T
)



# Set up data to be passed into TMB
gh_tmb_data <- model_data(fishy_dat = gh_fishy_dat, modType = "nonlinear", tempType = "noTemp", n = n)



# Set up parameters estimated by TMB
gh_tmb_params <- model_param(tmb_data = gh_tmb_data, lbeta = log(0.5), lchi = log(0.5))



# Run model
obj <- MakeADFun(data = gh_tmb_data,
                 parameters = gh_tmb_params,
                 DLL = "SclNLFPM", 
                 random = c("iye"))

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(trace = 10, eval.max = 2000, iter.max = 1000), silent = TRUE)

# Save model outputs
gh_rep <- obj$report()
gh_sdrep <- sdreport(obj)
gh_sdrep
