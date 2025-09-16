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

ac_fishy_dat <- make_data(n = c(length(trawls$year), 
                                length(c(ac1$year, ac2$year)), 
                                length(c(ap1$year))),
                          pa = c(trawls$pa, 
                                 c(ac1$pa, ac2$pa), 
                                 c(ap1$pa)), 
                          year = c(trawls$year, 
                                   c(ac1$year, ac2$year), 
                                   c(ap1$year)), 
                          names = c("Trawl", "Atlantic cod", "American plaice"),
                          tempType = "noTemp",
                          onto = T
)



# Set up data to be passed into TMB
ac_tmb_data <- model_data(fishy_dat = ac_fishy_dat, modType = "nonlinear", tempType = "noTemp", n = n)



# Set up parameters estimated by TMB
ac_tmb_params <- model_param(tmb_data = ac_tmb_data, lbeta = log(0.5), lchi = log(0.5))



# Run model
obj <- MakeADFun(data = ac_tmb_data,
                 parameters = ac_tmb_params,
                 DLL = "SclNLFPM", 
                 random = c("iye"))

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(trace = 10, eval.max = 2000, iter.max = 1000), silent = TRUE)

# Save model outputs
ac_rep <- obj$report()
ac_sdrep <- sdreport(obj)
ac_sdrep
