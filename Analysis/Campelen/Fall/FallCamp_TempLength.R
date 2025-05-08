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

fishy_dat <- make_data(n = n,
                       pa = pa, 
                       year = year, 
                       names = names,
                       tempType = "scaledTemp",
                       onto = T,
                       temp = temp
)




# Set up data to be passed into TMB
tmb_data <-  model_data(fishy_dat = fishy_dat, modType = "nonlinear", tempType = "scaledTemp", n = n)


# Set up parameters estimated by TMB
tmb_params <- model_param(tmb_data = tmb_data, lbeta = log(0.5), lchi = log(0.5), lscl = log(0.75))


# Run model
obj <- MakeADFun(data = tmb_data,
                 parameters = tmb_params,
                 DLL = "SclNLFPM", 
                 random = c("iye"))

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(trace = 10, eval.max = 2000, iter.max = 1000), silent = TRUE)

# Save model outputs
rep <- obj$report()
sdrep <- sdreport(obj)
sdrep
2*opt$objective + 2*length(opt$par)

