# Set up data
source("./Analysis/Campelen/Fall/FallCamp-SetUp.R")




model <- 5

fishy_dat <- data.frame(pa = pa, 
                        year = year, 
                        names = rep(names, n),
                        idmod = rep(model, sum(n)),
                        idex = rep(id, n),
                        sppex = rep(spp_id, spp_n),
                        temp = temp)

tmb_data <- list(n = length(fishy_dat$pa),
                 nyrs = length(unique(fishy_dat$year)),
                 iyear = fishy_dat$year - min(fishy_dat$year),
                 k = 1,
                 pa = fishy_dat$pa,
                 idmod = fishy_dat$idmod,
                 ndex = length(unique(fishy_dat$idex)),
                 idex = fishy_dat$idex,
                 sppex = fishy_dat$sppex,
                 temp = (fishy_dat$temp-min(fishy_dat$temp))
)

param_list <- list(
  iye = seq(from = log(10), to = log(20), length.out = tmb_data$nyrs),
  logrw_var = log(1),
  lchi = c(rep(log(0.5), (length(n)-1))),
  lbeta = c(rep(log(0.5), (length(n)-1))),
  lmink = c(rep(log(1), (length(spp_n)-1))),
  slo = c(rep(2, (length(spp_n)-1))),
  lmidp = c(rep(log(3), (length(spp_n)-1)))
)


obj <- MakeADFun(data = tmb_data,
                 parameters = param_list,
                 DLL = "SppmonoTempNLFPM", 
                 random = c("iye"))

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(trace = 10, eval.max = 2000, iter.max = 1000), silent = TRUE)


rep <- obj$report()
sdrep <- sdreport(obj)
sdrep
2*opt$objective + 2*length(opt$par)
