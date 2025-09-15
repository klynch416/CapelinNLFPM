# Set up data
load(file = "./ModelSaves/Set_Up/SC-2019.RData")



## Set up for model
n <- c(length(trawls$year), 
       length(c(ac1$year, ac2$year)), 
       length(c(gh1$year, gh2$year)), 
       length(c(ap1$year)))
pa <- c(trawls$pa, 
        c(ac1$pa, ac2$pa), 
        c(gh1$pa, gh2$pa), 
        c(ap1$pa))
year <- c(trawls$year, 
          c(ac1$year, ac2$year), 
          c(gh1$year, gh2$year), 
          c(ap1$year))
names <- c("Trawl", "Atlantic cod", "Greenland halibut", "American plaice")
id = c(0,1,2,3)
model <- 0

sd_temp <- c(sd(trawls$bot.temp), 
             sd(c(ac1$bot.temp, ac2$bot.temp)), 
             sd(c(gh1$bot.temp, gh2$bot.temp)), 
             sd(c(ap1$bot.temp)))
topt <- c(median(trawls$bot.temp), 
          median(c(ac1$bot.temp, ac2$bot.temp)), 
          median(c(gh1$bot.temp, gh2$bot.temp)), 
          median(c(ap1$bot.temp)))
temp <- c(trawls$bot.temp, 
          c(ac1$bot.temp, ac2$bot.temp), 
          c(gh1$bot.temp, gh2$bot.temp), 
          c(ap1$bot.temp))

fishy_dat <- data.frame(pa = pa, 
                        year = year, 
                        names = rep(names, n),
                        idmod = rep(model, sum(n)),
                        idex = rep(id, n), 
                        temp = temp,
                        isd_temp = rep(sd_temp, n),
                        itopt = rep(topt, n))



omega <- data.frame(species = fishy_dat$names, value = rep(NA, length(fishy_dat$pa)))
for(i in 1:length(fishy_dat$pa)){
  omega[i,2] = (1/(fishy_dat$isd_temp[i]*sqrt(2*pi)))*exp((-0.5)*((fishy_dat$temp[i]-fishy_dat$itopt[i])/fishy_dat$isd_temp[i])^2)
}
omega <- omega %>% group_by(species) %>% mutate(min = min(value), max = max(value))


tmb_data <- list(n = length(fishy_dat$pa),
                 nyrs = length(unique(fishy_dat$year)),
                 iyear = fishy_dat$year - min(fishy_dat$year),
                 k = 1,
                 pa = fishy_dat$pa,
                 idmod = fishy_dat$idmod,
                 ndex = length(unique(fishy_dat$idex)),
                 idex = fishy_dat$idex,
                 omega = omega$value,
                 Omin = unique(omega$min)[2:length(n)],
                 Omax = unique(omega$max)[2:length(n)]
)

param_list <- list(
  iye = seq(from = log(10), to = log(20), length.out = tmb_data$nyrs),
  logrw_var = log(1),
  lchi = c(rep(log(0.5), (length(n)-1))),
  lbeta = c(rep(log(0.5), (length(n)-1))),
  lscl = c(rep(log(0.75), (length(n)-1)))
)



map <- list(lscl = as.factor(c(rep(NA, length(param_list$lscl)))))

obj <- MakeADFun(data = tmb_data, map=map,
                 parameters = param_list,
                 DLL = "SclNLFPM", 
                 random = c("iye"))

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(trace = 10, eval.max = 2000, iter.max = 1000), silent = TRUE)


rep <- obj$report()
sdrep <- sdreport(obj)
sdrep
# 2*opt$objective + 2*length(opt$par)
