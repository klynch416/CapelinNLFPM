# Set up data
source("./Analysis/Campelen/Spring/SprCamp-SetUp.R")

# Set up model env
n <- c(length(trawlf$year),
       length(ac1$year),length(ac2$year),
       length(ap1$year))
pa <- c(trawlf$pa, 
        ac1$pa, ac2$pa, 
        ap1$pa)
year <- c(trawlf$year, 
          ac1$year, ac2$year,
          ap1$year)
names <- c("Trawl", "ac1","ac2", "ap1")
id <- c(seq(0,(length(names)-1)))


sd_temp <- c(sd(trawlf$bot.temp),
             sd(ac1$bot.temp),sd(ac2$bot.temp),
             sd(ap1$bot.temp))
topt <- c(median(trawlf$bot.temp), 
          median(ac1$bot.temp), median(ac2$bot.temp),  
          median(ap1$bot.temp))
temp <- c(trawlf$bot.temp, 
          ac1$bot.temp, ac2$bot.temp,
          ap1$bot.temp)


spp_n <- c(length(trawlf$year), 
           sum(length(ac1$year),length(ac2$year)), 
           sum(length(ap1$year)))
spp_id <- c(0,1,2)

model <- 0

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
