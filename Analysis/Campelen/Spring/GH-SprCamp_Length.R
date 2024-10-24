# Set up data
source("./Analysis/Campelen/Spring/SprCamp-SetUp.R")

# Set up model env
n <- c(length(trawlf$year),
       length(gh1$year),length(gh2$year))
pa <- c(trawlf$pa, 
        gh1$pa, gh2$pa)
year <- c(trawlf$year, 
          gh1$year, gh2$year)
names <- c("Trawl", "gh1","gh2")
id <- c(seq(0,(length(names)-1)))


sd_temp <- c(sd(trawlf$bot.temp),
             sd(gh1$bot.temp),sd(gh2$bot.temp))
topt <- c(median(trawlf$bot.temp), 
          median(gh1$bot.temp), median(gh2$bot.temp))
temp <- c(trawlf$bot.temp, 
          gh1$bot.temp, gh2$bot.temp)


spp_n <- c(length(trawlf$year), 
           sum(length(gh1$year),length(gh2$year)))
spp_id <- c(0,1)

model <- 0

gh_fishy_dat <- data.frame(pa = pa, 
                        year = year, 
                        names = rep(names, n),
                        idmod = rep(model, sum(n)),
                        idex = rep(id, n),
                        temp = temp,
                        isd_temp = rep(sd_temp, n),
                        itopt = rep(topt, n))

omega <- data.frame(species = gh_fishy_dat$names, value = rep(NA, length(gh_fishy_dat$pa)))
for(i in 1:length(gh_fishy_dat$pa)){
  omega[i,2] = (1/(gh_fishy_dat$isd_temp[i]*sqrt(2*pi)))*exp((-0.5)*((gh_fishy_dat$temp[i]-gh_fishy_dat$itopt[i])/gh_fishy_dat$isd_temp[i])^2)
}
omega <- omega %>% group_by(species) %>% mutate(min = min(value), max = max(value))

gh_tmb_data <- list(n = length(gh_fishy_dat$pa),
                    nyrs = length(unique(gh_fishy_dat$year)),
                    iyear = gh_fishy_dat$year - min(gh_fishy_dat$year),
                    k = 1,
                    pa = gh_fishy_dat$pa,
                    idmod = gh_fishy_dat$idmod,
                    ndex = length(unique(gh_fishy_dat$idex)),
                    idex = gh_fishy_dat$idex,
                    omega = omega$value,
                    Omin = unique(omega$min)[2:length(n)],
                    Omax = unique(omega$max)[2:length(n)]
)

gh_param_list <- list(
  iye = seq(from = log(10), to = log(20), length.out = gh_tmb_data$nyrs),
  logrw_var = log(1),
  lchi = c(rep(log(0.5), (length(n)-1))),
  lbeta = c(rep(log(0.5), (length(n)-1))),
  lscl = c(rep(log(0.75), (length(n)-1)))
)





map <- list(lscl = as.factor(c(rep(NA, length(gh_param_list$lscl)))))

obj <- MakeADFun(data = gh_tmb_data, map=map,
                 parameters = gh_param_list,
                 DLL = "SclNLFPM", 
                 random = c("iye"))

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(trace = 10, eval.max = 2000, iter.max = 1000), silent = TRUE)


gh_rep <- obj$report()
gh_sdrep <- sdreport(obj)
# gh_sdrep
# 2*opt$objective + 2*length(opt$par)
