# Set up data
library(tidyverse)
library(TMB)
dyn.load("./TMBcode/SclNLFPM")

# load rstrap set data
load(file = "./Data/setdet.rda") 
setdet <- setdet %>% mutate(v.t.s = paste(vessel,trip,set, sep = ".")) %>% 
  filter(spec == "889"|spec == "438"|spec == "892") %>%
  mutate(lat_dd = as.numeric(format(trimws(lat.start), digits = 5)),
         long_dd = as.numeric(format(trimws(long.start), digits = 5)))


# load and format trawl and stomach contents data
camp_fall <- read.csv("./Data/dfo_camp_dat.csv") %>% mutate(v.t.s = paste(VESSEL,TRIP,SET, sep = ".")) %>%
  rename(year = YEAR, trawl_pa = Capelin_PA)
camp_fall <- left_join(camp_fall, 
                       setdet %>% select(v.t.s, set.depth.max, data.series, bot.temp), by = c("v.t.s"), multiple = "any")
camp_fall <- camp_fall %>% filter(data.series == "Campelen" & season == "FALL"& is.na(bot.temp) == FALSE) %>%
  mutate(lat_dd = as.numeric(trimws(lat_dd)),
         long_dd = as.numeric(trimws(long_dd)))


call_fall <- read.csv("./Data/dfo_call_dat.csv")
call_fall <- left_join(call_fall, 
                       camp_fall %>% select(v.t.s, lat_dd, long_dd, set.depth.max, bot.temp), by = c("v.t.s"))
call_fall <- call_fall %>% filter(data.series == "Campelen" & season == "fall" & is.na(bot.temp) == FALSE)%>% 
  mutate(content = rep("called", length(data.series))) %>%
  select(year, month, content, alt.name, pa, NAFOdiv, length, bot.temp, lat_dd, long_dd)
call_fall <- call_fall  %>%
  filter((alt.name == "Atlantic cod" & length > 17)|(alt.name == "Greenland halibut" & length > 19)|(alt.name == "American plaice" & length > 29))


full_fall <- read.csv("./Data/NAFC_diet_capelin_COD_TURBOT_PLAICE_2J3KL.csv") %>% mutate(v.t.s = paste(VESSEL,TRIP,SET, sep = "."))
full_fall <- left_join(full_fall, 
                       camp_fall %>% select(v.t.s, data.series, set.depth.max, bot.temp), by = c("v.t.s"))
full_fall <- full_fall %>% mutate(alt.name = PRED_COMM_NAME,
                                  alt.name = ifelse(alt.name == "AMERICAN PLAICE", "American plaice", alt.name),
                                  alt.name = ifelse(alt.name == "COD,ATLANTIC", "Atlantic cod", alt.name),
                                  alt.name = ifelse(alt.name == "TURBOT", "Greenland halibut", alt.name),
                                  pa = ifelse(PREY_CAP == "Capelin", 1, 0),
                                  content = rep("full", length(PRED_COMM_NAME))) %>% 
  rename(year = YEAR, month = MONTH, length = LENGTH, NAFOdiv = NAFO) %>%
  filter(data.series == "Campelen" & SEASON == "Fall" & is.na(length) == FALSE & length < 400 & is.na(bot.temp) == FALSE) %>%
  select(year, month, alt.name, content, pa, NAFOdiv, length, bot.temp, lat_dd, long_dd)
full_fall <- full_fall %>% 
  filter((alt.name == "Atlantic cod" & length > 17)|(alt.name == "Greenland halibut" & length > 19)|(alt.name == "American plaice" & length > 29))

# Merge species called and full stomach contents
stom_fall <- rbind(call_fall, full_fall)



#Remove data with no capelin
# RasterLayer with change degree parameters
x <- raster::raster(xmn = -61, xmx = -42.5, ymn = 46, ymx = 56)
# Change resolution
terra::res(x) = c(0.5,0.5)
#create spatial dataframe of positive catches
cords <- data.frame(long = -camp_fall$long_dd, lat = camp_fall$lat_dd) 
anomaly_spdf <- sp::SpatialPointsDataFrame(cords, data.frame(pa = camp_fall$trawl_pa))

# Mean catch
mean_catch <- raster::dropLayer(raster::rasterize(anomaly_spdf, x, fun = mean, na.rm = TRUE), 1) 
mean_catch <- terra::rast(mean_catch)
# Get raster values
rasValue <- raster::extract(mean_catch, cords)

# Join lat & long values 
combinePointValue <- cbind(cords, rasValue)
# Filter good raster values
combinePointValue <- combinePointValue %>% filter(pa > 0) ##### modify this

good_coords <- data.frame("lat_dd" = combinePointValue$lat, "long_dd" = combinePointValue$long)
# Remove coords in original df with bad pa levels
camp_fall_good <- camp_fall %>% filter(lat_dd %in% good_coords$lat_dd & long_dd %in% -good_coords$long_dd)
stom_fall_good <- stom_fall %>% filter(lat_dd %in% good_coords$lat_dd & long_dd %in% -good_coords$long_dd)


# Separate data by species
trawlf <- camp_fall_good %>% rename(pa = trawl_pa) %>% select(year, pa, bot.temp)
AC_stof <- stom_fall_good %>% filter(alt.name == "Atlantic cod") %>% select(year, pa, length, content, bot.temp)
GH_stof <- stom_fall_good %>% filter(alt.name == "Greenland halibut") %>% select(year, pa, length, content, bot.temp)
AP_stof <- stom_fall_good %>% filter(alt.name == "American plaice") %>% select(year, pa, length, content, bot.temp) %>% filter(year != 2014)


# Set up different length bins
ac1 <- AC_stof %>% filter(length > 17 & length <= 45)
ac2 <- AC_stof %>% filter(length > 45) %>% filter(year != 1995)
ap1 <- AP_stof %>% filter(length > 29) 
gh1 <- GH_stof %>% filter(length > 19 & length <= 40)
gh2 <- GH_stof %>% filter(length > 40)



# Set up model env
n <- c(length(trawlf$year),
       length(ac1$year),length(ac2$year),
       length(gh1$year),length(gh2$year),
       length(ap1$year))
pa <- c(trawlf$pa, 
        ac1$pa, ac2$pa, 
        gh1$pa, gh2$pa,
        ap1$pa)
year <- c(trawlf$year, 
          ac1$year, ac2$year,
          gh1$year, gh2$year,
          ap1$year)
names <- c("Trawl", "ac1","ac2", "gh1","gh2", "ap1")
id <- c(seq(0,(length(names)-1)))


sd_temp <- c(sd(trawlf$bot.temp),
             sd(ac1$bot.temp),sd(ac2$bot.temp),
             sd(gh1$bot.temp),sd(gh2$bot.temp),
             sd(ap1$bot.temp))
topt <- c(median(trawlf$bot.temp), 
          median(ac1$bot.temp), median(ac2$bot.temp),  
          median(gh1$bot.temp), median(gh2$bot.temp),  
          median(ap1$bot.temp))
temp <- c(trawlf$bot.temp, 
          ac1$bot.temp, ac2$bot.temp,
          gh1$bot.temp, gh2$bot.temp,
          ap1$bot.temp)


spp_n <- c(length(trawlf$year), 
           sum(length(ac1$year),length(ac2$year)), 
           sum(length(gh1$year),length(gh2$year)),
           sum(length(ap1$year)))
spp_id <- c(0,1,2,3)








# 
# table(ac1$pa, ac1$year)
# table(ac2$pa, ac2$year)
# table(gh1$pa, gh1$year)
# table(gh2$pa, gh2$year)
# table(ap1$pa, ap1$year)
# 
# table(ac1$content)
# table(ac2$content)
# table(gh1$content)
# table(gh2$content)
# table(ap1$content)
# 
# 
# # %empty stomachs - have to change
# length(filter(ac1, pa == 0)$pa)/length(ac1$pa)*100
# length(filter(ac2, pa == 0)$pa)/length(ac2$pa)*100
# length(filter(gh1, pa == 0)$pa)/length(gh1$pa)*100
# length(filter(gh2, pa == 0)$pa)/length(gh2$pa)*100
# length(filter(ap1, pa == 0)$pa)/length(ap1$pa)*100
# 
# 
# 
# table(AC_stof$content)
# table(GH_stof$content)
# table(AP_stof$content)
# 
# 
# length(filter(AC_stof, pa == 0)$pa)/length(AC_stof$pa)*100
# length(filter(GH_stof, pa == 0)$pa)/length(GH_stof$pa)*100
# length(filter(AP_stof, pa == 0)$pa)/length(AP_stof$pa)*100
