# load rstrap set data
load(file = "./Data/setdet.rda") 
setdet <- setdet %>% mutate(v.t.s = paste(vessel,trip,set, sep = ".")) %>% 
  filter(spec == "889"|spec == "438"|spec == "892") %>%
  mutate(lat_dd = as.numeric(format(trimws(lat.start), digits = 5)), # check that lat and longs are the same format in all data sources
         long_dd = as.numeric(format(trimws(long.start), digits = 5))) # trimws() needed due to data entry errors


# Trawl data
eng_fall <- read.csv("./Data/camp_capelin_2J3KL_FULL.csv") %>% rename(year = YEAR, trawl_pa = Capelin_PA) %>% 
  mutate(v.t.s = paste(VESSEL,TRIP,SET, sep = "."),
         season = ifelse(MONTH == 3|MONTH == 4|MONTH == 5|MONTH == 6|MONTH == 7, "spring", "fall"),
         lat_dd = as.numeric(format(trimws(lat_dd), digits = 5)), # check that lat and longs are the same format in all data sources
         long_dd = as.numeric(format(trimws(long_dd), digits = 5))) # trimws() needed due to data entry errors

eng_fall <- left_join(eng_fall, 
                      setdet %>% select(v.t.s, set.depth.max, data.series, bot.temp), by = c("v.t.s"), multiple = "any")
eng_fall <- eng_fall %>% filter(data.series == "Engel" & 
                                  season == "fall" & 
                                  is.na(bot.temp) == FALSE &
                                  year >= 1984)


# Called stomach data
load(file = "./Data/ag.rda")

call_fall <- ag %>% select(-c(source.file, oedge:gutremvol)) %>%  
  filter((spec == "889" | spec == "438" | spec == "892") & 
           (NAFOdiv == "2J" | NAFOdiv == "3K" |NAFOdiv == "3L") &
           which.survey == "multispecies") %>% 
  mutate(v.t.s = paste(vessel, trip, set, sep = "."),
         alt.name = spec,
         alt.name = ifelse(alt.name == "889", "American plaice", alt.name),
         alt.name = ifelse(alt.name == "438", "Atlantic cod", alt.name),
         alt.name = ifelse(alt.name == "892", "Greenland halibut", alt.name),
         pa = ifelse(prey1 == "Capelin"|prey2 == "Capelin", 1, 0), 
         content = rep("called", length(data.series)),
         season = ifelse(month == 3|month == 4|month == 5|month == 6|month == 7, "spring", "fall"))
call_fall <- left_join(call_fall, 
                       setdet %>% select(v.t.s, lat_dd, long_dd, set.depth.max, bot.temp), by = c("v.t.s"), multiple = "any")
call_fall <- call_fall %>%
  filter(data.series == "Engel" & 
           season == "fall" & 
           is.na(length) == FALSE & 
           is.na(empty) == FALSE & 
           is.na(bot.temp) == FALSE) %>%
  select(year, month, season, alt.name, content, pa, empty, NAFOdiv, length, bot.temp, lat_dd, long_dd)



# Full stomach data
full_fall <- read.csv("./Data/NAFC_diet_capelin_COD_TURBOT_PLAICE_2J3KL.csv") %>% mutate(v.t.s = paste(VESSEL,TRIP,SET, sep = "."))
full_fall <- left_join(full_fall, 
                       setdet %>% select(v.t.s, set.depth.max, data.series, bot.temp), by = c("v.t.s"), multiple = "any")
full_fall <- full_fall %>% mutate(alt.name = PRED_COMM_NAME,
                                  alt.name = ifelse(alt.name == "AMERICAN PLAICE", "American plaice", alt.name),
                                  alt.name = ifelse(alt.name == "COD,ATLANTIC", "Atlantic cod", alt.name),
                                  alt.name = ifelse(alt.name == "TURBOT", "Greenland halibut", alt.name),
                                  pa = ifelse(PREY_CAP == "Capelin", 1, 0),
                                  empty = ifelse(PREY_CAP == "Empty", TRUE, FALSE),
                                  content = rep("full", length(PRED_COMM_NAME)),
                                  SEASON = ifelse(MONTH == 3|MONTH == 4|MONTH == 5|MONTH == 6|MONTH == 7, "spring", "fall")) %>% 
  rename(year = YEAR, month = MONTH, length = LENGTH, NAFOdiv = NAFO, season = SEASON) %>%
  filter(data.series == "Engel" & 
           season == "fall" & 
           length < 400 & 
           is.na(length) == FALSE & 
           is.na(bot.temp) == FALSE) %>%
  select(year, month, season, alt.name, content, pa, empty, NAFOdiv, length, bot.temp, lat_dd, long_dd)




# Merge species called and full stomach contents
stom_fall <- rbind(call_fall, full_fall) %>% filter(year >= 1984) 



#Remove data with no capelin
# RasterLayer with change degree parameters
x <- raster::raster(xmn = -61, xmx = -42.5, ymn = 46, ymx = 56)
# Change resolution
terra::res(x) = c(0.5,0.5)
#create spatial dataframe of positive catches
cords <- data.frame(long = -eng_fall$long_dd, lat = eng_fall$lat_dd) 
anomaly_spdf <- sp::SpatialPointsDataFrame(cords, data.frame(pa = eng_fall$trawl_pa))

# Mean catch
mean_catch <- raster::dropLayer(raster::rasterize(anomaly_spdf, x, fun = mean, na.rm = TRUE), 1) 
mean_catch <- terra::rast(mean_catch)
# Get raster values
rasValue <- raster::extract(mean_catch, cords)

# Join lat & long values 
combinePointValue <- cbind(cords, rasValue) %>% 
  filter(pa > 0) # filter good raster values, mean pa in cell is be above 0


good_coords <- data.frame("lat_dd" = combinePointValue$lat, "long_dd" = combinePointValue$long)
# Remove coords in original df with bad pa levels
eng_fall_good <- eng_fall %>% filter(lat_dd %in% good_coords$lat_dd & long_dd %in% -good_coords$long_dd)
stom_fall_good <- stom_fall %>% filter(lat_dd %in% good_coords$lat_dd & long_dd %in% -good_coords$long_dd)


# Separate data by species
trawlf <- eng_fall_good %>% rename(pa = trawl_pa) %>% select(year, pa, bot.temp) %>% mutate(geartype = "engel")
AC_stof <- stom_fall_good %>% filter(alt.name == "Atlantic cod") %>% select(year, pa, length, season, empty, content, bot.temp) %>% mutate(geartype = "engel")
GH_stof <- stom_fall_good %>% filter(alt.name == "Greenland halibut") %>% select(year, pa, length, season, empty, content, bot.temp) %>% mutate(geartype = "engel")
AP_stof <- stom_fall_good %>% filter(alt.name == "American plaice") %>% select(year, pa, length, season, empty, content, bot.temp) %>% mutate(geartype = "engel")

# Set up different length bins
ac1 <- AC_stof %>% filter(length > 17 & length <= 45)
ac2 <- AC_stof %>% filter(length > 45) %>% filter(year != 1994) # removed years with fewer than 10 overall observations
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
