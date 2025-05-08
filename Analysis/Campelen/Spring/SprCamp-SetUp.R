# load rstrap set data
load(file = "./Data/setdet.rda") 
setdet <- setdet %>% mutate(v.t.s = paste(vessel,trip,set, sep = ".")) %>% 
  filter(spec == "889"|spec == "438"|spec == "892") %>%
  mutate(lat_dd = as.numeric(format(trimws(lat.start), digits = 5)), # check that lat and longs are the same format in all data sources
         long_dd = as.numeric(format(trimws(long.start), digits = 5))) # trimws() needed due to data entry errors


# Trawl data
camp_spr <- read.csv("./Data/camp_capelin_2J3KL_FULL.csv") %>% rename(year = YEAR, trawl_pa = Capelin_PA) %>% 
  mutate(v.t.s = paste(VESSEL,TRIP,SET, sep = "."),
         season = ifelse(MONTH == 3|MONTH == 4|MONTH == 5|MONTH == 6|MONTH == 7, "spring", "spr"),
         lat_dd = as.numeric(format(trimws(lat_dd), digits = 5)), # check that lat and longs are the same format in all data sources
         long_dd = as.numeric(format(trimws(long_dd), digits = 5))) # trimws() needed due to data entry errors
camp_spr <- left_join(camp_spr, 
                      setdet %>% select(v.t.s, set.depth.max, data.series, bot.temp), by = c("v.t.s"), multiple = "any")
camp_spr <- camp_spr %>% filter(data.series == "Campelen" & 
                                  season == "spring"& 
                                  is.na(bot.temp) == FALSE &
                                  year <= 2020) 


# Called stomach data
load(file = "./Data/ag.rda")

call_spr <- ag %>% select(-c(source.file, oedge:gutremvol)) %>%  
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
call_spr <- left_join(call_spr, 
                       setdet %>% select(v.t.s, lat_dd, long_dd, set.depth.max, bot.temp), by = c("v.t.s"), multiple = "any")
call_spr <- call_spr %>%
  filter(data.series == "Campelen" & 
           season == "spring" & 
           is.na(length) == FALSE & 
           is.na(empty) == FALSE & 
           is.na(bot.temp) == FALSE) %>%
  select(year, month, season, alt.name, content, pa, empty, NAFOdiv, length, bot.temp, lat_dd, long_dd)



# Full stomach data
full_spr <- read.csv("./Data/NAFC_diet_capelin_COD_TURBOT_PLAICE_2J3KL.csv") %>% mutate(v.t.s = paste(VESSEL,TRIP,SET, sep = "."))
full_spr <- left_join(full_spr, 
                       setdet %>% select(v.t.s, set.depth.max, data.series, bot.temp), by = c("v.t.s"), multiple = "any")
full_spr <- full_spr %>% mutate(alt.name = PRED_COMM_NAME,
                                  alt.name = ifelse(alt.name == "AMERICAN PLAICE", "American plaice", alt.name),
                                  alt.name = ifelse(alt.name == "COD,ATLANTIC", "Atlantic cod", alt.name),
                                  alt.name = ifelse(alt.name == "TURBOT", "Greenland halibut", alt.name),
                                  pa = ifelse(PREY_CAP == "Capelin", 1, 0),
                                  empty = ifelse(PREY_CAP == "Empty", TRUE, FALSE),
                                  content = rep("full", length(PRED_COMM_NAME)),
                                  SEASON = ifelse(MONTH == 3|MONTH == 4|MONTH == 5|MONTH == 6|MONTH == 7, "spring", "fall")) %>% 
  rename(year = YEAR, month = MONTH, length = LENGTH, NAFOdiv = NAFO, season = SEASON) %>%
  filter(data.series == "Campelen" & 
           season == "spring" & 
           length < 400 & 
           is.na(length) == FALSE & 
           is.na(bot.temp) == FALSE) %>%
  select(year, month, season, alt.name, content, pa, empty, NAFOdiv, length, bot.temp, lat_dd, long_dd)




# Merge species called and full stomach contents
stom_spr <- rbind(call_spr, full_spr)



#Remove data with no capelin
# RasterLayer with change degree parameters
x <- raster::raster(xmn = -61, xmx = -42.5, ymn = 46, ymx = 56)
# Change resolution
terra::res(x) = c(0.5,0.5)
#create spatial dataframe of positive catches
cords <- data.frame(long = -camp_spr$long_dd, lat = camp_spr$lat_dd) 
anomaly_spdf <- sp::SpatialPointsDataFrame(cords, data.frame(pa = camp_spr$trawl_pa))

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
camp_spr_good <- camp_spr %>% filter(lat_dd %in% good_coords$lat_dd & long_dd %in% -good_coords$long_dd)
stom_spr_good <- stom_spr %>% filter(lat_dd %in% good_coords$lat_dd & long_dd %in% -good_coords$long_dd)


# Separate data by species
trawls <- camp_spr_good %>% rename(pa = trawl_pa) %>% select(year, pa, bot.temp) %>% mutate(geartype = "campelen")
AC_stos <- stom_spr_good %>% filter(alt.name == "Atlantic cod") %>% select(year, pa, season, content, empty, length, bot.temp) %>% mutate(geartype = "campelen")
GH_stos <- stom_spr_good %>% filter(alt.name == "Greenland halibut") %>% select(year, pa, season, content, empty, length, bot.temp) %>% mutate(geartype = "campelen")
AP_stos <- stom_spr_good %>% filter(alt.name == "American plaice") %>% select(year, pa, season, content, empty, length, bot.temp) %>% mutate(geartype = "campelen")


# Set up different length bins
ac1 <- AC_stos %>% filter(length > 17 & length <= 45) %>% filter(year != 1998 & year != 2017) # removed years with fewer than 10 overall observations
ac2 <- AC_stos %>% filter(length > 45) %>% filter(year != 2015)
ap1 <- AP_stos %>% filter(length > 29) %>% filter(year != 2017)
gh1 <- GH_stos %>% filter(length > 19 & length <= 40) %>% filter(year != 2017)
gh2 <- GH_stos %>% filter(length > 40)

# Set up model env
n <- c(length(trawls$year),
       length(ac1$year),length(ac2$year),
       length(gh1$year),length(gh2$year),
       length(ap1$year))
pa <- c(trawls$pa, 
        ac1$pa, ac2$pa, 
        gh1$pa, gh2$pa,
        ap1$pa)
year <- c(trawls$year, 
          ac1$year, ac2$year,
          gh1$year, gh2$year,
          ap1$year)
names <- c("Trawl", "ac1","ac2", "gh1","gh2", "ap1")

