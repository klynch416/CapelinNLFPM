# load rstrap set data
load(file = "./Data/setdet.rda") 
setdet <- setdet %>% mutate(v.t.s = paste(vessel,trip,set, sep = ".")) %>% 
  mutate(lat_dd = as.numeric(format(trimws(lat.start), digits = 5)), # check that lat and longs are the same format in all data sources
         long_dd = as.numeric(format(trimws(long.start), digits = 5))) # trimws() needed due to data entry errors


# Trawl data
camp_fall <- read.csv("./Data/camp_capelin_2J3KL_FULL.csv") %>% 
  mutate(v.t.s = paste(VESSEL,TRIP,SET, sep = "."),
         season = ifelse(MONTH == 3|MONTH == 4|MONTH == 5|MONTH == 6|MONTH == 7, "spring", "fall"),
         lat_dd = as.numeric(format(trimws(lat_dd), digits = 5)), # check that lat and longs are the same format in all data sources
         long_dd = as.numeric(format(trimws(long_dd), digits = 5))) # trimws() needed due to data entry errors
camp_fall <- left_join(camp_fall, 
                       setdet %>% select(v.t.s, set.depth.max, data.series, bot.temp), by = c("v.t.s"), multiple = "any") #checked that it doesnt matter which match it picks
camp_fall <- camp_fall %>% rename(year = YEAR, trawl_pa = Capelin_PA) %>% 
  filter(data.series == "Campelen" &
           season == "fall" & 
           is.na(bot.temp) == FALSE &
           year <= 2020) 


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
  filter(data.series == "Campelen" & 
           season == "fall" & 
           is.na(pa) == FALSE &
           is.na(length) == FALSE & 
           is.na(bot.temp) == FALSE) %>%
  select(year, month, season, alt.name, content, pa, empty, NAFOdiv, length, bot.temp, lat_dd, long_dd)



# Full stomach data
# full_fall <- read.csv("./Data/NAFC_diet_capelin_COD_TURBOT_PLAICE_2J3KL_2025UPDATE_allprey.csv") %>%
#   filter(YEAR >= 1984 & YEAR  <= 2020)
# 
# prey1 <- full_fall %>% group_by(ID_PRED) %>% slice(which.max(MASS))
# #retain only elements in stomach not equal to elements in prey1, rows from df1 which are not in df2
# prey2 <- setdiff(full_fall, prey1) %>% group_by(ID_PRED) %>% 
#   slice(which.max(MASS)) #find new highest prey mass
# 
# uni_full_fall <- full_fall %>% select(!c(TOT_PREY_NUMBER:MASS)) %>% unique() %>% # create unique pred obs dataframe
#   mutate(prey1 = NA, # add prey1 and prey2 columns to fill
#          prey2 = NA)
# 
# # match prey1 and prey2 for each ID_PRED, this step takes a while
# pb <- txtProgressBar(min = 0,max = nrow(uni_full_fall),style = 3,width = 50,char = "=")
# for(i in 1:nrow(uni_full_fall)){  
#   uni_full_fall[i, 'prey1'] <- prey1[prey1$ID_PRED==uni_full_fall$ID_PRED[i], ]$PREY_GR_RAP
#   # not all preds have top prey2. match which ones do, NA for those that dont
#   if(uni_full_fall$ID_PRED[i] %in% unique(prey2$ID_PRED)){
#     uni_full_fall[i, 'prey2'] <- prey2[prey2$ID_PRED==uni_full_fall$ID_PRED[i], ]$PREY_GR_RAP
#   } else(uni_full_fall[i, 'prey2'] <- NA)
#   setTxtProgressBar(pb, i)
# }
# close(pb) # Close the connection

load(file = "./Data/uni_full_stom.RData")

full_fall <- uni_full_stom %>% rename(year = YEAR, month = MONTH, length = LENGTH, NAFOdiv = NAFO, season = SEASON) %>%
  mutate(v.t.s = paste(VESSEL,TRIP,SET, sep = "."),
         pa = ifelse(prey1 == "Capelin"|prey2 == "Capelin", 1, 0)) %>% 
  mutate(pa = ifelse(prey1 != "Capelin" & is.na(prey2) == T, 0, pa)) %>%
  mutate(alt.name = SPECIES,
         alt.name = ifelse(alt.name == "889", "American plaice", alt.name),
         alt.name = ifelse(alt.name == "438", "Atlantic cod", alt.name),
         alt.name = ifelse(alt.name == "892", "Greenland halibut", alt.name),
         content = rep("full", length(PRED_COMM_NAME)),
         empty = ifelse(prey1 == "Empty", 1, 0),
         season = ifelse(month == 3|month == 4|month == 5|month == 6|month == 7, "spring", "fall"))
full_fall <- left_join(full_fall, 
                       setdet %>% select(v.t.s, set.depth.max, data.series, bot.temp), by = c("v.t.s"), multiple = "any") %>%
  filter(data.series == "Campelen" & 
           season == "fall" & 
           length < 200 & 
           is.na(pa) == FALSE &
           is.na(length) == FALSE &
           is.na(bot.temp) == FALSE) %>%
  select(year, month, season, alt.name, content, pa, empty, NAFOdiv, length, bot.temp, lat_dd, long_dd)




# Merge species called and full stomach contents
stom_fall <- rbind(call_fall, full_fall) %>% filter(year <= 2020)



# Remove data with no capelin
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
combinePointValue <- cbind(cords, rasValue) %>% 
  filter(pa > 0) # filter good raster values, mean pa in cell must be above 0


good_coords <- data.frame("lat_dd" = combinePointValue$lat, "long_dd" = combinePointValue$long)
# Remove coords in original df with bad pa levels
camp_fall_good <- camp_fall %>% filter(lat_dd %in% good_coords$lat_dd & long_dd %in% -good_coords$long_dd)
stom_fall_good <- stom_fall %>% filter(lat_dd %in% good_coords$lat_dd & long_dd %in% -good_coords$long_dd)


# Separate data by species
trawlf <- camp_fall_good %>% rename(pa = trawl_pa) %>% select(year, pa, bot.temp) %>% mutate(geartype = "campelen")
AC_stof <- stom_fall_good %>% filter(alt.name == "Atlantic cod") %>% select(year, pa, empty, length, season, content, bot.temp) %>% mutate(geartype = "campelen")
GH_stof <- stom_fall_good %>% filter(alt.name == "Greenland halibut") %>% select(year, pa, empty, length, season, content, bot.temp) %>% mutate(geartype = "campelen")
AP_stof <- stom_fall_good %>% filter(alt.name == "American plaice") %>% select(year, pa, empty, length, season, content, bot.temp) %>% mutate(geartype = "campelen")


# Set up different length bins
ac1 <- AC_stof %>% filter(length > 17 & length <= 45)
ac2 <- AC_stof %>% filter(length > 45) %>% filter(year != 1995) # removed years with fewer than 10 overall observations
ap1 <- AP_stof %>% filter(length > 29) %>% filter(year != 2014)
gh1 <- GH_stof %>% filter(length > 19 & length <= 40)
gh2 <- GH_stof %>% filter(length > 40)


# Set up model data
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


temp <- c(trawlf$bot.temp, 
          ac1$bot.temp, ac2$bot.temp,
          gh1$bot.temp, gh2$bot.temp,
          ap1$bot.temp)


spp_n <- c(length(trawlf$year), 
           sum(length(ac1$year),length(ac2$year)), 
           sum(length(gh1$year),length(gh2$year)),
           sum(length(ap1$year)))
spp_id <- c(0,1,2,3)

