

#import full stomach data with mass and count for each prey item in stomachs, multiple pred obs for preds with >1 item in stomach
full_stom <- read.csv("./Data/NAFC_diet_capelin_COD_TURBOT_PLAICE_2J3KL_2025UPDATE_allprey.csv") %>%
  filter(YEAR >= 1984 & YEAR <= 2020)

prey1 <- full_stom %>% group_by(ID_PRED) %>% slice(which.max(MASS))
#retain only elements in stomach not equal to elements in prey1, rows from df1 which are not in df2
prey2 <- setdiff(full_stom, prey1) %>% group_by(ID_PRED) %>% 
  slice(which.max(MASS)) #find new highest prey mass

uni_full_stom <- full_stom %>% select(!c(TOT_PREY_NUMBER:MASS)) %>% unique() %>% # create unique pred obs dataframe
  mutate(prey1 = NA, # add prey1 and prey2 columns to fill
         prey2 = NA)

# match prey1 and prey2 for each ID_PRED, this step takes a while
pb <- txtProgressBar(min = 0,max = nrow(uni_full_stom),style = 3,width = 50,char = "=")
for(i in 1:nrow(uni_full_stom)){  
  uni_full_stom[i, 'prey1'] <- prey1[prey1$ID_PRED==uni_full_stom$ID_PRED[i], ]$PREY_GR_RAP
  # not all preds have top prey2. match which ones do, NA for those that dont
  if(uni_full_stom$ID_PRED[i] %in% unique(prey2$ID_PRED)){
    uni_full_stom[i, 'prey2'] <- prey2[prey2$ID_PRED==uni_full_stom$ID_PRED[i], ]$PREY_GR_RAP
  } else(uni_full_stom[i, 'prey2'] <- NA)
  setTxtProgressBar(pb, i)
}
close(pb) # Close the connection

save(uni_full_stom, file = "./Data/uni_full_stom.RData")
# load(file = "./Data/uni_full_stom.RData")
