library(tidyverse)
library(janitor)

#read data
lov_hplc <- read_delim("Data/lov_hplc/PigmentsDB-LOVonly_190617.txt", 
                                                    "\t", escape_double = FALSE, trim_ws = TRUE)
lov_hplc <- clean_names(lov_hplc)
map <- read_csv("Data/map_vec")

colnames(lov_hplc)

#We select columns that we need

lov_hplc_short <- lov_hplc %>%
  select(project, cruise, "day" = date, "month" = x9, "year" = x10, sdy, lat, lon, depth, peri, but, fuco, hex, allo, dv_chla, chla)

lov_hplc_short$profile_id <- paste(lov_hplc_short$day, lov_hplc_clean$year, lov_hplc_short$month, round(lov_hplc_short$lon, 2), sep  = "/") 

length(unique(lov_hplc_short$profile_id))
    #compute the sum of my pigments

lov_hplc_short <- lov_hplc_short %>% mutate(pigsum = chla + dv_chla + fuco + peri + but + hex + allo)

table(is.na(lov_hplc_short$pigsum))

#remove the NA

lov_hplc_clean <- lov_hplc_short[- which(is.na(lov_hplc_short$pigsum)),]

lov_hplc_clean$profile_id <- paste(lov_hplc_clean$day, lov_hplc_clean$year, lov_hplc_clean$month, round(lov_hplc_clean$lon, 2), sep  = "/") 

length(unique(lov_hplc_clean$profile_id))

LON_LAT_MONTH <- unique(lov_hplc_short[,c("lon","lat","month")]) %>% mutate(nprof = seq(1, nrow(.), by = 1)) #create a num of profile

hplc2 <- left_join(lov_hplc_short, LON_LAT_MONTH) #add the num of the profile on hplc2
short_prof <- filter(hplc2, nprof %in% which(table(hplc2$nprof)<4)) #create a df with profiles that have less than 4 data
length(unique(short_prof$nprof))

#plot the data

ggplot()+
  geom_point(data = lov_hplc_clean, aes(x = lon, y = lat), shape = 21, fill = "#1c9099", colour = "black")+
  geom_polygon(data = map, aes(x = long, y = lat, group = group), fill = "Grey")+
  coord_equal()+
  theme_minimal()+
  ggtitle("Map of HPLC sample from lov database")

#write_csv(lov_hplc_clean, "DB_climato/Data/lov_hplc.csv")  
