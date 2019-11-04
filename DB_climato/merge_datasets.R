library(tidyverse)
library(lubridate)

#read the df that have been shaped in the lov_hplc.R and maredate_data.R scripts
lov <- read_csv("DB_climato/Data/lov_hplc.csv")
maredat <- read_csv("DB_climato/Data/maredat_hplc.csv")
map <- read_csv("Data/map_vec")


#create a date column for both dataset
lov$date <- ymd(paste(lov$year, lov$month, lov$day, "-"))
maredat$date <- ymd(paste(maredat$year, maredat$month, maredat$day, "-"))

#select the same columns in each dataset
lov <- lov %>% select(date, lon, lat, depth, peri, but, fuco, hex, allo, dv_chla, chla)
maredat <- maredat %>% select(date, "lon" = long, lat, depth, peri, but, fuco, hex, allo, dv_chla = dvchla, chla)

#create one unique dataframe
combine_hplc <- bind_rows(lov, maredat)

#keep the info of the database
combine_hplc$dataset <- c(rep("lov", nrow(lov)), rep("maredat", nrow(maredat)))

ggplot()+
  geom_point(aes(x = lon, y = lat, colour = dataset), alpha = 1/10, data = combine_hplc)+
  geom_polygon(aes(x = long, y = lat, group = group), data = map)+
  coord_equal()

#Now we are looking for samples that are present in each dataset
#create an id that correspond either to the sample ("test") or to the profile ("profile_id")
combine_hplc <- combine_hplc %>% mutate(lon = round(lon, digit = 3),
                                        lat = round(lat, 3),
                                        depth = round(depth, 1),
                                        test = paste(date, lon, lat, depth, sep = "/"),
                                        profile_id = paste(date, lon, lat, sep = "/"))

table(combine_hplc$profile_id)


duplicate <- combine_hplc[which(duplicated(combine_hplc$test)),]

ggplot()+
  geom_point(aes(x = lon, y = lat, colour = dataset), data = combine_hplc)+
  geom_point(aes(x = lon, y = lat), data = duplicate, size = 2)+
  geom_polygon(aes(x = long, y = lat, group = group), data = map)+
  coord_equal()

combine_hplc_clean <- combine_hplc[- which(duplicated(combine_hplc$test)),]

#we compute the week of the year to match with climato data
combine_hplc_clean$week <- week(combine_hplc_clean$date)

#write_csv(combine_hplc_clean, "DB_climato/Data/global_hplc")




