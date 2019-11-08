library(tidyverse)
library(lubridate)

hplc_db <- read_csv("DB_climato/Data/database_final")
argo <- read_csv("Data/merged_argo")
map <- read_csv("Data/map_vec")

#plot the cooccurence of the two dataset
ggplot()+
  geom_point(aes(x = lon.y, y = lat.y), data = argo, col = "red")+
  geom_point(aes(x = lon, y = lat), data = hplc_db)+
 geom_polygon(aes(x = long, y = lat, group = group), data = map)

#create an id for each profile
argo$profile_id <- paste(round(argo$lon.x, 2), round(argo$lat.x, 2), month(argo$date), year(argo$date), sep = "/")
hplc_db$profile_id <- paste(hplc_db$lon, hplc_db$lat, hplc_db$month, year(hplc_db$date), sep = "/")

#remove the profiles with same id in hplc database
hplc_clean <- filter(hplc_db, ! profile_id %in% argo$profile_id) %>% select(- profile_id)

#write_csv(hplc_db, "DB_climato/Data/database_final")





