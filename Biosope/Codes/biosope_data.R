library(tidyverse)
library(data.table)
library(readxl)
library(lubridate)

biosoperaw <- read_excel("Data/Raw/Full/BDD_HPLC_fluo_cruises/Biosope_odv_btl.xlsx")
names(biosoperaw) <- tolower(colnames(biosoperaw))

biosope <- biosoperaw %>% select(campagne, date, longitude, latitude, btle, "p[dbar]", "t[its90]", salinite,
                                 "ox[ml/l]", "fluo[v]", "f[urel]", "isus[v]", nitrates)
names(biosope) <- c("cruise", "date", "lon", "lat", "btl_nbr", "pressure", "temp", "salinity", "ox_ml_l",
                    "fluo_v", "fluo_urel", "isus_v", "nitrates")
biosope$lat <- substr(biosope$lat, 1, 7)
biosope$lat <- as.numeric(sub(",", ".", biosope$lat, fixed = TRUE))
biosope$lon <- substr(biosope$lon, 1, 8)

biosope$lon <- as.numeric(sub(",", ".", biosope$lon, fixed = TRUE))
biosope$date <- as_date(biosope$date, tz = "UTC", format = "%m/%d/%Y")

biosope$lat <- round(biosope$lat, 2)
biosope$lon <- round(biosope$lon, 2)

map_vec <- read_csv("Data/map_vec")

ggplot() + 
  geom_polygon( data = filter(map_vec, long < -30 & long > -160 & lat > -60 & lat < 10) , aes(x = long,y = lat, group = group), fill = "Grey")+
  geom_point(aes(lon, lat), data = biosope, size = 2)+
  coord_quickmap()+
  theme_bw(base_size = 20)+
  ylab("Latitude")+ xlab("Longitude")+
  scale_color_viridis_c()

write_csv(biosope, "Data/Process/Biosope/biosope")