library(tidyverse)
library(janitor)

#read data
lov_hplc <- read_delim("Data/lov_hplc/PigmentsDB-LOVonly_190617.txt", 
                                                    "\t", escape_double = FALSE, trim_ws = TRUE)
lov_hplc <- clean_names(lov_hplc)
map <- read_csv("Data/map_vec")

colnames(lov_hplc)

#We select columns that we need

lov_hplc_short <- lov_hplc %>% select(project, cruise, "day" = date, "month" = x9, "year" = x10, sdy, lat, lon, depth, peri, but, fuco, hex, allo, dv_chla, chla)

#compute the sum of my pigments

lov_hplc_short <- lov_hplc_short %>% mutate(pigsum = chla + dv_chla + fuco + peri + but + hex + allo)

table(is.na(lov_hplc_short$pigsum))

#remove the NA

lov_hplc_clean <- lov_hplc_short[- which(is.na(lov_hplc_short$pigsum)),]

#plot the data

ggplot()+
  geom_point(data = lov_hplc_clean, aes(x = lon, y = lat), shape = 21, fill = "#1c9099", colour = "black")+
  geom_polygon(data = map, aes(x = long, y = lat, group = group), fill = "Grey")+
  coord_equal()+
  theme_minimal()+
  ggtitle("Map of HPLC sample from lov database")
  
