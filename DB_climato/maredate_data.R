library(tidyverse)
library(readxl)
library(janitor)

#open maredat data and clean names
maredat <- read_excel("Data/maredat/MAREDAT_pigments_master_file_Updated210313.xls")
maredat <- clean_names(maredat)
map <- read_csv("Data/map_vec")

colnames(maredat)

#we select columns that we need, convert the concentration in mg.m-3
#we use the merged depth which is a merge between pressure and depth in meter, where pressure given priority if both given

maredat_short <- maredat %>% select(database_num, sample_num, exp, ctd, bottle_num, lat, long, month, day, year, merged_depth_m,
                                    number_of_pigments, total_chla_mg_m3, chla_ng_l, dv_chla_ng_l, x19hex_ng_l, x19but_ng_l, fucox_ng_l, perid_ng_l, allox_ng_l) %>% 
  mutate(depth = merged_depth_m,
         tchla = total_chla_mg_m3,
         chla = chla_ng_l * 1000,
         dvchla = dv_chla_ng_l * 1000,
         hex = x19hex_ng_l * 1000,
         but = x19but_ng_l * 1000,
         fuco = fucox_ng_l * 1000,
         peri = perid_ng_l * 1000,
         allo = allox_ng_l * 1000) %>% 
  select(- merged_depth_m, - chla_ng_l, - dv_chla_ng_l, - x19hex_ng_l, - x19but_ng_l, - fucox_ng_l, - perid_ng_l, - allox_ng_l)

ggplot()+
  geom_point(data = maredat_short, aes(x = long, y = lat, colour = exp))+
  geom_polygon(data = map, aes(x = long, y = lat, group = group))+
  coord_equal()

#compute the sum of my pigment, if Na return, then unusable sample

maredat_short <- maredat_short %>% mutate(pigsum = chla + dvchla + hex + but + fuco + peri + allo)
table(is.na(maredat_short$pigsum))

#I remove the NA
maredat_clean <- maredat_short[- which(is.na(maredat_short$pigsum)),]

#plot it
ggplot()+
  geom_point(data = maredat_clean, aes(x = long, y = lat), shape = 21, fill = "#1c9099", colour = "black")+
  geom_polygon(data = map, aes(x = long, y = lat, group = group), fill = "Grey")+
  coord_equal()+
  theme_minimal()+
  ggtitle("Map of Maredat samples we can use")

#We have 16 276 samples ! 

