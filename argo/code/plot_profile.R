library(tidyverse)
source("functions/profile_numb.r")

pigments <- c("fuco", "peri", "hex", "but", "allo", "tchlb", "zea")

hplc <- read_csv("Data/argo/hplc_argo")
ref <- read_csv("Data/argo/ref.csv")
ref_bis <- read_csv("Data/argo/ref_bis")
map_vec <- read_csv("Data/map_vec")
first_profiles <- read_csv("Data/argo/first_profiles")
first_profiles <- first_profiles[-1,]
first_profiles$date <- date(first_profiles$date)
ref <- bind_rows(ref, ref_bis)

first_profiles$depth <- round(first_profiles$pres)

hplc <- left_join(hplc, ref, by = c("id" = "lovbio"))
first_profiles$number <- as.numeric(first_profiles$id)
first_profiles <- left_join(first_profiles, ref)
first_profiles <- left_join(first_profiles, ref_bis)

merged <- read_csv("Data/merged_argo")

merged <- select(merged, lon.y, lat.y, chla, tchla, depth, lovbio)

profiles <- left_join(first_profiles, merged, by = c("chla", "lon" = "lon.y", "lat" = "lat.y", "depth"))


group_1 <- c("lovbio043b", "lovbio057b")
group_2 <- c("lovbio050b", "lovbio082b")
group_3 <- c("lovbio100c", "lovbio107c")

ggplot()+
  geom_path(aes(x = chla_adjusted, y = -pres, colour = lovbio), data = filter(first_profiles, lovbio %in% group_1))+
  geom_point(aes(x = tchla, y = - depth), data = filter(merged, lovbio %in% group_1))+
  ylim(-250,0)

ggplot()+
  geom_path(aes(x = chla_adjusted, y = -pres, colour = lovbio), data = filter(first_profiles, lovbio %in% group_2))+
  geom_point(aes(x = tchla, y = - depth), data = filter(merged, lovbio %in% group_2))+
  ylim(-250,0)

ggplot()+
  geom_path(aes(x = chla_adjusted, y = -pres, colour = lovbio), data = filter(first_profiles, lovbio %in% group_3))+
  geom_point(aes(x = tchla, y = - depth), data = filter(merged, lovbio %in% group_3))


ggplot()+
  geom_polygon(aes(x = long, y = lat, group = group), data = map_vec)+
  geom_text_repel(aes(x = lon, y = lat, label = lovbio), data = first_profiles, colour = "Red")+
  coord_quickmap()

ggplot()+
  geom_point(aes(x = tchla, y = chla), data = merged)+
  geom_point(aes(x = tchla, y = chla), data = filter(merged, lovbio == i), colour = "Red")
