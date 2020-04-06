library(tidyverse)
library(readxl)
library(janitor)
library(lubridate)

biosope <- read_csv("Biosope/Data/biosope")
spectre <- read_excel("Biosope/Data/Spectres_annick.xlsx")
spectre <- clean_names(spectre)
spectre440<- filter(spectre, lambda == 440)
spectre470 <- filter(spectre, lambda == 470)


biosope <- biosope %>% mutate(photo_440 = peri * spectre440$peri + but * spectre440$x19_bf + fuco * spectre440$fuco +  chla * spectre440$chl_a + hex * spectre440$x19_hf + dvchlb * spectre440$dv_chlb + dvchla * spectre440$dv_chla,
                              protect_440 = zea * spectre440$zea,
                              photo_470 = peri * spectre470$peri + but * spectre470$x19_bf + fuco * spectre470$fuco + chla * spectre470$chl_a + hex * spectre440$x19_hf + dvchlb * spectre440$dv_chlb + dvchla * spectre440$dv_chla,
                              protect_470 = zea * spectre470$zea, 
                              ratio_photo = photo_440/fluo_urel,
                              ratio_protect = protect_440/fluo_urel)

ggplot(filter(biosope, depth < 20), aes(y = fluo_urel))+
  geom_point(aes(x = photo_440), colour = "green")+
  geom_point(aes(x = protect_440), colour = "black")+
  ylim(0,1) + xlim(0,0.025)

ggplot(filter(biosope, depth < 20), aes(y = fluo_urel))+
  geom_point(aes(x = photo_440), colour = "green")+
  geom_point(aes(x = photo_470), colour = "black")+
  ylim(0,1) + xlim(0,0.025)

ggplot(biosope)+
  geom_point(aes(x = lon, y = -depth))+
  geom_point(aes(x = lon, y = -ze), colour = "red")

ggplot(filter(biosope, site == "St8"))+
  geom_point(aes(x = photo_440/photo_470, y = -depth))+
  theme_bw()

ggplot(filter(biosope, site == "St8"))+
  geom_point(aes(x = tchla, y = -depth))+
  theme_bw()

ggplot(filter(biosope, site %in% c("GYR2" , "GYR3" , "GYR4")))+
  geom_point(aes(x = photo_440/photo_470, y = -depth))+
  theme_bw()

ggplot(filter(biosope, site %in% c("GYR2" , "GYR3" , "GYR4")))+
  geom_point(aes(x = tchla, y = -depth))+
  theme_bw()

ggplot(filter(biosope, depth < 20))+
  geom_point(aes(x = tchla, y = photo_440/photo_470))+
  theme_bw()

ggplot(filter(biosope, depth < 20))+
  geom_point(aes(x = lon, y = photo_440/photo_470))+
  theme_bw()

ggplot(filter(biosope, depth < 20), aes(x = fluo_urel))+
  geom_point(aes(y = photo_440/photo_470))+
  theme_bw()

ggplot(filter(biosope, depth < 20), aes(x = tchla))+
  geom_point(aes(y = photo_440/photo_470))+
  theme_bw()


argo <- read_csv("argo/Data/merged_argo")

argo <- argo %>% mutate(photosynthetic = peri * spectre$peri + but * spectre$x19_bf + fuco * spectre$fuco + allo * spectre$allox + tchla * spectre$chl_a,
                              protect = zea * spectre$zea)

ggplot(filter(argo), aes(x = chla))+
  geom_point(aes(y = photosynthetic), colour = "green")+
  geom_point(aes(y = protect), colour = "black")+
  xlim(0,1) + ylim(0,0.025)

boussole <- read_csv("Boussole/Data/boussole.csv")

boussole <- boussole %>% mutate(photo_440 = peri * spectre440$peri + but * spectre440$x19_bf + fuco * spectre440$fuco + allo * spectre440$allox + tchla * spectre440$chl_a + ,
                              protect_440 = zea * spectre440$zea,
                              photo_470 = peri * spectre470$peri + but * spectre470$x19_bf + fuco * spectre470$fuco + allo * spectre470$allox + tchla * spectre470$chl_a,
                              protect_470 = zea * spectre470$zea, 
                              ratio_photo = photo_440/fluo,
                              ratio_protect = protect_440/fluo)

ggplot(filter(boussole, depth < 20), aes(y = fluo))+
  geom_point(aes(x = photo_440), colour = "green")

ggplot(filter(boussole, depth < 20))+
  geom_point(aes(x = date, y = ratio_photo))+
  geom_point(aes(x = date, y = zea), colour = "red")

ggplot(filter(boussole, depth < 20))+
  geom_point(aes(x = date, y = ratio_photo))


boussole$year <- year(boussole$date)

ggplot(filter(boussole, depth < 15))+
  geom_point(aes(x = date, y = photo_440/photo_470))+
  facet_wrap(facets = "year", ncol = 1, scales = "free_x")+
  ylab("Rapport Absorbance/Fluo")

ggplot(filter(boussole, depth < 20))+
  geom_point(aes(x = date, y = ratio_protect))
