library(tidyverse)
library(readxl)
library(janitor)
library(lubridate)
library(gridExtra)

#open biosope data
biosope <- read_csv("Biosope/Data/biosope")

#open pigments absorbtion
spectre <- read_excel("Biosope/Data/Spectres_annick.xlsx")
spectre <- clean_names(spectre)

#select the absorbtion at 440 and 470nm
spectre440<- filter(spectre, lambda == 440)
spectre470 <- filter(spectre, lambda == 470)

#create columns that correspond to the total photosynthetic absorbance and non photosynthetic absorbance at 440 and 470. Create also a ratio between the two photosynthetic absorbtion
biosope <- biosope %>% mutate(photo_440 = peri * spectre440$peri + but * spectre440$x19_bf + hex * spectre440$x19_hf + fuco * spectre440$fuco + allo * spectre440$allox + chla * spectre440$chl_a + dvchla * spectre440$dv_chla,
                              protect_440 = zea * spectre440$zea,
                              photo_470 = peri * spectre470$peri + but * spectre470$x19_bf + hex * spectre470$x19_hf + fuco * spectre470$fuco + allo * spectre470$allox + chla * spectre470$chl_a + dvchla * spectre470$dv_chla,
                              protect_470 = zea * spectre470$zea, 
                              ratio = photo_440/photo_470)

#Create a df to ease the visualisation
biosope_visu <- select(biosope, depth, site, lon, tchla, ratio, photo_440, photo_470, ratio)

#see the ratio on a profile
ggplot(filter(biosope_visu, site == "St8"))+
  geom_point(aes(x = ratio, y = -depth))+
  theme_bw()
#see the chla on this profile
ggplot(filter(biosope_visu, site == "St8"))+
  geom_point(aes(x = tchla, y = -depth))+
  theme_bw()

#see the ratio as function of chla on the surface of the campain
ggplot(filter(biosope, depth < 20))+
  geom_point(aes(x = tchla, y = ratio))+
  theme_bw()

#see the ratio as function of fluorescence on the surface of the campain
ggplot(filter(biosope, depth < 20), aes(x = fluo_urel))+
  geom_point(aes(y = ratio))+
  theme_bw()

#create an empty df of biosope_visu with the same columns
biosope_first_point <- biosope_visu[1,]
biosope_first_point <- biosope_first_point[-1,]

#loop to select only the surface point in each station and bind all of that in the precedent df
for(i in unique(round(biosope$lon, 1))){
  t <- filter(biosope_visu, round(lon, 1) == i)
  t <- filter(t, depth== min(t$depth))
  biosope_first_point <- bind_rows(biosope_first_point, t)
}

#visualisation of the ratio over the longitude on the surface and the chla

gratio <- ggplot(biosope_first_point)+
  geom_point(aes(x = lon, y = ratio), colour = "#fdae6b")+
  theme_bw()+
  ylim(1,3)

gchla <- ggplot(biosope_first_point)+
  geom_point(aes(x = lon, y = tchla), colour = "#a1d99b")+
  theme_bw()+
  ylim(0,0.75)

grid.arrange(gchla, gratio)


  # argo <- read_csv("argo/Data/merged_argo")
# 
# argo <- argo %>% mutate(photosynthetic = peri * spectre$peri + but * spectre$x19_bf + fuco * spectre$fuco + allo * spectre$allox + tchla * spectre$chl_a,
#                               protect = zea * spectre$zea)
# 
# ggplot(filter(argo), aes(x = chla))+
#   geom_point(aes(y = photosynthetic), colour = "green")+
#   geom_point(aes(y = protect), colour = "black")+
#   xlim(0,1) + ylim(0,0.025)
# 
