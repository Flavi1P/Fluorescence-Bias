library(tidyverse)
library(readxl)
library(janitor)
library(lubridate)
library(gridExtra)
library(Metrics)

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

#analyse of the variance bewteen the two aPS

#create a df to use geom_boxplot

absorbtion <- data.frame("lambda" = as.factor(c(rep(440, 433), rep(470, 433))), "aps" = c(biosope$photo_440, biosope$photo_470))

ggplot(absorbtion)+
  geom_boxplot(aes(x = lambda, y = aps))+
  theme_classic()+
  ylab("aPS")


#Estimation of the error of linear reconstruction of absorbtion

#read the absorbtion data

biosope_ap <- read_excel("Biosope/Data/biosope_ap_spectro_98_2005_09_12.xls")
names(biosope_ap) <- tolower(names(biosope_ap))

biosope_nap <- read_excel("Biosope/Data/biosope_aNAP_spectro_98_2005_09_12.xls")
names(biosope_nap) <- tolower(names(biosope_nap))

#select only 440 and 470nm wavelength
biosope_ap <- select(biosope_ap, station, ctd, "bottle no.", depth, "440", "470")
biosope_nap <- select(biosope_nap, station, ctd, "bottle no.", depth, "440", "470")

#rename some columns to match with biosope df
biosope_ap <- rename(biosope_ap, btl_nbr = "bottle no.", site = station, aap470 = "470", aap440 = "440")
biosope_nap <- rename(biosope_nap, btl_nbr = "bottle no.", site = station, anp470 = "470", anp440 = "440")

#make biosope nap 440 et 470 numeric
biosope_nap$anp470 <- as.numeric(biosope_nap$anp470)
biosope_nap$anp440 <- as.numeric(biosope_nap$anp440)


#make btl nbr numeric to enable left join
biosope_ap$btl_nbr <- as.numeric(biosope_ap$btl_nbr)
biosope_nap$btl_nbr <- as.numeric(biosope_nap$btl_nbr)

#make station names the same
biosope_ap$site <- gsub("^STB", "St", biosope_ap$site)
biosope_nap$site <- gsub("^STB", "St", biosope_nap$site)

#merge all dataset
biosope_absorbtion <- left_join(biosope, biosope_ap, by = c("site", "btl_nbr"))
biosope_absorbtion <- left_join(biosope_absorbtion, biosope_nap, by = c("site", "btl_nbr"))

biosope_absorbtion <- biosope_absorbtion %>% mutate(atot_440 = photo_440 + protect_440,
                                                    atot_470 = photo_470 + protect_470,
                                                    aph_440 = aap440 - anp440,
                                                    aph_470 = aap470 - anp470)


#plot the estimation of absorbtion vs the real absorbtion
ggplot(biosope_absorbtion)+
  geom_point(aes(x = aph_440, y = atot_440))

ggplot(biosope_absorbtion)+
  geom_point(aes(x = aph_470, y = atot_470))

#compute the mean absolute percentage error
biosope_absorbtion <- biosope_absorbtion[-which(is.na(biosope_absorbtion$aph_470)),]
mape(biosope_absorbtion$aph_470, biosope_absorbtion$atot_470)

mape(biosope_absorbtion$aph_440, biosope_absorbtion$atot_440)

#une forte erreur due à une pente de 1.7 sur la régression, mais un r² de 0.91

