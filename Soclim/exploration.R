library(tidyverse)
library(janitor)
library(readxl)

somerged <- read_csv("Data/Final/somerged")

soclim_btl_v3 <- read_delim("Data/Raw/soclim/soclim_btl_V3.txt", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)
soclim_btl_v9 <- read_delim("Data/Raw/soclim/soclim_btl_V9.txt", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)


ggplot(somerged)+
  geom_point(aes(x = fluo_eco_rfu, y = tchla))

spectre <- read_excel("Biosope/Data/Spectres_annick.xlsx")
spectre <- clean_names(spectre)
spectre440<- filter(spectre, lambda == 440)
spectre470 <- filter(spectre, lambda == 470)


somerged <- somerged %>% mutate(photo_440 = peri * spectre440$peri + but * spectre440$x19_bf + fuco * spectre440$fuco + allo * spectre440$allox + chla * spectre440$chl_a,
                              protect_440 = zea * spectre440$zea,
                              photo_470 = peri * spectre470$peri + but * spectre470$x19_bf + fuco * spectre470$fuco + allo * spectre470$allox + chla * spectre470$chl_a,
                              protect_470 = zea * spectre470$zea, 
                              ratio_photo = photo_440/fluo_eco_rfu,
                              ratio_protect = protect_440/fluo_eco_rfu)

ggplot(somerged)+
  geom_point(aes(x = tchla, y = photo_440))+
  geom_point(aes(x = tchla, y = photo_470), colour = "red")

ggplot(somerged)+
  geom_point(aes(x = fluo_eco_rfu, y = photo_440))+
  geom_point(aes(x = fluo_eco_rfu, y = photo_470), colour = "red")

summary(lm(photo_440~fluo_eco_rfu, data = somerged))
summary(lm(photo_470~fluo_eco_rfu, data = somerged))


