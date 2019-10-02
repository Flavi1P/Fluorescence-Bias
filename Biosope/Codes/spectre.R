library(tidyverse)
library(readxl)
library(janitor)

biosope <- read_csv("Biosope/Data/biosope")
spectre <- read_excel("Biosope/Data/Spectres_annick.xlsx")
spectre <- clean_names(spectre)
spectre <- filter(spectre, lambda == 440)

biosope <- biosope %>% mutate(photosynthetic = peri * spectre$peri + but * spectre$x19_bf + fuco * spectre$fuco + allo * spectre$allox + chla * spectre$chl_a,
                              protect = diad * spectre$diad + zea * spectre$zea)

ggplot(filter(biosope, depth < 50), aes(x = fluo_urel))+
  geom_point(aes(y = photosynthetic), colour = "green")+
  geom_point(aes(y = protect), colour = "black")+
  xlim(0,1) + ylim(0,0.025)


summary(lm(biosope$photosynthetic~biosope$fluo_urel))


argo <- read_csv("argo/Data/merged_argo")

argo <- argo %>% mutate(photosynthetic = peri * spectre$peri + but * spectre$x19_bf + fuco * spectre$fuco + allo * spectre$allox + tchla * spectre$chl_a,
                              protect = zea * spectre$zea)

ggplot(filter(argo, depth < 50), aes(x = chla))+
  geom_point(aes(y = photosynthetic), colour = "green")+
  geom_point(aes(y = protect), colour = "black")+
  xlim(0,1) + ylim(0,0.025)
