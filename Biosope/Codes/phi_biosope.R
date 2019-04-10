library(tidyverse)
library(zoo)
library(vegan)

source("functions/phi_lm.R")
source("functions/phi_boot.R")

source("functions/profile_numb.R")
source("functions/zeu_moma.R")

pigments <- c("fuco", "peri", "hex", "but", "allo", "tchlb", "zea")

biosope <- read_csv("Biosope/Data/biosope")
biosope$optical_layer <- round(biosope$optical_layer)

biosope <- filter(biosope, optical_layer < 4)
phi_biosope <- phi_lm(biosope, "fluo_urel")

ggplot(phi_biosope, aes(x=size, y = phi, fill = as.factor(optical_layer))) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin = phi-se, ymax = phi+se),position=position_dodge())+
  scale_fill_viridis_d( name = "optical layer")+
  ylab("Phi")+ xlab("size classe")+
  ggtitle("Fluorescent yield, Biosope")

phi_biosope <- phi_biosope %>% select(phi, optical_layer, size) %>% spread(key = size, value = phi)
names(phi_biosope) <- c("optical_layer", "phi_micro", "phi_nano", "phi_pico") 

biosope_calibration <- biosope %>% select(micro, nano, pico, amicro, anano, apico, tchla, fluo_urel, optical_layer) %>% 
  left_join(phi_biosope) %>% 
  mutate(fluo_calibrate = (fluo_urel/(micro*amicro * phi_micro + nano*anano * phi_nano + pico*apico * phi_pico)))

a <- rmse(biosope_calibration$tchla, biosope_calibration$fluo_calibrate)
b <- rmse(biosope_calibration$tchla, biosope_calibration$fluo_urel)
a/b

ggplot(biosope_calibration)+
  geom_point(aes(x = tchla, y = fluo_calibrate), colour = "Green")+
  geom_point(aes(x = tchla, y = fluo_urel))

biosope_calibration <- biosope_calibration %>% mutate(ratio_model = micro * phi_micro + nano * phi_nano + pico * phi_pico)

ggplot(biosope_calibration)+
  geom_density(aes(x =ratio_model))


