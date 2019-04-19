library(tidyverse)
library(zoo)
library(vegan)

source("functions/phi_simple.R")

pigments <- c("fuco", "peri", "hex", "but", "allo", "tchlb", "zea")

biosope <- read_csv("Biosope/Data/biosope")
biosope$optical_layer <- round(biosope$optical_layer)

biosope <- filter(biosope, optical_layer < 4) %>% mutate(microquanti = micro * tchla,
                                                         picoquanti = pico * tchla,
                                                         nanoquanti = nano * tchla)
phi_biosope <- phi_simple(biosope, "fluo_urel")

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
  mutate(fluo_calibrate = (fluo_urel/(micro* phi_micro + nano* phi_nano + pico* phi_pico)))

a <- rmse(biosope_calibration$tchla, biosope_calibration$fluo_calibrate)
b <- rmse(biosope_calibration$tchla, biosope_calibration$fluo_urel)
a/b

ggplot(biosope_calibration)+
  geom_point(aes(x = tchla, y = fluo_calibrate), colour = "Red")+
  geom_point(aes(x = tchla, y = fluo_urel, colour = micro))+
  scale_color_viridis_c()


#test is on sampled dataset####

biosope <- read_csv("Biosope/Data/biosope")

biosope$optical_layer <- round(biosope$optical_layer)

biosope <- filter(biosope, optical_layer < 4) %>% mutate(microquanti = micro * tchla,
                                                         nanoquanti = nano * tchla,
                                                         picoquanti = pico * tchla)

biosope_train <- sample_frac(biosope, 0.8)
biosope_test <- sample_frac(biosope, 0.2)

phi_biosope <- phi_simple(biosope_train, "fluo_urel")

ggplot(phi_biosope, aes(x=size, y = phi, fill = as.factor(optical_layer))) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin = phi-se, ymax = phi+se),position=position_dodge())+
  scale_fill_viridis_d( name = "optical layer")+
  ylab("Phi")+ xlab("size classe")+
  ggtitle("Fluorescent yield, Biosope")

phi_biosope <- phi_biosope %>% select(phi, optical_layer, size) %>% spread(key = size, value = phi)
names(phi_biosope) <- c("optical_layer", "phi_micro", "phi_nano", "phi_pico") 

biosope_calibration <- biosope_test %>% select(micro, nano, pico, amicro, anano, apico, tchla, fluo_urel, optical_layer) %>% 
  left_join(phi_biosope) %>% 
  mutate(fluo_calibrate = (fluo_urel/(micro * phi_micro + nano * phi_nano + pico * phi_pico)))

a <- rmse(biosope_calibration$tchla, biosope_calibration$fluo_calibrate)
b <- rmse(biosope_calibration$tchla, biosope_calibration$fluo_urel)
a/b

rmse_computed <- data.frame("run" = c(1:10000), "rmse" = 0)

for (i in rmse_computed$run){
  
  biosope_train <- sample_frac(biosope, 0.8)
  biosope_test <- sample_frac(biosope, 0.2)
  
  phi_biosope <- phi_simple(biosope_train, "fluo_urel")
  coeff_lm <- summary(lm(fluo_urel~tchla, data = biosope_train))$coefficient[2,1]
  
  biosope_test$fluo <- biosope_test$fluo_urel*coeff_lm
  
  phi_biosope <- phi_biosope %>% select(phi, optical_layer, size) %>% spread(key = size, value = phi)
  names(phi_biosope) <- c("optical_layer", "phi_micro", "phi_nano", "phi_pico") 
  
  biosope_calibration <- biosope_test %>% select(micro, nano, pico, amicro, anano, apico, tchla, fluo_urel, optical_layer, fluo) %>% 
    left_join(phi_biosope, by = "optical_layer") %>% 
    mutate(fluo_calibrate = (fluo_urel/(micro* phi_micro + nano * phi_nano + pico * phi_pico)))
  
  a <- rmse(biosope_calibration$tchla, biosope_calibration$fluo_calibrate)
  b <- rmse(biosope_calibration$tchla, biosope_calibration$fluo)
  rmse_computed$rmse[i] <- a/b
  
}
mean(rmse_computed$rmse)
sd(rmse_computed$rmse)
