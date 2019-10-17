library(tidyverse)
library(zoo)
library(vegan)

source("functions/phi_lm.R")
source("functions/phi_boot.R")

pigments <- c("fuco", "peri", "hex", "but", "allo", "tchlb", "zea")

biosope <- read_csv("Biosope/Data/biosope")
biosope$optical_layer <- round(biosope$optical_layer)

biosope <- filter(biosope, optical_layer < 4)
phi_biosope <- phi_lm(biosope, "fluo_urel")
phi_biosope_tchla <- phi_lm(biosope, "tchla")

yield_ratio <- phi_biosope$phi/phi_biosope_tchla$phi
phi_biosope$yield_ratio <- yield_ratio

ggplot(phi_biosope, aes(x=size, y = phi, fill = as.factor(optical_layer))) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin = phi-se, ymax = phi+se),position=position_dodge())+
  scale_fill_viridis_d( name = "optical layer")+
  ylab("Phi")+ xlab("size classe")+
  ggtitle("Fluorescent yield, Biosope")

ggplot(phi_biosope, aes(x=size, y = yield_ratio, fill = as.factor(optical_layer))) +
  geom_bar(position=position_dodge(), stat="identity") +
  scale_fill_viridis_d( name = "optical layer")+
  ylab("Ratio")+ xlab("size classe")+
  ggtitle("Fluorescent yield ratio, Biosope")

phi_biosope <- phi_biosope %>% select(phi, optical_layer, size) %>% spread(key = size, value = phi)
names(phi_biosope) <- c("optical_layer", "phi_micro", "phi_nano", "phi_pico") 

biosope_calibration <- biosope %>% select(micro, nano, pico, amicro, anano, apico, tchla, fluo_urel, optical_layer) %>% 
  left_join(phi_biosope) %>% 
  mutate(fluo_calibrate = (fluo_urel/(micro*amicro * phi_micro + nano*anano * phi_nano + pico*apico * phi_pico)))

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

biosope <- filter(biosope, optical_layer < 4)

biosope_train <- sample_frac(biosope, 0.8)
biosope_test <- sample_frac(biosope, 0.2)

phi_biosope <- phi_boot(biosope_train, "fluo_urel")

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
  mutate(fluo_calibrate = (fluo_urel/(micro*amicro * phi_micro + nano*anano * phi_nano + pico*apico * phi_pico)))

a <- rmse(biosope_calibration$tchla, biosope_calibration$fluo_calibrate)
b <- rmse(biosope_calibration$tchla, biosope_calibration$fluo_urel)
a/b

rmse_computed <- data.frame("run" = c(1:100), "rmse" = 0)

for (i in rmse_computed$run){
  
  biosope_train <- sample_frac(biosope, 0.8)
  biosope_test <- sample_frac(biosope, 0.2)
  
  phi_biosope <- phi_boot(biosope_train, "fluo_urel")
  coeff_lm <- summary(lm(fluo_urel~tchla, data = biosope_train))$coefficient[2,1]
  
  biosope_test$fluo <- biosope_test$fluo_urel*coeff_lm
  
  phi_biosope <- phi_biosope %>% select(phi, optical_layer, size) %>% spread(key = size, value = phi)
  names(phi_biosope) <- c("optical_layer", "phi_micro", "phi_nano", "phi_pico") 
  
  biosope_calibration <- biosope_test %>% select(micro, nano, pico, amicro, anano, apico, tchla, fluo_urel, optical_layer, fluo) %>% 
    left_join(phi_biosope, by = "optical_layer") %>% 
    mutate(fluo_calibrate = (fluo_urel/(micro*amicro * phi_micro + nano*anano * phi_nano + pico*apico * phi_pico)))
  
  a <- rmse(biosope_calibration$tchla, biosope_calibration$fluo_calibrate)
  b <- rmse(biosope_calibration$tchla, biosope_calibration$fluo)
  rmse_computed$rmse[i] <- a/b
  
}
mean(rmse_computed$rmse)
sd(rmse_computed$rmse)
