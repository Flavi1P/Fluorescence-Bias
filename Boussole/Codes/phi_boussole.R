library(tidyverse)
library(nnls)
library(Metrics)
theme_set(theme_bw())
source("functions/phi_lm.R")
source("functions/phi_boot.R")


boussole <- read_csv("Boussole/Data/boussole.csv")
boussole$optical_layer <- round(boussole$optical_layer)
boussole <- filter(boussole, optical_layer < 4) %>% mutate(fluo_precise = fluo * 3 - 0.42)

mod <- lm(fluo~ tchla, micro, nano, pico, data = boussole)
cooksd <- cooks.distance(mod)

plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")
abline(h = 4*mean(cooksd, na.rm=T), col="red") 
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")

influential <- as.numeric(names(cooksd)[(cooksd > 4*mean(cooksd, na.rm=T))]) 
head(boussole[influential, ])

boussole <- boussole[-influential,]

phi_bouss_lm <- phi_lm(boussole, "fluo")
phi_bouss_nnls <- phi_boot(boussole, "fluo")

phi_bouss_nnls$se <- ifelse(phi_bouss_nnls$phi > phi_bouss_nnls$se, phi_bouss_nnls$se, phi_bouss_nnls$phi)

ggplot(phi_bouss_nnls, aes(x=size, y = phi, fill = as.factor(optical_layer))) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin = phi-se, ymax = phi+se),position=position_dodge())+
  scale_fill_viridis_d( name = "optical layer")+
  ylab("Phi")+ xlab("size classe")

phi_bouss_lm <- phi_bouss_lm %>% select(phi, optical_layer, size) %>% spread(key = size, value = phi)
names(phi_bouss_lm) <- c("optical_layer", "phi_micro", "phi_nano", "phi_pico") 


phi_bouss_nnls <- phi_bouss_nnls %>% select(phi, optical_layer, size) %>% spread(key = size, value = phi)
names(phi_bouss_nnls) <- c("optical_layer", "phi_micro", "phi_nano", "phi_pico") 


boussole_calibration <- boussole %>% select(micro, nano, pico, amicro, anano, apico, tchla, fluo, optical_layer, fluo_precise) %>% 
  left_join(phi_bouss_nnls) %>% 
  mutate(fluo_calibrate = (fluo/(micro*amicro * phi_micro + nano*anano * phi_nano + pico*apico * phi_pico)))

a <- rmse(boussole_calibration$tchla, boussole_calibration$fluo_calibrate)
b <- rmse(boussole_calibration$tchla, boussole_calibration$fluo_precise)
a/b

phi_bouss_nnls <- phi_bouss_nnls %>% select(phi, optical_layer, size) %>% spread(key = size, value = phi)
names(phi_bouss_nnls) <- c("optical_layer", "phi_micro", "phi_nano", "phi_pico") 

boussole_calibration_nnls <- boussole %>% select(micro, nano, pico, amicro, anano, apico, tchla, fluo_precise, optical_layer, fluo_precise) %>% 
  left_join(phi_bouss_nnls) %>% 
  mutate(fluo_calibrate = (fluo_precise/(micro*amicro * phi_micro + nano*anano * phi_nano + pico*apico * phi_pico)))

c <- rmse(boussole_calibration_nnls$tchla, boussole_calibration_nnls$fluo_calibrate)
d <- rmse(boussole_calibration_nnls$tchla, boussole_calibration_nnls$fluo_precise)
c/d

ggplot(filter(boussole_calibration_nnls, fluo_precise>0))+
  geom_point(aes(x = tchla, y = fluo_precise), colour = "Blue")+
  geom_point(aes(x = tchla, y = fluo_calibrate), colour = "Red")+
  geom_line(aes(x = tchla, y = tchla))+
  coord_trans(x = "log", y = "log")

summary(lm(tchla~fluo, data = boussole_calibration))
