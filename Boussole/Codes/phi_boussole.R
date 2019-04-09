library(tidyverse)

source("functions/phi_lm.R")

boussole <- read_csv("Boussole/Data/boussole.csv")
boussole$optical_layer <- round(boussole$optical_layer)
boussole <- filter(boussole, optical_layer < 4)
phi_bouss <- phi_lm(boussole, "fluo")

ggplot(phi_bouss, aes(x=size, y = phi, fill = as.factor(optical_layer))) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin = phi-se, ymax = phi+se),position=position_dodge())+
  scale_fill_viridis_d( name = "optical layer")+
  ylab("Phi")+ xlab("size classe")
