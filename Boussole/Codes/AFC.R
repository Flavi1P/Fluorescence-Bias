library(tidyverse)
library(FactoMineR)

pigments <- c("fuco", "peri", "hex", "but", "allo", "tchlb", "zea")

boussole <- read_csv("Boussole/Data/boussole.csv")

boussole_ca <- boussole %>% mutate(pigsum = rowSums(select(boussole, pigments))) %>% filter(pigsum > 0 & ratio < 10) %>% select(pigments, ratio, micro, nano, pico)
boussole_ca <- na.omit(boussole_ca)


AFC <- cca(select(boussole_ca, pigments))

scores <- data.frame(scores(AFC, choices = c(1,2,3), display = "site"))
boussole_ca <- bind_cols(boussole_ca, scores)

pigscore <- data.frame(scores(AFC, choices = c(1,2,3), display = "species"))


fitscore <- envfit(AFC, select(boussole_ca, micro, nano, pico, ratio))
fitarrow <- as.data.frame(fitscore$vectors$arrows)

ggplot(boussole_ca)+
  geom_point(aes(x = CA1, y = CA2, colour = nano))+
  geom_segment(aes(x = 0, xend = CA1, y = 0, yend = CA2), data = pigscore)+
  geom_text(aes(x = CA1, y = CA2, label = rownames(pigscore)), data = pigscore)+
  geom_segment(aes(x = 0, y = 0, xend = CA1*1.7, yend = CA2*1.7), data = fitarrow, colour = "#33a02c")+
  geom_text(aes(x = CA1*1.7, y = CA2*1.7, label=rownames(fitarrow), fontface = 2), data = fitarrow)+
  scale_color_viridis_c()+
  xlab("CA1 39%")+
  ylab("CA2 22%")+
  ggtitle("Correspondance Analysis of boussole pigment (2013:2015)")+
  coord_equal()

boussole$optical_layer <- round(boussole$optical_layer)
boussole <- filter(boussole, optical_layer < 4)
summary(lm(ratio~(fuco + peri + allo + but + hex + zea+ tchlb)*press, data = boussole))
