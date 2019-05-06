library(tidyverse)
library(FactoMineR)
library(lubridate)

pigments <- c("fuco", "peri", "hex", "but", "allo", "tchlb", "zea")

boussole <- read_csv("Boussole/Data/boussole.csv")

boussole_ca <- boussole %>% mutate(pigsum = rowSums(select(boussole, pigments))) %>% filter(pigsum > 0 & ratio < 10) %>% select(pigments, ratio, micro, nano, pico, date)
boussole_ca <- na.omit(boussole_ca)

boussole_ca$month <- month(boussole_ca$date, label = TRUE, abbr = FALSE)
boussole_ca$month <- substr(boussole_ca$month, 0, 3)


AFC <- cca(select(boussole_ca, pigments))

scores <- data.frame(scores(AFC, choices = c(1,2,3), display = "site"))
boussole_ca <- bind_cols(boussole_ca, scores)

pigscore <- data.frame(scores(AFC, choices = c(1,2,3), display = "species"))


fitscore <- envfit(AFC, select(boussole_ca, micro, nano, pico, ratio))
fitarrow <- as.data.frame(fitscore$vectors$arrows)

ggplot(boussole_ca)+
  geom_point(aes(x = CA1, y = CA2, colour = month), size = 1.5)+
  geom_segment(aes(x = 0, xend = CA1*2, y = 0, yend = CA2*2), data = pigscore)+
  geom_text_repel(aes(x = CA1*2, y = CA2*2, label = rownames(pigscore)), data = pigscore, size = 6)+
  geom_segment(aes(x = 0, y = 0, xend = CA1*3, yend = CA2*3), data = fitarrow, colour = "#33a02c")+
  geom_text(aes(x = CA1*3, y = CA2*3, label=rownames(fitarrow), fontface = 2), data = fitarrow, size = 6)+
  scale_color_viridis_d(name = "mois")+
  xlab("CA1 39%")+
  ylab("CA2 22%")+
  ggtitle("Boussole")+
  #ggtitle("Correspondance Analysis of boussole pigment (2013:2015)")+
  ylim(-4,2.5)+
  coord_equal()+
  theme_bw(base_size = 18)
ggsave("Boussole/plots/afc_boussole.png", scale = 1)

boussole$optical_layer <- round(boussole$optical_layer)
boussole <- filter(boussole, optical_layer < 4)
summary(lm(ratio~(fuco + peri + allo + but + hex + zea+ tchlb)*press, data = boussole))
