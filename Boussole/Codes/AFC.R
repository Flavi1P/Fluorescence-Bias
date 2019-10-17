library(tidyverse)
library(FactoMineR)
library(lubridate)
library(zoo)
library(castr)
source("functions/outliers.R")

pigments <- c("fuco", "peri", "hex", "but", "allo", "tchlb", "zea")

boussole <- read_csv("Boussole/Data/boussole.csv")

influential <- outliers(select(boussole, fluo, tchla))

boussole <- boussole[-influential,]

boussole_ca <- boussole %>% mutate(pigsum = rowSums(select(boussole, pigments))) %>% filter(pigsum > 0 & ratio < 10) %>% select(pigments, ratio, micro, nano, pico, date, depth)
boussole_ca <- na.omit(boussole_ca)

boussole_ca$month <- month(boussole_ca$date, label = TRUE)
boussole_ca$month <- substr(boussole_ca$month, 0, 3)
boussole_ca$year <- year(boussole_ca$date)

boussole_ts <- filter(boussole_ca, depth < 10 & ratio <6)
smooth_ratio <- boussole_ts %>% group_by(date) %>% select(date,ratio) %>% summarize_all(mean) %>% ungroup()
names(smooth_ratio) <- c("date", "smooth")

boussole_ts <- left_join(boussole_ts, smooth_ratio)
boussole_ts$smooth_ratio <- smooth(boussole_ts$ratio, k = 2, n = 2)

ggplot(boussole_ts)+
  geom_point(aes(x = date, y = ratio))+
  geom_smooth(aes(x = date, y = ratio), se = FALSE)+
  facet_wrap(facets = "year", ncol = 1, scales = "free_x")+
  ylab("rapport [ chla]fluo/[chla]hplc")+
  theme_bw(base_size = 20)
  
#ggsave("Boussole/Plots/ts.png", scale = 2)

library(ggtern)

seasons <- data.frame("saison" = rep(c("winter","spring","summer","fall"), each = 3), "month" = c("jan", "fév", "mar", "avr", "mai", "jui", "jui", "aoû", "sep", "oct", "nov", "déc"))

boussole_ts <- left_join(boussole_ts, seasons)

ggplot(boussole_ts)+
  coord_tern()+
  stat_interpolate_tern(geom = "polygon", formula = value~x+y,
                        method = lm, n = 50,
                        breaks = seq(0,4, by = 0.5),
                        aes(x = micro, y = nano, z = pico, value = ratio, fill =..level..), expand = 1)+
  geom_point(aes(x= micro, y = nano, z = pico, colour = saison), size = 3)+
  scale_color_brewer(palette = "Set1", name = "Saison")+
  theme_bw(base_size = 20)+
  weight_percent()+
  scale_fill_gradient(name = "Rapport Fluo/Chla", low = "lightgrey", high = "gray36")

ggplot(boussole_ts)+
  geom_smooth(aes(x = date, y = micro, colour = "micro"), se = FALSE)+
  geom_smooth(aes(x = date, y = nano, colour = "nano"), se = FALSE)+
  geom_smooth(aes(x = date, y = pico, colour = "pico"), se = FALSE)+
  facet_wrap(facets = "year", ncol = 1, scales = "free_x")+
  ylab("fmicro")+
  theme_bw(base_size = 20)

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
