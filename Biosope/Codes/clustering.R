library(tidyverse)
library(zoo)
library(vegan)

source("functions/profile_numb.R")
source("functions/zeu_moma.R")

pigments <- c("fuco", "peri", "hex", "but", "allo", "tchlb", "zea")

biosope <- read_csv("Data/Final/biomerged")
biosope$depth <- abs(biosope$depth)

biosope$profile <- profile_numb(biosope$depth, "upward")

table(biosope$profile)

momadata <- filter(biosope, profile %in% which(table(biosope$profile) > 1)) %>%  select(depth, tchla, profile)
momadata$ze <- NA
for (i in momadata$profile){
  dat <- filter(momadata, profile == i)
  momadata$ze[which(momadata$profile == i)] <- Zeu_moma(dat$tchla, dat$depth)
}

momadata <- momadata %>% group_by(profile) %>% summarize_all(mean)
biosope <- left_join(biosope, select(momadata, ze, profile), by ="profile")


biosope$ze  = approx(x = momadata$profile, y = momadata$ze, xout = biosope$profile)$y #linear interpolmation of ze depth through profiles to avoid NA


biosope$ze <- rollmean(biosope$ze, 20, fill = c(65.9, NA, 43.9))
biosope$ze <- rollmean(biosope$ze, 20, fill = c(65.9, NA, 43.9))
biosope$ze <- rollmean(biosope$ze, 20, fill = c(65.9, NA, 43.9))
biosope$ze <- rollmean(biosope$ze, 20, fill = c(65.9, NA, 43.9))



ggplot(biosope)+
  geom_point(aes(x = lon, y = -depth))+
  geom_path(aes(x = lon, y = -ze))


biosope <- biosope %>% mutate(wdp = 1.56 * fuco + 0.92 * peri + 4.8 * allo + 1.02 * but + 1.12 * hex + 1.51 * zea + 0.69 * tchlb,
                              micro = (1.56 * fuco + 0.92 * peri)/wdp,
                              nano = (4.8 * allo + 1.02 * but + 1.51 * hex)/wdp,
                              pico = (1.51 * zea + 0.69 * tchlb)/wdp,
                              amicro = 0.0164 * exp(0.79*(depth/ze)),
                              anano = 0.0876 * exp(-0.45*(depth/ze)),
                              apico = 0.1393 * exp(-0.69*(depth/ze)),
                              microfluo = amicro * micro * tchla,
                              nanofluo = anano * nano * tchla,
                              picofluo = apico * pico * tchla,
                              optical_layer = (depth/ze)/0.50001 + 1,
                              ratio = fluo_urel / tchla) %>% filter(tchla > 0.02)

ggplot(biosope)+
  geom_point(aes(x = lon, y = -depth, colour = ratio))+
  scale_color_viridis_c()

biosope$pigsum <- rowSums(select(biosope, pigments))
biosope <-  filter(biosope, pigsum > 0 & fluo_urel != "NA" & optical_layer < 4)

mod <- lm(tchla~fluo_urel + micro + nano + pico, data = biosope)
cooksd <- cooks.distance(mod)

plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")
abline(h = 4*mean(cooksd, na.rm=T), col="red") 
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")

influential <- as.numeric(names(cooksd)[(cooksd > 4*mean(cooksd, na.rm=T))]) 
biosope <- biosope[-influential,]

AFC <- cca(select(biosope, pigments))

scores <- data.frame(scores(AFC, choices = c(1,2,3), display = "site"))
biosope <- bind_cols(biosope, scores)

pigscore <- data.frame(scores(AFC, choices = c(1,2,3), display = "species"))



fitscore <- envfit(AFC, select(biosope, micro, nano, pico, ratio))
fitarrow <- as.data.frame(fitscore$vectors$arrows)

ggplot(biosope)+
  geom_point(aes(x = CA1, y = CA2, colour = pico))+
  geom_segment(aes(x = 0, xend = CA1, y = 0, yend = CA2), data = pigscore)+
  geom_text(aes(x = CA1, y = CA2, label = rownames(pigscore)), data = pigscore)+
  geom_segment(aes(x = 0, y = 0, xend = CA1*1.7, yend = CA2*1.7), data = fitarrow, colour = "#33a02c")+
  geom_text(aes(x = CA1*1.7, y = CA2*1.7, label=rownames(fitarrow), fontface = 2), data = fitarrow)+
  scale_color_viridis_c()

distbio <- dist(select(biosope, CA1, CA2))
biosope$group <- as.factor(cutree(hclust(distbio, method = "ward.D"), k = 4))

ggplot(biosope)+
  geom_point(aes(x = CA1, y = CA2, colour = group), size = 1.5)+
  geom_segment(aes(x = 0, xend = CA1, y = 0, yend = CA2), data = pigscore)+
  geom_text(aes(x = CA1, y = CA2, label = rownames(pigscore)), data = pigscore)+
  geom_segment(aes(x = 0, y = 0, xend = CA1*1.7, yend = CA2*1.7), data = fitarrow, colour = "#33a02c")+
  geom_text(aes(x = CA1*1.7, y = CA2*1.7, label=rownames(fitarrow), fontface = 2), data = fitarrow)+
  xlab("CA1 (63%)")+ ylab("CA2 (22%)")+
  theme_bw(base_size = 14)+
  coord_equal()+
  scale_color_viridis_d()

ggsave("Biosope/Plots/afc_biosope.png", scale = 1)

resume_clust <- biosope %>% select(group, micro, nano, pico, ratio, tchla, fluo_urel) %>% group_by(group) %>% 
  summarize_all(c(mean, sd)) %>% ungroup()
names(resume_clust) <- c("group","micro_mean", "nano_mean", "pico_mean", "ratio_mean", "tchla_mean", "fluo_mean",
                         "micro_sd", "nano_sd", "pico_sd", "ratio_sd", "tchla_sd", "fluo_sd")

ggplot(resume_clust)+
  geom_col(aes(x = group, y = ratio_mean))+
  geom_errorbar(aes(x = group, ymin = ratio_mean - ratio_sd, ymax = ratio_mean + ratio_sd))

tall_cluster <- gather(resume_clust, key = group, "mean")
tall_cluster$cluster = c(1,2,3,4)

ggplot(filter(tall_cluster, group %in% c("micro_mean", "nano_mean", "pico_mean")))+
  geom_col(aes(x = cluster, y = mean, fill = group),position = "fill")

AFC_detrend <- decorana(select(biosope, pigments))
AFC_detrend

plot(AFC_detrend)
scores_detrend <- data.frame(scores(AFC_detrend, choices = c(1,2,3), display = "site"))
biosope <- bind_cols(biosope, scores_detrend)

pigscore_detrend <- data.frame(scores(AFC_detrend, choices = c(1,2,3), display = "species"))

fitscore_detrend <- envfit(AFC_detrend, select(biosope, micro, nano, pico, ratio, tchla))
fitarrow_detrend <- as.data.frame(fitscore_detrend$vectors$arrows)

distbio_detrend <- dist(select(biosope, DCA1, DCA2))
biosope$group_detrend <- as.factor(cutree(hclust(distbio_detrend, method = "ward.D"), k = 4))

plot(hclust(distbio_detrend, method = "ward.D"))
ggplot(biosope)+
  geom_point(aes(x = DCA1, y = DCA2, colour = group_detrend))+
  geom_segment(aes(x = 0, xend = DCA1, y = 0, yend = DCA2), data = pigscore_detrend)+
  geom_text(aes(x = DCA1, y = DCA2, label = rownames(pigscore_detrend)), data = pigscore_detrend)+
  geom_segment(aes(x = 0, y = 0, xend = DCA1, yend = DCA2), data = fitarrow_detrend, colour = "#33a02c")+
  geom_text(aes(x = DCA1, y = DCA2, label=rownames(fitarrow_detrend), fontface = 2), data = fitarrow_detrend)+
  scale_color_viridis_d(name = "cluster")+
  coord_equal()+
  xlab("DCA1 52%")+ylab("DCA2 15%")+
  ggtitle("Detrend Correspondance analysis on Biosope HPLC data")

library(ggsci)
ggplot(biosope)+
  geom_point(aes(x = lon, y = -depth, colour = group_detrend), size = 4)+
  geom_path(aes(x = lon, y = -ze), colour = "Black")+
  ylab("Profondeur")+
  scale_color_futurama(name = "Environnement", labels = c("Profond", "Upwelling", "DCM", "Surface"))+
  xlab("longitude")+
  theme_bw(base_size = 20)

ggsave("Biosope/Plots/transect.png", scale = 2)

world_data <- map_data("world") %>% fortify()


ggplot(biosope)+
  geom_point(aes(x = lon, y = lat))+
  geom_polygon(aes(x = long, y = lat, group = group), data = world_data)+
  coord_quickmap(xlim = c(-160, -50), ylim = c(-60, 0))+
  xlab("longitude")+ylab("latitude")+
  theme_bw(base_size = 20)

ggsave("Biosope/Plots/cruise_map.png")

summary(lm(ratio~(microfluo+nanofluo+picofluo) * optical_layer, data = biosope))

#write_csv(biosope, "Biosope/Data/biosope")
# resume_clust_detrend <- biosope %>% select(group_detrend, micro, nano, pico, ratio, tchla, fluo_urel) %>% group_by(group_detrend) %>% 
#   summarize_all(c(mean, sd)) %>% ungroup()
# names(resume_clust_detrend) <- c("group","micro_mean", "nano_mean", "pico_mean", "ratio_mean", "tchla_mean", "fluo_mean",
#                          "micro_sd", "nano_sd", "pico_sd", "ratio_sd", "tchla_sd", "fluo_sd")
# 
# ggplot(resume_clust_detrend)+
#   geom_col(aes(x = reorder(group, ratio_mean), y = ratio_mean))+
#   geom_errorbar(aes(x = reorder(group, ratio_mean), ymin = ratio_mean - ratio_sd, ymax = ratio_mean + ratio_sd))
# 
# tall_cluster <- gather(resume_clust, key = group, "mean")
# tall_cluster$cluster = c(1,2,3)
# 
# ggplot(filter(tall_cluster, group %in% c("micro_mean", "nano_mean", "pico_mean")))+
#   geom_col(aes(x = cluster, y = mean, fill = group),position = "fill")

