library(tidyverse)
library(Metrics)
library(vegan)
library(nnls)
library(sf)
library(FactoMineR)
library(readxl)
library(janitor)
source("functions/phi_lm.R")
source("functions/outliers.R")
source("functions/phi_boot.R")
path = "Data/Longhurst"
pigments <- c("fuco", "peri", "hex", "but", "allo", "tchlb", "zea")
NAT_IRS_list <- c("lovbio059c", "lovbio045b", "lovbio024c", "lovbio044b", "lovbio031c", "lovbio027b", "lovbio040b", "lovbio026c")

map_vec <- read_csv("Data/map_vec")
merged_argo <- read_csv("Data/merged_argo")
merged_argo <- filter(merged_argo, lovbio != "lovbio067c" & lovbio != "lovbio083d" & lovbio != "lovbio085d" & lovbio != "lovbio090d") #filter profile with only 1 match
merged_argo <- filter(merged_argo, !(lovbio %in% NAT_IRS_list))


#a_sol_ph calculation####

spectre <- read_excel("Biosope/Data/Spectres_annick.xlsx")
spectre <- clean_names(spectre)

spectre470 <- filter(spectre, lambda == 470)

merged_argo <- merged_argo %>% mutate(abs_470 = peri * spectre470$peri + but * spectre470$x19_bf + hex * spectre470$x19_hf + fuco * spectre470$fuco + allo * spectre470$allox + chla * spectre470$chl_a + spectre470$zea * zea + spectre470$chl_b * tchlb,
                        protect_470 =  spectre470$zea * zea,
                        a_ps_470 = abs_470 - protect_470,
                        fluo = chl_smooth * 2,
                        ratio = fluo/tchla,
                        phi_app = ratio/abs_470) %>% filter(ratio < 100)

ggplot(merged_argo)+
  geom_point(aes(y = chl_smooth, x = a_ps_470))+
  geom_point(aes(x = a_ps_470, y = tchla, colour = 'green'))

table_pig <- select(merged_argo, phi_app, pigments, tchla, fluo, ratio, abs_470, lovbio) %>%
  mutate(pigsum = rowSums(select(., pigments), na.rm = TRUE)) %>% 
  pivot_longer(., pigments, names_to = 'pigment', values_to = 'concentration')

ggplot(filter(table_pig, concentration > 0 & phi_app > 0 & phi_app < 5000))+
  geom_point(aes(y  = phi_app, x = concentration, colour = pigment), method = 'lm')+
  coord_trans(x = 'log', y = 'log')+
  facet_wrap(.~ pigment)

table_pig <- filter(table_pig,concentration > 0 & phi_app > 0 & phi_app < 5000) %>% 
  mutate(phi_app = phi_app/max(.$phi_app))

ggplot(table_pig)+
  geom_point(aes(y = phi_app, x = tchla, colour = ratio))+
  coord_trans(x = 'log', y = 'log')+
  scale_color_viridis_c(name = 'fluo/chla')


table_pig <- mutate(table_pig, aspe = (abs_470/tchla),
                    phi_app2 = ratio/aspe)

ggplot(table_pig)+
  geom_point(aes(y = aspe, x = tchla))+
  coord_trans(x = 'log', y = 'log')

ggplot(table_pig)+
  geom_point(aes(y = phi_app2, x = tchla, colour = ratio))+
  coord_trans(x = 'log', y = 'log')+
  scale_color_viridis_c(name = 'fluo/chla')


ggplot(merged_argo)+
  geom_point(aes(x = a_ps_470, y = ratio, colour = lovbio), colour = )

ggplot(problem)+
  geom_point(aes(x = a_ps_470, y = ratio, colour = lovbio))+
  geom_smooth(aes(x = a_ps_470, y = ratio), method = 'lm')

summary(lm(a_ps_470~ratio, data = problem))

ggplot(problem)+
  geom_point(aes(y = phi_app, x = hex))

#phi####


influential <- outliers(select(merged_argo, fluo, tchla, ratio, micro, nano, pico))

merged_argo <- merged_argo[-influential,]

ggplot(merged_argo)+
  geom_density(aes(x = ratio))
merged_argo <- filter(merged_argo, optical_layer < 4)

phi_argo <- phi_boot(merged_argo, variable = "fluo")
phi_argo$se <- ifelse(phi_argo$se > phi_argo$phi, phi_argo$phi, phi_argo$se)
#write_csv(phi_argo, "argo/Data/phi_argo")


phi_argo_tchla <- phi_boot(merged_argo, variable = "tchla")

yield_ratio <- phi_argo$phi/phi_argo_tchla$phi

phi_argo$yield_ratio <- yield_ratio

ggplot(phi_argo, aes(x=size, y = phi, fill = as.factor(optical_layer))) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin = phi-se, ymax = phi+se),position=position_dodge())+
  scale_fill_viridis_d( name = "optical layer")+
  ylab("Phi")+ xlab("size classe")

ggplot(phi_argo, aes(x=size, y = yield_ratio, fill = as.factor(optical_layer))) +
  geom_bar(position=position_dodge(), stat="identity") +
  scale_fill_viridis_d( name = "optical layer")+
  ylab("Yield ratio")+ xlab("size classe")+
  geom_errorbar(aes(size, ymax = 1, ymin = 1),
                size=0.5, linetype = "longdash", inherit.aes = F, width = 1)+
  ggtitle("Fluorescent yield ratio, argo")


phi_argo <- phi_argo %>% select(phi, optical_layer, size) %>% spread(key = size, value = phi)
names(phi_argo) <- c("optical_layer", "phi_micro", "phi_nano", "phi_pico") 


argo_calibration <- left_join(merged_argo, phi_argo) 

argo_calibration <- argo_calibration %>% mutate(predict_fluo = tchla*(micro*amicro * phi_micro + nano*anano * phi_nano + pico*apico * phi_pico),
                                                calibrate_fluo = (fluo/(micro*amicro * phi_micro + nano*anano * phi_nano + pico*apico * phi_pico))) 

ggplot(argo_calibration)+
  geom_violin(aes(y = chla_adjusted/tchla, x = optical_layer))
longhurst_sf <- read_sf(dsn = path.expand(path), quiet = TRUE)

names(longhurst_sf) <- c("code", "region", "geometry")


pnts_sf <- do.call("st_sfc",c(lapply(1:nrow(argo_calibration),
                                     function(i) {st_point(as.numeric(argo_calibration[i,c("lon.y", "lat.y") ]))}), list("crs" = 4326))) 
pnts_trans <- st_transform(pnts_sf, 4326)
longhurst_trans <- st_transform(longhurst_sf, 4326)  
argo_calibration$code <- apply(st_intersects(longhurst_trans, pnts_trans, sparse = FALSE), 2, 
                   function(col) { 
                     longhurst_trans[which(col), ]$code
                   })

argo_mean <- argo_calibration %>% mutate(ratio = calibrate_fluo/tchla, ratio2 = fluo/tchla) %>% select(code, ratio, ratio2) %>% 
  group_by(code) %>% 
  na.omit(.) %>% 
  summarise_all(c(mean, sd))

argo_mean$ratio2_fn1[argo_mean$code == "BPLR"] <- argo_mean$ratio2_fn1[argo_mean$code == "BPLR"] +0.5
argo_mean$ratio2_fn1[argo_mean$code == "Arct"] <- argo_mean$ratio2_fn1[argo_mean$code == "Arct"] +1


ggplot(argo_mean)+
  geom_col(aes(x = reorder(code, ratio2_fn1), y = ratio_fn1, fill = code))+
  geom_errorbar(aes(x = code, ymin = ratio_fn1-ratio_fn2, ymax = ratio_fn1 + ratio_fn2))+
  scale_fill_brewer(palette = "Set1")+
  geom_errorbar(aes(code, ymax = 1, ymin = 1),
                size= 1, linetype = "longdash", inherit.aes = F, width = 1)+
  guides(fill = FALSE)+
  ylab("Rapport Fluo/Chla")+
  xlab("Province de Longhurst")+
  ylim(0,9)+
  theme_bw(base_size = 20)+
  ggtitle("Corrigé")

ggplot(argo_mean)+
  geom_col(aes(x = reorder(code, ratio2_fn1), y = ratio2_fn1, fill = code))+
  geom_errorbar(aes(x = code, ymin = ratio2_fn1-ratio2_fn2, ymax = ratio2_fn1 + ratio2_fn2))+
  scale_fill_brewer(palette = "Set1")+
  geom_errorbar(aes(code, ymax = 1, ymin = 1),
                size= 1, linetype = "longdash", inherit.aes = F, width = 1)+
  guides(fill = FALSE)+
  ylab("Rapport Fluo/Chla")+
  xlab("Province de Longhurst")+
  ylim(0,9)+
  theme_bw(base_size = 20)+
  ggtitle("Brut")

grid.arrange(g2,g1, ncol = 1)

#ggsave("argo/Plots/histo_calibre.png")

argo_mean_roesler <- argo_calibration %>% mutate(ratio = chla_adjusted/tchla) %>% select(code, ratio) %>% 
  group_by(code) %>% 
  na.omit(.) %>% 
  summarise_all(c(mean, sd))

ggplot(argo_mean_roesler)+
  geom_col(aes(x = reorder(code, fn1), y = fn1, fill = code))+
  geom_errorbar(aes(x = code, ymin = fn1-fn2, ymax = fn1 + fn2))+
  scale_fill_brewer(palette = "Set1")+
  geom_errorbar(aes(code, ymax = 2, ymin = 2),
                size=0.5, linetype = "longdash", inherit.aes = F, width = 1)+
  geom_errorbar(aes(code, ymax = 1, ymin = 1),
                size=0.5, linetype = "longdash", inherit.aes = F, width = 1)+
  guides(fill = FALSE)+
  ylim(0,8)+
  theme_bw()

ggplot(filter(argo_calibration, optical_layer < 4))+
  geom_point(aes(x = tchla , y = calibrate_fluo, colour = "calibrate"), size = 1.5)+
  geom_point(aes(x = tchla , y = chla_adjusted*2, colour = "fluo"), size = 1.5)+
  geom_line(aes(x = tchla, y = tchla))+
  scale_color_brewer(palette = "Set1", name = "", label = c("Chla calibrée", "Chla fluo"))+
  ylab("Concentration en Chla estimée")+xlab("Concentration en Chla mesurée")+
  theme_bw(base_size = 16)+
  coord_trans(x = "log", y = "log")


#ggsave("argo/Plots/calibration.png")
  
g1 <- ggplot(filter(argo_calibration, optical_layer < 4))+
  geom_point(aes(x = tchla , y = chla_adjusted, colour = "fluo"), size = 1.5)+
  geom_smooth(aes(x = tchla, y = chla_adjusted), method = "lm", se = F, colour = "blue")+
  geom_line(aes(x = tchla, y = tchla))+
  ylab("Concentration en Chla estimée")+xlab("concentration en Chla mesurée")+
  theme_bw(base_size = 16)+
  ylim(0,3.5)

g2 <- ggplot(filter(argo_calibration, optical_layer < 4))+
  geom_point(aes(x = tchla , y = calibrate_fluo, colour = "calibrate"), size = 1.5)+
  geom_smooth(aes(x = tchla, y = calibrate_fluo), method = "lm", se = F, colour = "red")+
  geom_line(aes(x = tchla, y = tchla))+
  scale_color_brewer(palette = "Set1", name = "", label = c("Chla calibrée", "Chla fluo"))+
  ylab("Concentration en Chla estimée")+xlab("concentration en Chla mesurée")+
  theme_bw(base_size = 16)+
  ylim(0,3.5)

grid.arrange(g1,g2, ncol = 2)  



argo_calibration <- filter(argo_calibration, is.na(calibrate_fluo) == FALSE)
a <- rmse(argo_calibration$tchla, argo_calibration$calibrate_fluo)
b <- rmse(argo_calibration$tchla, argo_calibration$chla_adjusted)
a/b

library(ggtern)

argo_calibration <- argo_calibration %>% mutate(ratio = chla_adjusted*2/tchla)

ggplot(argo_calibration)+
  coord_tern()+
  stat_interpolate_tern(geom = "polygon", formula = value~x+y,
                        method = lm, n = 50,
                        breaks = seq(0,7, by = 1),
                        aes(x = micro, y = nano, z = pico, value = ratio, fill =..level..), expand = 1)+
  geom_point(aes(x= micro, y = nano, z = pico, colour = code), size = 3)+
  scale_color_brewer(palette = "Set1", name = "Province océanique")+
  theme_bw(base_size = 20)+
  weight_percent()+
  scale_fill_gradient(name = "Rapport Fluo/Chla", low = "lightgrey", high = "gray45")

ggsave("argo/Plots/ggtern.png")
