library(tidyverse)
library(Metrics)
library(vegan)
library(nnls)
library(FactoMineR)
source("functions/phi_lm.R")
source("functions/outliers.R")
source("functions/phi_simple.R")
source("functions/normalize.R")

pigments <- c("fuco", "peri", "hex", "but", "allo", "tchlb", "zea")
NAT_IRS_list <- c("lovbio059c", "lovbio045b", "lovbio024c", "lovbio044b", "lovbio031c", "lovbio027b", "lovbio040b", "lovbio026c")

map_vec <- read_csv("Data/map_vec")
merged_argo <- read_csv("Data/merged_argo")
merged_argo <- filter(merged_argo, lovbio != "lovbio067c" & lovbio != "lovbio083d" & lovbio != "lovbio085d" & lovbio != "lovbio090d") #filter profile with only 1 match
merged_argo <- filter(merged_argo, !(lovbio %in% NAT_IRS_list))


#phi####
merged_argo$fluo <- merged_argo$chla_adjusted * 2
merged_argo$ratio <- merged_argo$fluo/merged_argo$tchla

influential <- outliers(select(merged_argo, fluo, tchla, ratio, micro, nano, pico))

merged_argo <- merged_argo[-influential,]

ggplot(merged_argo)+
  geom_density(aes(x = ratio))
merged_argo <- filter(merged_argo, optical_layer < 4)

phi_argo <- phi_simple(merged_argo, variable = "fluo")
write_csv(phi_argo, "argo/Data/phi_argo")
phi_argo$phi <- ifelse(phi_argo$size == "pico", 0, phi_argo$phi)
phi_argo$se <- ifelse(phi_argo$se > phi_argo$phi, phi_argo$phi, phi_argo$se)



ggplot(phi_argo, aes(x=size, y = phi, fill = as.factor(optical_layer))) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin = phi-se, ymax = phi+se),position=position_dodge())+
  scale_fill_viridis_d( name = "Couche optique")+
  ylab("Phi")+ xlab("classe de taille") +
  theme_bw(base_size = 18)

ggsave("argo/Plots/phi_argo.png")

  phi_argo <- phi_argo %>% select(phi, optical_layer, size) %>% spread(key = size, value = phi)
names(phi_argo) <- c("optical_layer", "phi_micro", "phi_nano", "phi_pico") 


argo_calibration <- left_join(merged_argo, phi_argo) 

argo_calibration <- argo_calibration %>% mutate(predict_fluo = tchla*(micro*amicro * phi_micro + nano*anano * phi_nano + pico*apico * phi_pico),
                                                calibrate_fluo = (fluo/(micro* phi_micro + nano * phi_nano + pico * phi_pico))) 

ggplot(argo_calibration)+
  geom_point(aes(x = tchla , y = calibrate_fluo, colour = "calibrate"))+
  geom_point(aes(x = tchla , y = chla_adjusted, colour = "fluo"))+
  geom_line(aes(x = tchla, y = tchla))+
  scale_colour_brewer(palette = "Set1")+
  theme_bw()
argo_calibration <- filter(argo_calibration, is.na(calibrate_fluo) == FALSE)
a <- rmse(argo_calibration$tchla, argo_calibration$calibrate_fluo)
b <- rmse(argo_calibration$tchla, argo_calibration$chla_adjusted)
a/b
