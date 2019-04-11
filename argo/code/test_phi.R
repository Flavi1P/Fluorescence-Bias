library(tidyverse)
library(Metrics)
library(vegan)
library(nnls)
library(FactoMineR)
source("functions/phi_lm.R")
source("functions/outliers.R")
source("functions/phi_boot.R")
pigments <- c("fuco", "peri", "hex", "but", "allo", "tchlb", "zea")
NAT_IRS_list <- c("lovbio059c", "lovbio045b", "lovbio024c", "lovbio044b", "lovbio031c", "lovbio027b", "lovbio040b", "lovbio026c")

map_vec <- read_csv("Data/map_vec")
merged_argo <- read_csv("Data/merged_argo")
merged_argo <- filter(merged_argo, lovbio != "lovbio067c" & lovbio != "lovbio083d" & lovbio != "lovbio085d" & lovbio != "lovbio090d") #filter profile with only 1 match
merged_argo <- filter(merged_argo, !(lovbio %in% NAT_IRS_list))

merged_argo$fluo <- merged_argo$chla_adjusted * 2
merged_argo$ratio <- merged_argo$fluo/merged_argo$tchla

influential <- outliers(select(merged_argo, fluo, tchla, micro, nano, pico))

merged_argo <- merged_argo[-influential,] %>% filter(optical_layer<4)

lines <- sample(1:233, 200)
merged_argo_train <- merged_argo[lines,]
merged_argo_test <- filter(merged_argo, !(rownames(merged_argo) %in% lines))

#phi####

phi_argo <- phi_boot(merged_argo_train, variable = "fluo")
phi_argo$se <- ifelse(phi_argo$se > phi_argo$phi, phi_argo$phi, phi_argo$se)

phi_argo_tchla <- phi_boot(merged_argo_train, variable = "tchla")

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
  ggtitle("Fluorescent yield ratio, argo")


phi_argo <- phi_argo %>% select(phi, optical_layer, size) %>% spread(key = size, value = phi)
names(phi_argo) <- c("optical_layer", "phi_micro", "phi_nano", "phi_pico") 


argo_calibration_test <- left_join(merged_argo_test, phi_argo) 

argo_calibration_test <- argo_calibration_test %>% mutate(predict_fluo = tchla*(micro*amicro * phi_micro + nano*anano * phi_nano + pico*apico * phi_pico),
                                                calibrate_fluo = (fluo/(micro*amicro * phi_micro + nano*anano * phi_nano + pico*apico * phi_pico))) 

ggplot(argo_calibration_test)+
  geom_point(aes(x = tchla , y = calibrate_fluo, colour = "calibrate"))+
  geom_point(aes(x = tchla , y = chla_adjusted, colour = "fluo"))+
  geom_line(aes(x = tchla, y = tchla))

argo_calibration_test <- filter(argo_calibration_test, is.na(calibrate_fluo) == FALSE)
a <- rmse(argo_calibration_test$tchla, argo_calibration_test$calibrate_fluo)
b <- rmse(argo_calibration_test$tchla, argo_calibration_test$chla_adjusted)
a/b

argo_calibration_test$ratio <- argo_calibration_test$calibrate_fluo/argo_calibration_test$tchla
#afc####

afc_table <- na.omit(select(argo_calibration_test, pigments, micro, nano, pico,  ratio))

afc_argo <- cca(na.omit(select(afc_table, pigments)))
test <- envfit(afc_argo, select(afc_table, micro, nano, pico, ratio))

env_arrow <- as.data.frame(test$vectors$arrows) #On récupère les donnes du modele sur les variables environnementales

argo_score <- as.data.frame(scores(afc_argo, choices = c(1,2,3,4,5), display = c("site"))) #this is a dataframe with the score on the 5 axes of each sample

afc_table <- bind_cols(afc_table, argo_score)

pig_score <- as.data.frame(scores(afc_argo, choices = c(1,2,3,4,5), display = c("species")))


ggplot()+
  geom_point(aes(x = CA1, y = CA2, colour = ratio), size = 2, data = afc_table)+
  geom_segment(aes(x = 0, xend = CA1 *1.5, y = 0, yend = CA2*1.5), data = pig_score)+
  geom_text_repel(aes(x = CA1*1.5, y = CA2*1.5, label = rownames(pig_score)), data = pig_score)+
  geom_segment(aes(x = 0, y = 0, xend = CA1*1.7, yend = CA2*1.7), data = env_arrow, colour = "#33a02c")+
  geom_text(aes(x = CA1*1.7, y = CA2*1.7, label=rownames(env_arrow), fontface = 2), data = env_arrow)+
  coord_equal()

