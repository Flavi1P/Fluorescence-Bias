library(tidyverse)
library(Metrics)
library(vegan)
library(nnls)
library(FactoMineR)
source("functions/phi_lm.R")
source("functions/outliers.R")
source("functions/phi_simple.R")
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

lines <- sample(1:nrow(merged_argo), 200)
merged_argo_train <- merged_argo[lines,]
merged_argo_test <- filter(merged_argo, !(rownames(merged_argo) %in% lines))

#phi####

phi_argo <- phi_simple(merged_argo_train, variable = "fluo")
phi_argo$se <- ifelse(phi_argo$se > phi_argo$phi, phi_argo$phi, phi_argo$se)

phi_argo_tchla <- phi_simple(merged_argo_train, variable = "tchla")

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

argo_calibration_test <- argo_calibration_test %>% mutate(predict_fluo = tchla*(micro* phi_micro + nano * phi_nano + pico * phi_pico),
                                                calibrate_fluo = (fluo/(micro * phi_micro + nano * phi_nano + pico * phi_pico))) 

ggplot(argo_calibration_test)+
  geom_point(aes(x = tchla , y = calibrate_fluo, colour = "calibrate"))+
  geom_point(aes(x = tchla , y = chla_adjusted, colour = "fluo"))+
  geom_line(aes(x = tchla, y = tchla))

argo_calibration_test <- filter(argo_calibration_test, is.na(calibrate_fluo) == FALSE)
a <- rmse(argo_calibration_test$tchla, argo_calibration_test$calibrate_fluo)
b <- rmse(argo_calibration_test$tchla, argo_calibration_test$chla_adjusted)
a/b

argo_calibration_test$ratio <- argo_calibration_test$calibrate_fluo/argo_calibration_test$tchla

rmse_value <- c(1:10000)
rmse_commu <- c(1:10000)
rmse_lineaire <- c(1:10000)

for(i in 1:10000){
  lines <- sample(1:nrow(merged_argo), 200)
  merged_argo_train <- merged_argo[lines,]
  merged_argo_test <- filter(merged_argo, !(rownames(merged_argo) %in% lines))
  
  #phi####
  
  phi_argo <- phi_simple(merged_argo_train, variable = "fluo")
  phi_argo$se <- ifelse(phi_argo$se > phi_argo$phi, phi_argo$phi, phi_argo$se)
  
  
  phi_argo <- phi_argo %>% select(phi, optical_layer, size) %>% spread(key = size, value = phi)
  names(phi_argo) <- c("optical_layer", "phi_micro", "phi_nano", "phi_pico") 
  
  
  argo_calibration_test <- left_join(merged_argo_test, phi_argo, by = "optical_layer") 
  
  argo_calibration_test <- argo_calibration_test %>% mutate(predict_fluo = tchla*(micro* phi_micro + nano * phi_nano + pico * phi_pico),
                                                            calibrate_fluo = (fluo/(micro * phi_micro + nano * phi_nano + pico * phi_pico))) 
  
  argo_calibration_test <- filter(argo_calibration_test, is.na(calibrate_fluo) == FALSE)
  a <- rmse(argo_calibration_test$tchla, argo_calibration_test$calibrate_fluo)
  b <- rmse(argo_calibration_test$tchla, argo_calibration_test$chla_adjusted)
  a/b
  
  argo_calibration_test$ratio <- argo_calibration_test$calibrate_fluo/argo_calibration_test$tchla
  print(i)
  
  rmse_value[i] <- a/b
  rmse_commu[i] <- a
  rmse_lineaire[i] <- b
}

mean(rmse_value)
median(rmse_value)
sd(rmse_value)

mean(rmse_lineaire)
median(rmse_lineaire)
sd(rmse_lineaire)

mean(rmse_commu)
median(rmse_commu)
sd(rmse_commu)
