library(tidyverse)
library(Metrics)
library(vegan)
library(nnls)
library(FactoMineR)
source("functions/phi_lm.R")
source("functions/phi_stat.R")
source("functions/phi_boot.R")
pigments <- c("fuco", "peri", "hex", "but", "allo", "tchlb", "zea")
NAT_IRS_list <- c("lovbio059c", "lovbio045b", "lovbio024c", "lovbio044b", "lovbio031c", "lovbio027b", "lovbio040b", "lovbio026c")

map_vec <- read_csv("Data/map_vec")
merged_argo <- read_csv("Data/merged_argo")
merged_argo <- filter(merged_argo, lovbio != "lovbio067c" & lovbio != "lovbio083d" & lovbio != "lovbio085d" & lovbio != "lovbio090d") #filter profile with only 1 match
merged_argo <- filter(merged_argo, !(lovbio %in% NAT_IRS_list))

ggplot(merged_argo)+
  geom_point(aes(x = lon.x, y = lat.x, colour = "hplc"), size = 3)+
  geom_point(aes(x = lon.y, y = lat.y, colour = "argo"))+
  geom_polygon(aes(x = long, y = lat, group = group), data = map_vec)+
  coord_quickmap()

ggplot(merged_argo)+
  geom_point(aes(x = lon.y, y = lat.y, colour = ze))+
  geom_polygon(aes(x = long, y = lat, group = group), data = map_vec)+
  coord_quickmap()+
  scale_color_viridis_c()

#phi####
merged_argo$fluo <- merged_argo$chla_adjusted * 2
merged_argo$ratio <- merged_argo$fluo/merged_argo$tchla

ggplot(merged_argo)+
  geom_density(aes(x = ratio))
merged_argo <- filter(merged_argo, optical_layer < 4)

phi_argo <- phi_boot(merged_argo, variable = "fluo")
phi_argo$se <- ifelse(phi_argo$se > phi_argo$phi, phi_argo$phi, phi_argo$se)

ggplot(phi_argo, aes(x=size, y = phi, fill = as.factor(optical_layer))) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin = phi-se, ymax = phi+se),position=position_dodge())+
  scale_fill_viridis_d( name = "optical layer")+
  ylab("Phi")+ xlab("size classe")


phi_argo <- phi_argo %>% select(phi, optical_layer, size) %>% spread(key = size, value = phi)
names(phi_argo) <- c("optical_layer", "phi_micro", "phi_nano", "phi_pico") 


argo_calibration <- left_join(merged_argo, phi_argo) 

argo_calibration <- argo_calibration %>% mutate(predict_fluo = tchla*(micro*amicro * phi_micro + nano*anano * phi_nano + pico*apico * phi_pico),
                                                calibrate_fluo = (fluo/(micro*amicro * phi_micro + nano*anano * phi_nano + pico*apico * phi_pico))) 

ggplot(filter(argo_calibration, optical_layer < 4))+
  geom_point(aes(x = tchla , y = calibrate_fluo, colour = "calibrate"))+
  geom_point(aes(x = tchla , y = chla_adjusted, colour = "fluo"))+
  geom_line(aes(x = tchla, y = tchla))

argo_calibration <- filter(argo_calibration, is.na(calibrate_fluo) == FALSE)
a <- rmse(argo_calibration$tchla, argo_calibration$calibrate_fluo)
b <- rmse(argo_calibration$tchla, argo_calibration$chla_adjusted)
a/b

#model ####
t <- data.frame("micro" = NA, "nano" = NA, "pico" = NA, "tchla" = NA)
new_data <- data.frame("micro" = NA, "nano" = NA, "pico" = NA, "tchla" = NA)
for (i in seq(0,1, 0.05)){
  t$micro <- i 
  for(j in seq(0,1, 0.05)){
    tot <- i+j
    if(tot < 1){
      t$nano <- j
      t$pico <- 1-tot
      for(l in seq(0,1.6, 0.2)){
        t$tchla <- l
        new_data <- bind_rows(new_data, t)
      }}
    if (tot == 1){
      t$pico <- 0
      for(l in seq(0,1.6, 0.2)){
        t$tchla <- l
        new_data <- bind_rows(new_data, t)
      }
    }
  }
}
new_data$sum <- rowSums(select(new_data, micro, nano, pico))
new_data <- na.omit(new_data)
new_data <- filter(new_data, sum == 1)

new_data_duplicate <- rbind(new_data, new_data[rep(c(1:1899),2),])
new_data_duplicate$optical_layer <- c(rep(1,1899), rep(2,1899), rep(3,1899))
new_data_phi <- left_join(new_data_duplicate, phi_argo)

new_data_phi <- new_data_phi %>% mutate(zze = optical_layer*0.5,
                                        amicro = 0.0164 * exp(0.79*(zze)),
                                        anano = 0.0876 * exp(-0.45*(zze)),
                                        apico = 0.1393 * exp(-0.69*(zze)),
                                        fluo_predict = tchla * (micro * amicro * phi_micro + nano * anano * phi_nano + pico * phi_pico * apico),
                                        fluo_calibrate = fluo_predict / (micro * amicro * phi_micro + nano * anano * phi_nano + pico * phi_pico * apico))

ggplot(filter(new_data_phi, tchla >0.1 & optical_layer == 1))+
  geom_point(aes(x = tchla, y = fluo_predict, colour = micro))+
  geom_line(aes(x = tchla, y = tchla))+
  scale_colour_viridis_c()

ggplot(filter(new_data_phi, tchla >0.1 & optical_layer == 1))+
  geom_point(aes(x = tchla, y = fluo_predict, colour = nano))+
  geom_line(aes(x = tchla, y = tchla))+
  scale_colour_viridis_c()

ggplot(filter(new_data_phi, tchla >0.1 & optical_layer == 1))+
  geom_point(aes(x = tchla, y = fluo_predict, colour = pico))+
  geom_line(aes(x = tchla, y = tchla))+
  scale_colour_viridis_c()

ggplot(filter(new_data_phi, tchla >0.1 & optical_layer == 2))+
  geom_point(aes(x = tchla, y = fluo_predict, colour = micro))+
  geom_line(aes(x = tchla, y = tchla))+
  scale_colour_viridis_c()

ggplot(filter(new_data_phi, tchla >0.1 & optical_layer == 2))+
  geom_point(aes(x = tchla, y = fluo_predict, colour = nano))+
  geom_line(aes(x = tchla, y = tchla))+
  scale_colour_viridis_c()

ggplot(filter(new_data_phi, tchla >0.1 & optical_layer == 2))+
  geom_point(aes(x = tchla, y = fluo_predict, colour = pico))+
  geom_line(aes(x = tchla, y = tchla))+
  scale_colour_viridis_c()

ggplot(filter(new_data_phi, tchla >0.1 & optical_layer == 3))+
  geom_point(aes(x = tchla, y = fluo_predict, colour = micro))+
  geom_line(aes(x = tchla, y = tchla))+
  scale_colour_viridis_c()

ggplot(filter(new_data_phi, tchla >0.1 & optical_layer == 3))+
  geom_point(aes(x = tchla, y = fluo_predict, colour = nano))+
  geom_line(aes(x = tchla, y = tchla))+
  scale_colour_viridis_c()

ggplot(filter(new_data_phi, tchla >0.1 & optical_layer == 3))+
  geom_point(aes(x = tchla, y = fluo_predict, colour = pico))+
  geom_line(aes(x = tchla, y = tchla))+
  scale_colour_viridis_c()


#AFC ####
biosope <- filter(dataset_filtered, id == "biosope")
merged_argo[is.na(merged_argo)] <- 0
afc_argo <- cca(select(biosope, pigments))
scores <- scores(afc_argo, choices = c(1,2,3), display = "site")

merged_argo <- cbind(biosope, scores)

ggplot(merged_argo)+
  geom_point(aes(x = CA1, y = CA2, colour = nano))+
  scale_color_viridis_c()
dist <- dist(select(merged_argo, starts_with("CA")))
merged_argo$cluster <- as.factor(cutree(hclust(dist, method = "ward.D2"), k = 3))

ggplot(merged_argo)+
  geom_point(aes(x = CA1, y = CA2, colour = cluster))+
  scale_color_viridis_d()
cluster_data <- merged_argo %>% 
  select(micro, nano, pico, tchla, chla, chla_adjusted, calibrate_fluo, cluster, optical_layer) %>% 
  mutate(ratio = chla/tchla) %>% 
  group_by(cluster, optical_layer) %>% 
  summarise_all(funs(mean, sd)) 

ggplot(cluster_data)+
  geom_col(aes(x = cluster, y =ratio_mean))+
  geom_errorbar(aes(x = cluster, ymin = ratio_mean - ratio_sd, ymax = ratio_mean + ratio_sd))+
  facet_grid(vars(optical_layer))

#model fsc####

argo_calibration <- argo_calibration %>% mutate(pico_model = 0.107 * (1 - exp(-6.801 * tchla)),
                                                pico_nano_model = 1.057 * (1 - exp(-0.851 * tchla)),
                                                nano_model = pico_nano_model - pico_model,
                                                micro_model = tchla - (pico_nano_model),
                                                micro_f = micro_model / tchla,
                                                nano_f = micro_model /tchla,
                                                pico_f = pico_model / tchla, 
                                                fluo_model = tchla * (micro_f * amicro * phi_micro + nano_f * anano * phi_nano + pico_f * phi_pico * apico))

for (i in unique(argo_calibration$lovbio)){
  argo_calibration$kd[argo_calibration$lovbio == i] <- 0.0166 + 0.07242 * first(argo_calibration$tchla[argo_calibration$lovbio == i]) * exp(0.68955)
}
                            
ggplot(argo_calibration)+
  geom_density(aes(x= kd))

argo_calibration$tau <- argo_calibration$depth*argo_calibration$kd

ggplot(argo_calibration)+
  geom_density(aes(x= tau))

argo_calibration <- argo_calibration %>% mutate(pico_nano_depth = (0.977 * exp(0.121 * tau))*(1-exp(-(0.910*exp(-0.122*tau))*tchla)),
                                                pico_depth = (0.095 * exp(0.142 * tau))*(1-exp(-(0.910*exp(-7.822*tau))*tchla)),
                                                nano_depth = pico_nano_depth - pico_depth,
                                                micro_depth = tchla - pico_nano_depth,
                                                micro_f_depth = micro_depth / tchla,
                                                nano_f_depth = micro_depth /tchla,
                                                pico_f_depth = pico_depth / tchla, 
                                                fluo_model_depth = tchla * (micro_f_depth * amicro * phi_micro + nano_f_depth * anano * phi_nano + pico_f_depth * phi_pico * apico))

ggplot(argo_calibration)+
  geom_line(aes(x = tchla, y = micro_f, colour = "micro"))+
  geom_line(aes(x = tchla, y = nano_f, colour = "nano"))+
  geom_line(aes(x = tchla, y = pico_f, colour = "pico"))+
  geom_point(aes(x = tchla, y = fluo_model, colour = "fluo", shape = as.factor(optical_layer)))
  geom_point(aes(x = tchla, y = chla), alpha = 3/4, colour = "Grey")

  ggplot(argo_calibration)+
    geom_line(aes(x = tchla, y = pico_f, colour = "micro"))+
    geom_point(aes(x = tchla, y = pico), colour = "Grey", alpha = 3/4)  
  
  ggplot(argo_calibration)+
    geom_density(aes(x = micro_f))
  
  ggplot(argo_calibration)+
    geom_line(aes(x = tchla, y = micro_f_depth, colour = "micro"))+
    geom_line(aes(x= tchla, y = nano_f_depth, colour = "nano"))+
    geom_line(aes(x = tchla, y = pico_f_depth, colour = "pico"))
  
#my model ####
dataset_filtered <- read_csv("Data/Final/data_phi.csv")
biosope <- filter(dataset_filtered, id == "biosope")
  
    
phi_biosope <- phi_lm(biosope, variable = "fluo")
phi_biosope <- phi_biosope %>% select(phi, optical_layer, size) %>% spread(key = size, value = phi)
names(phi_biosope) <- c("optical_layer", "phi_micro", "phi_nano", "phi_pico") 
se_biosope <- phi_lm(biosope, variable = "fluo") %>% select(se, optical_layer, size) %>% spread(key = size, value = se)
names(se_biosope) <- c("optical_layer", "se_micro", "se_nano", "se_pico") 

phi_argo <- phi_boot(merged_argo, variable = "fluo")
phi_argo$se <- ifelse(phi_argo$se > phi_argo$phi, phi_argo$phi, phi_argo$se)
se_argo <- phi_argo %>% select(se, optical_layer, size) %>% spread(key = size, value = se)
names(se_argo) <- c("optical_layer", "se_micro", "se_nano", "se_pico") 
phi_argo <- phi_argo %>% select(phi, optical_layer, size) %>% spread(key = size, value = phi)
names(phi_argo) <- c("optical_layer", "phi_micro", "phi_nano", "phi_pico") 

d_model <- tibble("tchla" = seq(0, 10, 0.05)) %>% mutate(pico_model = 0.107 * (1 - exp(-6.801 * tchla)),
                                                         pico_nano_model = 1.057 * (1 - exp(-0.851 * tchla)),
                                                         nano_model = pico_nano_model - pico_model,
                                                         micro_model = tchla - (pico_nano_model),
                                                         micro_f = micro_model / tchla,
                                                         nano_f = nano_model /tchla,
                                                         pico_f = pico_model / tchla)
ggplot(filter(d_model, tchla >0))+
  geom_line(aes(x = tchla, y = micro_f, colour = "micro"))+
  geom_line(aes(x = tchla, y = nano_f, colour = "nano"))+
  geom_line(aes(x = tchla, y = pico_f, colour = "pico"))+
  coord_trans(x = "log")+scale_x_continuous(breaks = c(0.01, 0.10, 1, 10))+
  theme_bw()

d_model <- bind_rows(d_model, d_model, d_model)

d_model$optical_layer <- c(rep(1, 201), rep(2,201), rep(3,201))
d_model <- na.omit(d_model)
d_model <- left_join(d_model, phi_argo)
d_model <- left_join(d_model, se_argo)

d_model <- d_model %>% mutate(amicro = 0.0164 * exp(0.79*(optical_layer*0.5)),
                              anano = 0.0876 * exp(-0.45*(optical_layer*0.5)),
                              apico = 0.1393 * exp(-0.69*(optical_layer*0.5)),
                              fluo_model = tchla * (micro_f * amicro * phi_micro + nano_f * anano * phi_nano + pico_f * phi_pico * apico),
                              fluo_min = tchla * (micro_f * amicro * (phi_micro - se_micro) + nano_f * anano * (phi_nano - se_nano) + pico_f * (phi_pico - se_pico) * apico),
                              fluo_max = tchla * (micro_f * amicro * (phi_micro + se_micro) + nano_f * anano * (phi_nano + se_nano) + pico_f * (phi_pico + se_pico) * apico),
                              fluo_f = fluo_model/tchla,
                              fluo_f_min = fluo_min/tchla,
                              fluo_f_max = fluo_max/tchla)

ggplot(filter(d_model, tchla >0))+
  geom_path(aes(x = tchla, y = fluo_f, colour = as.factor(optical_layer)))+
  geom_line(aes(x = tchla, y = fluo_f_min, colour = as.factor(optical_layer)), linetype = 2)+
  geom_line(aes(x = tchla, y = fluo_f_max, colour = as.factor(optical_layer)), linetype = 2)+
  coord_trans(x = "log")+scale_x_continuous(breaks = c(0.01, 0.10, 1, 10))+
  scale_color_brewer(palette = "Set1")+
  ylab("f/chla ratio")+ xlab("Chlorophylle a")+
  theme_bw()

b_model <- tibble("tchla" = seq(0, 10, 0.05)) %>% mutate(pico_model = 0.107 * (1 - exp(-6.801 * tchla)),
                                                           pico_nano_model = 1.057 * (1 - exp(-0.851 * tchla)),
                                                           nano_model = pico_nano_model - pico_model,
                                                           micro_model = tchla - (pico_nano_model),
                                                           micro_f = micro_model / tchla,
                                                           nano_f = nano_model /tchla,
                                                           pico_f = pico_model / tchla)
b_model <- bind_rows(b_model, b_model, b_model)

b_model$optical_layer <- c(rep(1, 201), rep(2,201), rep(3,201))
b_model <- na.omit(b_model)
b_model <- left_join(b_model, phi_biosope)
b_model <- left_join(b_model, se_biosope)

b_model <- b_model %>% mutate(amicro = 0.0164 * exp(0.79*(optical_layer*0.5)),
                              anano = 0.0876 * exp(-0.45*(optical_layer*0.5)),
                              apico = 0.1393 * exp(-0.69*(optical_layer*0.5)),
                              fluo_model = tchla * (micro_f * amicro * phi_micro + nano_f * anano * phi_nano + pico_f * phi_pico * apico),
                              fluo_min = tchla * (micro_f * amicro * (phi_micro - se_micro) + nano_f * anano * (phi_nano - se_nano) + pico_f * (phi_pico - se_pico) * apico),
                              fluo_max = tchla * (micro_f * amicro * (phi_micro + se_micro) + nano_f * anano * (phi_nano + se_nano) + pico_f * (phi_pico + se_pico) * apico),
                              fluo_f = fluo_model/tchla,
                              fluo_f_min = fluo_min/tchla,
                              fluo_f_max = fluo_max/tchla)

ggplot(filter(b_model, tchla >0))+
  geom_path(aes(x = tchla, y = fluo_f, colour = as.factor(optical_layer)))+
  geom_line(aes(x = tchla, y = fluo_f_min, colour = as.factor(optical_layer)), linetype = 2)+
  geom_line(aes(x = tchla, y = fluo_f_max, colour = as.factor(optical_layer)), linetype = 2)+
  coord_trans(x = "log")+scale_x_continuous(breaks = c(0.01, 0.10, 1, 10))+
  scale_color_brewer(palette = "Set1")+
  ylab("f/chla ratio")+ xlab("Chlorophylle a")+
  theme_bw()



#rda####
rda_table <- na.omit(select(argo_calibration, micro, nano, pico, ratio, zze, tchla))
rda_bio <- rda(select(rda_table, micro, nano, pico)~rda_table$ratio + rda_table$zze + rda_table$tchla)

summary(rda_bio)

plot(rda_bio)
anova(rda_bio, by = "axis")

summary(lm(ratio~micro+nano+pico, data = argo_calibration))

#pca####

pc <- function(.data, optical = 1){
  t <- filter(.data, optical_layer %in% optical)
  t <- na.omit(select(t, pigments, ratio, micro, nano, pico, zze, tchla))
  t <- t %>% mutate(fuco = fuco/tchla,
                             peri = peri/tchla,
                             hex = hex/tchla,
                             but = but/tchla,
                             allo = allo/tchla,
                             zea = zea/tchla,
                             tchlb = tchlb/tchla)
  pca_t <- PCA(t, quanti.sup = c(9, 10, 11, 12, 13), graph = FALSE)
  descriptors_scores = rbind(as.data.frame(pca_t$var$coord),as.data.frame(pca_t$quanti.sup$coord))
  descriptors_scores$index = c("In","In", "In", "In", "In", "In", "In","Supp", "Supp", "Supp", "Supp", "Supp", "Supp")
  sta_scores = as.data.frame(pca_t$ind$coord)
  print(summary(pca_t))
  pca_t <- bind_cols(t, sta_scores)
  
  g <- ggplot()+
    geom_point(aes(x = Dim.1, y = Dim.2, colour = ratio), data = pca_t)+
    geom_segment(aes(x = 0, y = 0 , xend = Dim.1*4.5, yend = Dim.2*4.5, linetype = index), size = 0.7, arrow = arrow(length=unit(0.2, "cm")), data = descriptors_scores)+
    geom_text(aes(x = Dim.1*6.2, y = Dim.2*4.5, label = rownames(descriptors_scores)), size = 4.5, data = descriptors_scores, nudge_y = -0.1)+
    xlim(-6.5,6.5) + ylim(-4.5,4.5)+
    labs(x = "PC1", y = "PC2", title = paste("optical layer :", optical, sep = ""))+
    scale_color_viridis_c(limits = c(0,8))
  return(g)
}
g1 <- pc(merged_argo, 1)
g1
g2 <- pc(merged_argo, 2)
g2
g3 <- pc(merged_argo, 3)
g3
g4 <- pc(merged_argo, c(1,2,3))
g4

gridExtra::grid.arrange(g1,g2,g3,g4, ncol = 2)


dataset_filtered <- read_csv("Scripts/Data/Final/data_phi.csv")
biosope <- filter(dataset_filtered, id == "biosope")
biosope$ratio <- biosope$fluo/biosope$tchla

g5 <- pc(biosope, 1)
g6 <- pc(biosope, 2)
g7 <- pc(biosope, 3)
g8 <- pc(biosope, c(1,2,3))
gridExtra::grid.arrange(g5,g6,g7,g8, ncol = 2)

g8

ggplot(biosope)+
  geom_point(aes(x = lon, y = -depth, colour = peri))

#redefine optical_layer####

merged_argo <- merged_argo %>% mutate(pd = ze/4.6)
merged_argo <- merged_argo[is.na(merged_argo$pd) == FALSE,]
for(i in 1:nrow(merged_argo)){
  if(merged_argo[i,]$depth <= merged_argo[i,]$pd){
    merged_argo[i,]$optical_layer = 1
  }
  if(merged_argo[i,]$depth <= merged_argo[i,]$ze & merged_argo[i,]$depth > merged_argo[i,]$pd){
    merged_argo[i,]$optical_layer = 2
  }
  if(merged_argo[i,]$depth <= 1.5 * merged_argo[i,]$ze & merged_argo[i,]$depth > merged_argo[i,]$ze){
    merged_argo[i,]$optical_layer = 3
  }
}

