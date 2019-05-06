library(tidyverse)
phi_argo <- read_csv("argo/Data/phi_argo")

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

d_model <- d_model %>% mutate(fluo_min = phi_argo$phi)

d_model <- bind_rows(d_model, d_model, d_model)

d_model$optical_layer <- c(rep(1, 201), rep(2,201), rep(3,201))
d_model <- na.omit(d_model)

phi_argo <- phi_argo %>% spread(key = size, value = c("phi", "se"))

se_argo <- phi_argo %>% select(se, optical_layer, size) %>% spread(key = size, value = se)
names(se_argo) <- c("optical_layer", "se_micro", "se_nano", "se_pico") 
phi_argo <- phi_argo %>% select(phi, optical_layer, size) %>% spread(key = size, value = phi)
names(phi_argo) <- c("optical_layer", "phi_micro", "phi_nano", "phi_pico") 

d_model <- left_join(d_model, phi_argo)
d_model <- left_join(d_model, se_argo)


d_model <- d_model %>% mutate(fluo_model = micro_model * phi_micro + nano_model * phi_nano + pico_model * phi_pico,
                              fluo_min = tchla * (micro_f * (phi_micro - se_micro) + nano_f * (phi_nano - se_nano) + pico_f * (phi_pico - se_pico)),
                              fluo_max = tchla * (micro_f * (phi_micro + se_micro) + nano_f * (phi_nano + se_nano) + pico_f * (phi_pico + se_pico)),
                              fluo_f = fluo_model/tchla,
                              fluo_f_min = fluo_min/tchla,
                              fluo_f_max = fluo_max/tchla)

ggplot(filter(d_model, tchla >0))+
  geom_path(aes(x = tchla, y = fluo_f, colour = as.factor(optical_layer)), size = 1.2)+
  geom_line(aes(x = tchla, y = fluo_f_min, colour = as.factor(optical_layer)), linetype = 2)+
  geom_line(aes(x = tchla, y = fluo_f_max, colour = as.factor(optical_layer)), linetype = 2)+
  coord_trans(x = "log")+scale_x_continuous(breaks = c(0.01, 0.10, 1, 10))+
  ylab("fluo/chla ratio")+ xlab("Chlorophylle a")+
  scale_color_viridis_d(name = "Couche optique")+
  theme_bw(base_size = 20)
ggsave("argo/Plots/modele_fluo.png")
