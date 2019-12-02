library(tidyverse)
library(FactoMineR)
library(patchwork)
library(janitor)
library(readxl)
library(wesanderson)
library(PNWColors)
library(sf)
source("functions/Zone_raph.R")


hplc <- read_csv("DB_climato/Data/lov_climato")
map <- read_csv("Data/map_vec")
spectre <- read_excel("Biosope/Data/Spectres_annick.xlsx")
spectre <- clean_names(spectre)

#select the absorbtion at 440 and 470nm
spectre430<- filter(spectre, lambda == 430)
spectre440 <- filter(spectre, lambda == 440)
spectre470 <- filter(spectre, lambda == 470)
spectre532<- filter(spectre, lambda == 532)


hplc <- hplc %>% mutate(photo_430 = peri * spectre430$peri + but * spectre430$x19_bf + hex * spectre430$x19_hf + fuco * spectre430$fuco + allo * spectre430$allox + chla * spectre430$chl_a + dv_chla * spectre430$dv_chla,
                              photo_532 = peri * spectre532$peri + but * spectre532$x19_bf + hex * spectre532$x19_hf + fuco * spectre532$fuco + allo * spectre532$allox + chla * spectre532$chl_a + dv_chla * spectre532$dv_chla,
                        abs_440 = peri * spectre440$peri + but * spectre440$x19_bf + hex * spectre440$x19_hf + fuco * spectre440$fuco + allo * spectre440$allox + chla * spectre440$chl_a + dv_chla * spectre440$dv_chla + spectre440$zea * zea + spectre440$chl_b * chlb + spectre440$dv_chlb * dv_chlb + spectre440$chlc12 * chlc1c2 + spectre440$a_car * a_caro + spectre440$diad * diad + spectre440$ss_car * b_caro,
                        abs_430 = peri * spectre430$peri + but * spectre430$x19_bf + hex * spectre430$x19_hf + fuco * spectre430$fuco + allo * spectre430$allox + chla * spectre430$chl_a + dv_chla * spectre430$dv_chla + spectre430$zea * zea + spectre430$chl_b * chlb + spectre430$dv_chlb * dv_chlb + spectre430$chlc12 * chlc1c2 + spectre430$a_car * a_caro + spectre430$diad * diad + spectre430$ss_car * b_caro,
                        abs_470 = peri * spectre470$peri + but * spectre470$x19_bf + hex * spectre470$x19_hf + fuco * spectre470$fuco + allo * spectre470$allox + chla * spectre470$chl_a + dv_chla * spectre470$dv_chla + spectre470$zea * zea + spectre470$chl_b * chlb + spectre470$dv_chlb * dv_chlb + spectre470$chlc12 * chlc1c2 + spectre470$a_car * a_caro + spectre470$diad * diad + spectre470$ss_car * b_caro,
                        abs_532 = peri * spectre532$peri + but * spectre532$x19_bf + hex * spectre532$x19_hf + fuco * spectre532$fuco + allo * spectre532$allox + chla * spectre532$chl_a + dv_chla * spectre532$dv_chla + spectre532$zea * zea + spectre532$chl_b * chlb + spectre532$dv_chlb * dv_chlb + spectre532$chlc12 * chlc1c2 + spectre532$a_car * a_caro + spectre532$diad * diad + spectre532$ss_car * b_caro,
                        protect_440 =  spectre440$zea * zea + spectre440$a_car * a_caro + spectre440$diad * diad + spectre440$ss_car * b_caro,
                        protect_470 =  spectre470$zea * zea + spectre470$a_car * a_caro + spectre470$diad * diad + spectre470$ss_car * b_caro,
                        protect_532 =  spectre532$zea * zea + spectre532$a_car * a_caro + spectre532$diad * diad + spectre532$ss_car * b_caro,
                        a_ps_440 = abs_440 - protect_440,
                        a_ps_470 = abs_470 - protect_470,
                        a_ps_532 = abs_532 - protect_532,
                        ratio = a_ps_440/a_ps_470,
                        ratio440_532 = a_ps_440/a_ps_532,
                        ratio470_532 = a_ps_470/a_ps_532
                        )



hplc <- mutate(hplc, tchla = chla + dv_chla)

hplc$system <- ifelse(hplc$ze_morel > hplc$mld, "Stratified", "Mixed")


# add the bioregion

path = "Data/Longhurst"
longhurst_sf <- read_sf(dsn = path.expand(path), quiet = TRUE)

names(longhurst_sf) <- c("code", "region", "geometry")

pnts_sf <- do.call("st_sfc",c(lapply(1:nrow(hplc),
                                     function(i) {st_point(as.numeric(hplc[i,c("lon", "lat") ]))}), list("crs" = 4326))) 
pnts_trans <- st_transform(pnts_sf, 4326)
longhurst_trans <- st_transform(longhurst_sf, 4326)  
hplc$code <- apply(st_intersects(longhurst_trans, pnts_trans, sparse = FALSE), 2, 
                    function(col) { 
                      longhurst_trans[which(col), ]$code
                    })


coeff430 <- summary(lm(photo_430~tchla, data = hplc))$coefficient[2,1]
coeff440 <- summary(lm(photo_440~tchla, data = hplc))$coefficient[2,1]
coeff470 <- summary(lm(photo_470~tchla, data = hplc))$coefficient[2,1]
coeff532 <- summary(lm(photo_532~tchla, data = hplc))$coefficient[2,1]

r430 <- summary(lm(photo_430~tchla, data = hplc))$adj.r.squared
r440 <- summary(lm(photo_440~tchla, data = hplc))$adj.r.squared
r470 <- summary(lm(photo_470~tchla, data = hplc))$adj.r.squared
r532 <- summary(lm(photo_532~tchla, data = hplc))$adj.r.squared






ggplot(hplc)+
  geom_point(aes(x = tchla, y = photo_440, colour = "440"))+
  geom_smooth(aes(x = tchla, y = photo_440, colour = "440"), se = FALSE, method = "lm", show.legend = FALSE)+
  geom_text(aes(x = 2, y = 0.5, label = paste("slope :", round(coeff440, 3), "; R² :", round(r440, 2), sep = " "), colour = "440"), show.legend = FALSE)+
  geom_point(aes(x = tchla, y = photo_470, colour = "470"))+
  geom_smooth(aes(x = tchla, y = photo_470, colour = "470"), se = FALSE, method = "lm", show.legend = FALSE)+
  geom_text(aes(x = 2, y = 0.47, label = paste("slope :", round(coeff430, 3),"; R² :", round(r430, 2), sep = " "), colour = "430"), show.legend = FALSE)+
  geom_point(aes(x = tchla, y = photo_430, colour = "430"))+
  geom_smooth(aes(x = tchla, y = photo_430, colour = "430"), se = FALSE, method = "lm", show.legend = FALSE)+
  geom_text(aes(x = 2, y = 0.44, label = paste("slope :", round(coeff470, 3), "; R² :", round(r470, 2), sep = " "), colour = "470"), show.legend = FALSE)+
  geom_point(aes(x = tchla, y = photo_532, colour = "532"))+
  geom_smooth(aes(x = tchla, y = photo_532, colour = "532"), se = FALSE, method = "lm", show.legend = FALSE)+
  geom_text(aes(x = 2, y = 0.41, label = paste("slope :", round(coeff532, 3), "; R² :", round(r532, 2), sep = " "), colour = "532"), show.legend = FALSE)+
  scale_color_manual(values = wes_palette("Darjeeling2"))+
  ylab("Photosynthetical absorbtion")+
  labs(color = "Wavelength")+
  theme_classic()


coeff430 <- summary(lm(abs_430~tchla, data = hplc))$coefficient[2,1]
coeff440 <- summary(lm(abs_440~tchla, data = hplc))$coefficient[2,1]
coeff470 <- summary(lm(abs_470~tchla, data = hplc))$coefficient[2,1]
coeff532 <- summary(lm(abs_532~tchla, data = hplc))$coefficient[2,1]

r430 <- summary(lm(abs_430~tchla, data = hplc))$adj.r.squared
r440 <- summary(lm(abs_440~tchla, data = hplc))$adj.r.squared
r470 <- summary(lm(abs_470~tchla, data = hplc))$adj.r.squared
r532 <- summary(lm(abs_532~tchla, data = hplc))$adj.r.squared

ggplot(hplc)+
  geom_point(aes(x = tchla, y = abs_440, colour = "440"))+
  geom_smooth(aes(x = tchla, y = abs_440, colour = "440"), se = FALSE, method = "lm", show.legend = FALSE)+
  geom_text(aes(x = 2, y = 0.5, label = paste("slope :", round(coeff440, 3), "; R² :", round(r440, 2), sep = " "), colour = "440"), show.legend = FALSE)+
  geom_point(aes(x = tchla, y = abs_470, colour = "470"))+
  geom_smooth(aes(x = tchla, y = abs_470, colour = "470"), se = FALSE, method = "lm", show.legend = FALSE)+
  geom_text(aes(x = 2, y = 0.47, label = paste("slope :", round(coeff430, 3),"; R² :", round(r430, 2), sep = " "), colour = "430"), show.legend = FALSE)+
  geom_point(aes(x = tchla, y = abs_430, colour = "430"))+
  geom_smooth(aes(x = tchla, y = abs_430, colour = "430"), se = FALSE, method = "lm", show.legend = FALSE)+
  geom_text(aes(x = 2, y = 0.44, label = paste("slope :", round(coeff470, 3), "; R² :", round(r470, 2), sep = " "), colour = "470"), show.legend = FALSE)+
  geom_point(aes(x = tchla, y = abs_532, colour = "532"))+
  geom_smooth(aes(x = tchla, y = abs_532, colour = "532"), se = FALSE, method = "lm", show.legend = FALSE)+
  geom_text(aes(x = 2, y = 0.41, label = paste("slope :", round(coeff532, 3), "; R² :", round(r532, 2), sep = " "), colour = "532"), show.legend = FALSE)+
  scale_color_manual(values = wes_palette("Darjeeling2"))+
  ylab("Photosynthetical absorbtion")+
  labs(color = "Wavelength")+
  theme_classic()


ggplot(hplc)+
  geom_point(aes(x = lon, y = lat, colour = project))+
  geom_polygon(aes(x = long, y = lat, group = group), data = map)+
  coord_map(projection = "gilbert")+
  theme_bw()

ggplot(filter(hplc, lat < 10 & lat > -10))+
  geom_point(aes(x = chla, y = -depth))+
  facet_wrap(.~ nprof)+
  theme_bw()

hplc[is.na(hplc)] <- 0

choice <- readline(prompt = "Mixed or Stratified ?") # possibility to compute plots for different system, ask user

abs_coef <- data.frame("wavelength" = spectre$lambda, "coef" = NA, "coef_tot" = NA, "coef_protect" = NA, "coef_tchla" = NA, "a_peri" = NA, "a_but" = NA, "a_hex" = NA, "a_fuco" = NA, "a_allo" = NA, "a_chla"= NA, "a_dvchla" = NA, "a_zea" = NA,
                       "a_chlb" = NA, "a_dvchlb" = NA, "a_chlc1c2" = NA, "a_acar" = NA, "a_diad" = NA, "a_bcar" = NA) #create empty df

if(choice == "Mixed" | choice == "Stratified"){
  for(i in spectre$lambda){
    t_spectre <- filter(spectre, lambda == i)
    t_spectre[is.na(t_spectre)] <- 0
    t_hplc <- hplc %>% filter(system == choice) %>%
      mutate(photo_t = peri * t_spectre$peri + but * t_spectre$x19_bf + hex * t_spectre$x19_hf + fuco * t_spectre$fuco + allo * t_spectre$allox + chla * t_spectre$chl_a + dv_chla * t_spectre$dv_chla,
             photo_tot = peri * t_spectre$peri + but * t_spectre$x19_bf + hex * t_spectre$x19_hf + fuco * t_spectre$fuco + allo * t_spectre$allox + chla * t_spectre$chl_a + dv_chla * t_spectre$dv_chla + t_spectre$zea * zea + t_spectre$chl_b * chlb + t_spectre$dv_chlb * dv_chlb + t_spectre$chlc12 * chlc1c2 + t_spectre$a_car * a_caro + t_spectre$diad * diad + t_spectre$ss_car * b_caro,
             protect =  t_spectre$zea * zea + t_spectre$a_car * a_caro + t_spectre$diad * diad + t_spectre$ss_car * b_caro,
             photo_chla = chla * t_spectre$chl_a + dv_chla * t_spectre$dv_chla,
             #compute the absorbtion of each pigment
             a_peri = peri * t_spectre$peri,
             a_but = but * t_spectre$x19_bf,
             a_hex = hex * t_spectre$x19_hf,
             a_fuco = fuco * t_spectre$fuco,
             a_allo = allo * t_spectre$allox,
             a_chla = chla * t_spectre$chl_a,
             a_dvchla = dv_chla * t_spectre$dv_chla, 
             a_zea = zea * t_spectre$zea,
             a_chlb = chlb * t_spectre$chl_b,
             a_dvchlb = dv_chlb * t_spectre$dv_chlb,
             a_chlc1c2 = chlc1c2 * t_spectre$chlc12,
             a_acar = a_caro * t_spectre$a_car,
             a_diad = diad * t_spectre$diad,
             a_bcar = b_caro * t_spectre$ss_car
      )
    a_peri = mean(t_hplc$a_peri)
    a_but = mean(t_hplc$a_but)
    a_hex = mean(t_hplc$a_hex)
    a_fuco = mean(t_hplc$a_fuco)
    a_allo = mean(t_hplc$a_allo)
    a_chla = mean(t_hplc$a_chla)
    a_dvchla = mean(t_hplc$a_dvchla)
    a_zea = mean(t_hplc$a_zea)
    a_chlb = mean(t_hplc$a_chlb)
    a_dvchlb = mean(t_hplc$a_dvchlb)
    a_chlc1c2 = mean(t_hplc$a_chlc1c2)
    a_acar = mean(t_hplc$a_acar)
    a_diad = mean(t_hplc$a_diad)
    a_bcar = mean(t_hplc$a_bcar)
    coef <-  summary(lm(photo_t~tchla, data = t_hplc))$coefficient[2,1]     
    coef_tot <-  summary(lm(photo_tot~tchla, data = t_hplc))$coefficient[2,1]
    coef_protect <-  summary(lm(protect~tchla, data = t_hplc))$coefficient[2,1]
    coef_chla <-  summary(lm(photo_chla~tchla, data = t_hplc))$coefficient[2,1]
    multivariate <- summary(lm(tchla~a_peri + a_but + a_hex + a_fuco + a_allo + a_chla + a_dvchla + a_zea + a_chlb + a_dvchlb + a_chlc1c2 + a_acar + a_diad + a_bcar, data = t_hplc))
    abs_coef[abs_coef$wavelength == i,]$coef <- coef
    abs_coef[abs_coef$wavelength == i,]$coef_tot <- coef_tot
    abs_coef[abs_coef$wavelength == i,]$coef_protect <- coef_protect
    abs_coef[abs_coef$wavelength == i,]$coef_tchla<- coef_chla
    abs_coef[abs_coef$wavelength == i,]$a_peri <- a_peri
    abs_coef[abs_coef$wavelength == i,]$a_but <- a_but
    abs_coef[abs_coef$wavelength == i,]$a_hex <- a_hex
    abs_coef[abs_coef$wavelength == i,]$a_fuco <- a_fuco
    abs_coef[abs_coef$wavelength == i,]$a_allo <- a_allo
    abs_coef[abs_coef$wavelength == i,]$a_chla <- a_chla
    abs_coef[abs_coef$wavelength == i,]$a_dvchla <- a_dvchla
    abs_coef[abs_coef$wavelength == i,]$a_zea <- a_zea
    abs_coef[abs_coef$wavelength == i,]$a_chlb <- a_chlb
    abs_coef[abs_coef$wavelength == i,]$a_dvchlb <- a_dvchlb
    abs_coef[abs_coef$wavelength == i,]$a_chlc1c2 <- a_chlc1c2
    abs_coef[abs_coef$wavelength == i,]$a_acar <- a_acar
    abs_coef[abs_coef$wavelength == i,]$a_diad <- a_diad
    abs_coef[abs_coef$wavelength == i,]$a_bcar <- a_bcar
  }} else{
    for(i in spectre$lambda){
      t_spectre <- filter(spectre, lambda == i)
      t_spectre[is.na(t_spectre)] <- 0
      t_hplc <- hplc %>%
        mutate(photo_t = peri * t_spectre$peri + but * t_spectre$x19_bf + hex * t_spectre$x19_hf + fuco * t_spectre$fuco + allo * t_spectre$allox + chla * t_spectre$chl_a + dv_chla * t_spectre$dv_chla,
               photo_tot = peri * t_spectre$peri + but * t_spectre$x19_bf + hex * t_spectre$x19_hf + fuco * t_spectre$fuco + allo * t_spectre$allox + chla * t_spectre$chl_a + dv_chla * t_spectre$dv_chla + t_spectre$zea * zea + t_spectre$chl_b * chlb + t_spectre$dv_chlb * dv_chlb + t_spectre$chlc12 * chlc1c2 + t_spectre$a_car * a_caro + t_spectre$diad * diad + t_spectre$ss_car * b_caro,
               protect =  t_spectre$zea * zea + t_spectre$a_car * a_caro + t_spectre$diad * diad + t_spectre$ss_car * b_caro,
               photo_chla = chla * t_spectre$chl_a + dv_chla * t_spectre$dv_chla,
               #compute the absorbtion of each pigment
               a_peri = peri * t_spectre$peri,
               a_but = but * t_spectre$x19_bf,
               a_hex = hex * t_spectre$x19_hf,
               a_fuco = fuco * t_spectre$fuco,
               a_allo = allo * t_spectre$allox,
               a_chla = chla * t_spectre$chl_a,
               a_dvchla = dv_chla * t_spectre$dv_chla, 
               a_zea = zea * t_spectre$zea,
               a_chlb = chlb * t_spectre$chl_b,
               a_dvchlb = dv_chlb * t_spectre$dv_chlb,
               a_chlc1c2 = chlc1c2 * t_spectre$chlc12,
               a_acar = a_caro * t_spectre$a_car,
               a_diad = diad * t_spectre$diad,
               a_bcar = b_caro * t_spectre$ss_car
        )
      a_peri = mean(t_hplc$a_peri)
      a_but = mean(t_hplc$a_but)
      a_hex = mean(t_hplc$a_hex)
      a_fuco = mean(t_hplc$a_fuco)
      a_allo = mean(t_hplc$a_allo)
      a_chla = mean(t_hplc$a_chla)
      a_dvchla = mean(t_hplc$a_dvchla)
      a_zea = mean(t_hplc$a_zea)
      a_chlb = mean(t_hplc$a_chlb)
      a_dvchlb = mean(t_hplc$a_dvchlb)
      a_chlc1c2 = mean(t_hplc$a_chlc1c2)
      a_acar = mean(t_hplc$a_acar)
      a_diad = mean(t_hplc$a_diad)
      a_bcar = mean(t_hplc$a_bcar)
      coef <-  summary(lm(photo_t~tchla, data = t_hplc))$coefficient[2,1]     
      coef_tot <-  summary(lm(photo_tot~tchla, data = t_hplc))$coefficient[2,1]
      coef_protect <-  summary(lm(protect~tchla, data = t_hplc))$coefficient[2,1]
      coef_chla <-  summary(lm(photo_chla~tchla, data = t_hplc))$coefficient[2,1]
      multivariate <- summary(lm(tchla~a_peri + a_but + a_hex + a_fuco + a_allo + a_chla + a_dvchla + a_zea + a_chlb + a_dvchlb + a_chlc1c2 + a_acar + a_diad + a_bcar, data = t_hplc))
      abs_coef[abs_coef$wavelength == i,]$coef <- coef
      abs_coef[abs_coef$wavelength == i,]$coef_tot <- coef_tot
      abs_coef[abs_coef$wavelength == i,]$coef_protect <- coef_protect
      abs_coef[abs_coef$wavelength == i,]$coef_tchla<- coef_chla
      abs_coef[abs_coef$wavelength == i,]$a_peri <- a_peri
      abs_coef[abs_coef$wavelength == i,]$a_but <- a_but
      abs_coef[abs_coef$wavelength == i,]$a_hex <- a_hex
      abs_coef[abs_coef$wavelength == i,]$a_fuco <- a_fuco
      abs_coef[abs_coef$wavelength == i,]$a_allo <- a_allo
      abs_coef[abs_coef$wavelength == i,]$a_chla <- a_chla
      abs_coef[abs_coef$wavelength == i,]$a_dvchla <- a_dvchla
      abs_coef[abs_coef$wavelength == i,]$a_zea <- a_zea
      abs_coef[abs_coef$wavelength == i,]$a_chlb <- a_chlb
      abs_coef[abs_coef$wavelength == i,]$a_dvchlb <- a_dvchlb
      abs_coef[abs_coef$wavelength == i,]$a_chlc1c2 <- a_chlc1c2
      abs_coef[abs_coef$wavelength == i,]$a_acar <- a_acar
      abs_coef[abs_coef$wavelength == i,]$a_diad <- a_diad
      abs_coef[abs_coef$wavelength == i,]$a_bcar <- a_bcar
    }
  }

ggplot(abs_coef)+
  geom_path(aes(x = wavelength, y = coef))+
  geom_path(aes(x = wavelength, y = coef_tot), colour = "Grey")+
  geom_path(aes(x = wavelength, y = coef_protect), colour = "Brown")+
  geom_path(aes(x = wavelength, y = coef_tchla), colour = "#addd8e")+
  geom_path(aes(x = wavelength, y = coef_tot - coef_protect), colour = "#43a2ca")+
  ylab("da/dlbd")+
  geom_label(aes(x = 600, y = 0.04, label = choice))+
  scale_x_continuous(breaks = c(400, 440, 470, 500, 532, 600 ,700))+
  theme_bw()

abs_coef$a_sum <- rowSums(abs_coef[,6:19])

abs_long <- pivot_longer(abs_coef, 6:19,  names_to = "absorbtion", names_prefix = "a_")

abs_long$value[abs_long$value < 0] <- 0
ggplot(abs_long)+
  geom_col(aes(x = wavelength, y = value, fill = absorbtion), position = "stack")+
  scale_x_continuous(breaks = c(400, 440, 470, 500, 532, 600, 700))+
  scale_fill_manual(values = pnw_palette("Bay", n = 14))+
  theme_bw()

ggplot(abs_long)+
  geom_area(aes(x = wavelength, y = value, fill = absorbtion))+
  geom_label(aes(x = 600, y = 0.012, label = choice))+
  scale_x_continuous(breaks = c(400, 440, 470, 500, 532, 600, 700))+
  scale_fill_manual(values = pnw_palette("Bay", n = 14))+
  theme_bw()

hplc <- mutate(hplc, zze = depth/ze_morel)
hplc$surf_chla <- NA
for(i in unique(hplc$nprof)){
  t <- filter(hplc, nprof == i)
  surf_chla <- t$tchla[t$depth == min(t$depth)]
  if(surf_chla < 0.05){
    troph <- "<0.05"
  }
  else{
    if(surf_chla < 0.1){
      troph <- "<0.1"
    }
    else{
      if(surf_chla < 0.15){
        troph <- "<0.15"
      }
      else{
        if(surf_chla < 0.2){
          troph <- "<0.2"
        }
        else{
          if(surf_chla < 0.3){
            troph <- "<0.3"
          }
          else{
            if(surf_chla < 0.5){
              troph <- "< 0.5"
            }
            else{
              if(surf_chla < 1){
                troph <- "<1"
              }
              else{
                troph <- "> 1"
              }
            }
          }
        }
      }
    }
  }
  hplc$surf_chla[hplc$nprof == i] <- troph
}

for(i in unique(hplc$nprof)){
  t <- filter(hplc, nprof == i)
  hplc$zone[hplc$nprof == i] <- Zone(t$lat, t$lon)
}

ggplot(filter(hplc, zze <2 & code != "MEDI"))+
  geom_point(aes(x = ratio, y = -zze, colour = code), show.legend = FALSE)+
  facet_wrap(.~surf_chla) +
  xlim(0,4)+
  theme_bw()

ggplot(filter(hplc, zze < 1.5))+
  geom_point(aes(x = chla, y = ratio, colour = zze))+
  scale_color_viridis_c()+
  xlim(0,1)+
  ylim(0,10)+
  theme_bw()

hplc$layer <- floor(hplc$zze)+1

ggplot(hplc)+
  geom_point(aes(x = lon, y = lat, colour = zone))+
  geom_polygon(aes(x = long, y = lat, group = group), data = map)+
  coord_map(projection = "gilbert")+
  theme_bw()

ggplot(filter(hplc, zze <2))+
  geom_point(aes(x = ratio440_532, y = -zze, colour = "440/532"))+
  geom_point(aes(x = ratio470_532, y = -zze, colour = "470/532"))+
  geom_point(aes(x = ratio, y = -zze, colour = "440/570"))+
  facet_wrap(.~zone, scales = "free_x") +
  theme_bw()


ggplot(filter(hplc, project == "EPOPE"))+
  geom_point(aes(x = ratio440_532, y = -zze, colour = "440/532", group = nprof))+
  geom_point(aes(x = ratio470_532, y = -zze, colour = "470/532", group = nprof))+
  geom_point(aes(x = ratio, y = -zze, colour = "440/570", group = nprof))+
  geom_point(aes(x = tchla*4, y = -zze, colour = "Chla", group = nprof))+
  scale_x_continuous(sec.axis = sec_axis(~. /4 , name = "Chla"))+
  xlab("Absorbtion")+
  scale_color_brewer(palette = "Paired")+
  theme_bw()

ggplot(filter(hplc, project == "EPOPE"))+
  geom_path(aes(x = tchla, y = -zze, colour = surf_chla, group = nprof), size = 2)+
  xlab("Absorbtion")+
  scale_color_brewer(palette = "Paired")+
  theme_bw()

ggplot(filter(hplc, project == "EPOPE"))+
  geom_point(aes(x = lon, y = lat, colour = nprof))+
  geom_polygon(aes(x = long, y = lat, group = group), data = map)+
  xlim(-160,-140)+
  ylim(-19,2)+
  coord_quickmap()+
  scale_color_viridis_c()+
  theme_bw()

ggplot(filter(hplc, project == "EPOPE"))+
  geom_path(aes(x = ratio440_532, y = -zze, colour = "440/532", group = nprof))+
  geom_path(aes(x = ratio470_532, y = -zze, colour = "470/532", group = nprof))+
  geom_path(aes(x = ratio, y = -zze, colour = "440/570", group = nprof))+
  geom_path(aes(x = tchla*10, y = -zze, colour = "Chla", group = nprof))+
  geom_path(aes(x = chlb/chla, y = -zze, colour = "chlb/chla", group = nprof))+
  scale_x_continuous(sec.axis = sec_axis(~. /10 , name = "Chla"))+
  xlab("Absorbtion")+
  scale_color_brewer(palette = "Paired")+
  theme_bw()+
  facet_wrap(~surf_chla, scales = "free_x")

ggplot(filter(hplc, project == "EPOPE"))+
  geom_path(aes(x = ratio470_532, y = -zze, colour = "ratio 470/532", group = nprof))+
  geom_path(aes(x = chlb/chla * 5, y = -zze, colour = "chlb/chla", group = nprof))+
  scale_color_viridis_d()+
  scale_x_continuous(sec.axis = sec_axis(~. /5 , name = "Chlb/Chla"))+
  xlab("Absorbtion")+
  facet_wrap(~surf_chla, scales = "free_x")+
  theme_bw()
  
  hplc_resumed <- mutate(hplc, zze = depth/ze_morel,
               zone = ifelse(zze < 0.5, "surface", "depth")) %>% 
  group_by(nprof, zone, system) %>% 
  summarise_if(is.numeric, mean, na.rm = TRUE) %>% 
  ungroup()

system_info <- hplc_resumed %>% group_by(system, zone) %>% 
  summarise_if(is.numeric, c(mean, sd)) %>% 
  ungroup()

ggplot(system_info)+
  geom_col(aes(x = zone, y = tchla_fn1, fill = system), position = "dodge")+
  geom_errorbar(aes(x = zone, ymin = tchla_fn1 - tchla_fn2, ymax = tchla_fn1 + tchla_fn2, fill = system), position = "dodge")+
  theme_classic()+
  scale_fill_manual(values = wes_palette("Darjeeling2"))


