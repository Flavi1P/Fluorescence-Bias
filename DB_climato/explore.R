library(tidyvers)
library(FactoMineR)
library(patchwork)
library(janitor)
library(readxl)
library(wesanderson)

hplc <- read_csv("DB_climato/Data/database_final")
map <- read_csv("Data/map_vec")
spectre <- read_excel("Biosope/Data/Spectres_annick.xlsx")
spectre <- clean_names(spectre)

#select the absorbtion at 440 and 470nm
spectre430<- filter(spectre, lambda == 430)
spectre532<- filter(spectre, lambda == 532)


hplc <- hplc %>% mutate(photo_430 = peri * spectre430$peri + but * spectre430$x19_bf + hex * spectre430$x19_hf + fuco * spectre430$fuco + allo * spectre430$allox + chla * spectre430$chl_a + dv_chla * spectre430$dv_chla,
                              photo_532 = peri * spectre532$peri + but * spectre532$x19_bf + hex * spectre532$x19_hf + fuco * spectre532$fuco + allo * spectre532$allox + chla * spectre532$chl_a + dv_chla * spectre532$dv_chla)



hplc <- mutate(hplc, tchla = chla + dv_chla)

hplc$system <- ifelse(hplc$ze_morel > hplc$mld, "Stratified", "Mixed")

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
  geom_text(aes(x = 2, y = 0.5, label = paste("slope :", round(coeff440, 2), "; R² :", round(r440, 2), sep = " "), colour = "440"), show.legend = FALSE)+
  geom_point(aes(x = tchla, y = photo_470, colour = "470"))+
  geom_smooth(aes(x = tchla, y = photo_470, colour = "470"), se = FALSE, method = "lm", show.legend = FALSE)+
  geom_text(aes(x = 2, y = 0.47, label = paste("slope :", round(coeff430, 2),"; R² :", round(r430, 2), sep = " "), colour = "430"), show.legend = FALSE)+
  geom_point(aes(x = tchla, y = photo_430, colour = "430"))+
  geom_smooth(aes(x = tchla, y = photo_430, colour = "430"), se = FALSE, method = "lm", show.legend = FALSE)+
  geom_text(aes(x = 2, y = 0.44, label = paste("slope :", round(coeff470, 2), "; R² :", round(r470, 2), sep = " "), colour = "470"), show.legend = FALSE)+
  geom_point(aes(x = tchla, y = photo_532, colour = "532"))+
  geom_smooth(aes(x = tchla, y = photo_532, colour = "532"), se = FALSE, method = "lm", show.legend = FALSE)+
  geom_text(aes(x = 2, y = 0.41, label = paste("slope :", round(coeff532, 2), "; R² :", round(r532, 2), sep = " "), colour = "532"), show.legend = FALSE)+
  scale_color_manual(values = wes_palette("Darjeeling2"))+
  ylab("Photosynthetical absorbtion")+
  labs(color = "Wavelength")+
  theme_classic()



ggplot(hplc)+
  geom_point(aes(x = lon, y = lat, colour = system))+
  geom_path(aes(x = long, y = lat, group = group), data = map)+
  coord_quickmap()

ggplot(hplc)+
  geom_point(aes(x = tchla, y = ratio))+
  geom_path(aes(x = tchla, y = 3), colour = "red")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  facet_wrap(.~system)

ggplot(hplc)+
  geom_point(aes(x = tchla, y = photo_440, colour = "aps 440"))+
  geom_point(aes(x = tchla, y = photo_470, colour = "aps 470"))+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  facet_wrap(.~system)

model <- lm(ratio~chla + peri + fuco + dv_chla + allo + hex + but, data = hplc)
summary(model)

pca_model <- PCA(select(hplc, chla, peri, fuco , allo , but , hex, dv_chla, ratio), scale.unit = TRUE)

#transform data to avoid the effect of "when there is a lot of pigment there is a lot of all"

hplc_pca <- hplc %>% select(chla, peri, fuco, allo, but, hex, dv_chla, ratio) %>% 
  mutate(pigsum = rowSums(.),
        chla = chla/pigsum,
        peri = peri/pigsum,
        fuco = fuco/pigsum,
        allo = allo/pigsum,
        but = but/pigsum,
        hex = hex/pigsum,
        dv_chla = dv_chla/pigsum) %>% 
  filter(pigsum > 0) %>% 
  select(-pigsum)

pca_model <- PCA(hplc_pca, scale.unit = TRUE)

hplc <- filter(hplc, chla != 0)
cca_model <- CA(select(hplc, chla, peri, fuco , allo , but , hex, dv_chla, ratio), col.sup = 8)

arrows <- as.data.frame(cca_model$col$coord)

ggplot(arrows)+
  geom_segment(aes(x = 0, xend = arrows$`Dim 1`, y =0, yend = arrows$`Dim 2`))+
  geom_text(aes(x = arrows$`Dim 1`, y = arrows$`Dim 2`, label = rownames(arrows))) 
