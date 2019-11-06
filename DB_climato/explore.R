library(tidyvers)
library(FactoMineR)

hplc <- read_csv("DB_climato/Data/database_final")

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
