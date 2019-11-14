library(tidyverse)
library(sf)
library(grid)

argo <- read_csv("DB_climato/Data/argo_fluo_corrected")
path = "Data/Longhurst"
longhurst_sf <- read_sf(dsn = path.expand(path), quiet = TRUE)

names(longhurst_sf) <- c("code", "region", "geometry")

#longhurst_sf %>% ggplot() + geom_sf(aes(fill = code))

pnts_sf <- do.call("st_sfc",c(lapply(1:nrow(argo),
                                     function(i) {st_point(as.numeric(argo[i,c("lon.y", "lat.y") ]))}), list("crs" = 4326))) 
pnts_trans <- st_transform(pnts_sf, 4326)
longhurst_trans <- st_transform(longhurst_sf, 4326)  
argo$code <- apply(st_intersects(longhurst_trans, pnts_trans, sparse = FALSE), 2, 
                   function(col) { 
                     longhurst_trans[which(col), ]$code
                   })

#lets mean the 7 profile of sarc because they all correspond to the same launch

sarc <- argo %>% filter(code == "SARC") %>% group_by(depth) %>% summarise_all(mean) %>% ungroup()
sarc$code <- "SARC"
sarc$lovbio <- "lovbio_sarc"

argo <- argo %>% filter(code != "SARC")
argo <- bind_rows(argo, sarc)

region_argo <- argo  %>% mutate(ratio_clean = corrected_fluo/tchla, ratio = fluo/tchla) %>%  group_by(code) %>% summarise_at(vars(c(ratio, ratio_clean, fitted_math)), c(mean, sd), na.rm = TRUE) %>% ungroup()
names(region_argo) <- c("code", "mean_ratio", "mean_ratio_clean", "mean_factor", "sd_ratio", "sd_ratio_clean", "sd_factor")

region_argo <- left_join(region_argo, resume_region)

region_argo$sd <- ifelse(region_argo$sd > region_argo$mean, region_argo$mean, region_argo$sd)


ggplot(region_argo)+
  geom_col(aes(x = reorder(code, mean_ratio), y = mean_ratio, fill = code))+
  geom_errorbar(aes(x = code, ymin = mean_ratio - sd_ratio, ymax = mean_ratio + sd_ratio))+
  xlab("Province océanique")+
  ylab("Rapport Fluo/Chla")+
  geom_errorbar(aes(code, ymax = 1, ymin = 1),
                size=1, linetype = "longdash", inherit.aes = F, width = 1)+
  guides(fill = FALSE)+
  scale_fill_brewer(palette = "Dark2")+
  theme_bw(base_size = 20)+
  ylim(0,9)

ggplot(region_argo)+
  geom_col(aes(x = reorder(code, mean_ratio), y = mean_ratio_clean, fill = code))+
  geom_errorbar(aes(x = code, ymin = mean_ratio_clean - sd_ratio_clean, ymax = mean_ratio_clean + sd_ratio_clean))+
  xlab("Province océanique")+
  ylab("Rapport Fluo/Chla")+
  geom_errorbar(aes(code, ymax = 1, ymin = 1),
                size=1, linetype = "longdash", inherit.aes = F, width = 1)+
  guides(fill = FALSE)+
  scale_fill_brewer(palette = "Dark2")+
  theme_bw(base_size = 20)+
  ylim(0,9)

ggplot(region_argo)+
  geom_col(aes(x = reorder(code, mean_ratio), y = mean_factor, fill = code))+
  geom_errorbar(aes(x = code, ymin = mean_factor - sd_factor, ymax = mean_factor + sd_factor))+
  xlab("Province océanique")+
  ylab("Rapport Fluo/Chla")+
  geom_errorbar(aes(code, ymax = 1, ymin = 1),
                size=1, linetype = "longdash", inherit.aes = F, width = 1)+
  guides(fill = FALSE)+
  scale_fill_brewer(palette = "Dark2")+
  theme_bw(base_size = 20)+
  ylim(0,9)
