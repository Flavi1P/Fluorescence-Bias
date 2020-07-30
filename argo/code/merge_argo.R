library(tidyverse)
library(fuzzyjoin)
library(rgeos)
library(plyr)
library(dplyr)
library(geosphere)
library(readxl)
library(lubridate)
library(nnls)
library(Metrics)
library(gridExtra)
source("functions/profile_numb.R")
source("functions/zeu_moma.R")
source("functions/phi_lm.R")
source("functions/phi_stat.R")
source("functions/phi_boot.R")




pigments <- c("fuco", "peri", "hex", "but", "allo", "tchlb", "zea")

#global lovbio ####

hplc <- read_csv("Data/argo/hplc_argo")
ref <- read_csv("Data/argo/ref.csv")
ref_bis <- read_csv("Data/argo/ref_bis")
map_vec <- read_csv("Data/map_vec")

#first_profiles <- read_csv("Data/argo/first_profiles")
first_profiles <- read_csv("Data/argo/data_argo_filter")

first_profiles <- first_profiles[-1,]
first_profiles$date <- date(first_profiles$date)

descent_profiles <- read_csv("Data/argo/descent_profiles")

ref <- bind_rows(ref, ref_bis)

hplc <- left_join(hplc, ref, by = c("id" = "lovbio"))
first_profiles$number <- as.numeric(first_profiles$id)
first_profiles <- left_join(first_profiles, ref)

first_profiles <- left_join(first_profiles, ref_bis)

descent_profiles$number <- as.numeric(descent_profiles$id)
descent_profiles <- left_join(descent_profiles, ref)


first_profiles$depth <- round(first_profiles$pres)
descent_profiles$depth <- round(descent_profiles$pres)
hplc$depth <- round(hplc$depth)


hplc$lat <- ifelse(hplc$id == "lovbio075b", - hplc$lat, hplc$lat)
hplc$lat <- ifelse(hplc$id %in% c("lovbio077b"), - hplc$lat, hplc$lat)
hplc$lon <- ifelse(hplc$id %in% c("lovbio077b"), - hplc$lon, hplc$lon)
hplc$lat <- ifelse(hplc$id %in% c("lovbio079b"), - hplc$lat, hplc$lat)
hplc$lon <- ifelse(hplc$id %in% c("lovbio079b"), - hplc$lon, hplc$lon)
first_profiles$date[first_profiles$lovbio == "lovbio057b"] <- date("2015-01-20")
first_profiles$date[first_profiles$lovbio == "lovbio043b"] <- date("2015-01-20")


#those two lines correct the date because the hplc sampling was done on this exact position at this day


hplc$profile <- profile_numb(hplc$depth, "upward")



lovbio <- unique(hplc$id)
lovbio <- lovbio[lovbio %in% unique(first_profiles$lovbio)]
lovbio <- lovbio[lovbio %in% unique(descent_profiles$lovbio)]

hplc$lon_round <- round(hplc$lon,1)
hplc$lat_round <- round(hplc$lat,1)
merged_dist <- data.frame(matrix(ncol = 33, nrow = 0))
for (i in lovbio){
  t1 <- filter(first_profiles, lovbio == i)
  #t1 <- filter(descent_profiles, lovbio == i)
  t2 <- filter(hplc, id == i)
  if(i == "lovbio079b"){
    t2 <- filter(t2, date != "2015-03-20")
  }
  if(!is.na(unique(t1$lat))){
  hplc_point <- cbind(t2$lon, t2$lat)
  argo_point <- cbind(rep(unique(t1$lon), length(hplc_point[,1])), rep(unique(t1$lat), length(hplc_point[,1])))
  distance <- distHaversine(hplc_point, argo_point)
  t2_lonround <- t2$lon_round[which(distance == min(distance))]
  t2_latround <- t2$lat_round[which(distance == min(distance))]
  t2 <- filter(t2, lon_round == unique(t2_lonround) & lat_round == unique(t2_latround))
  merged_t <- difference_left_join(t2,t1, by = c("depth"), max_dist = 1) %>%
    mutate(depth_abs_diff = abs(depth.x - pres)) %>% 
    filter(chla != "NA") %>% 
    arrange(depth.x, depth_abs_diff) %>% 
    group_by(depth.x) %>%
    slice(1) %>% 
    ungroup()
  names(merged_dist) <- colnames(merged_t)
  merged_dist <- bind_rows(merged_dist, merged_t)
  }
}




merged_dist$lag <- merged_dist$date.y - merged_dist$new_date
table(merged_dist$lag)

ggplot(filter(merged_dist, chla_qc < 4))+
  geom_point(aes(x = tchla, y = chl_smooth, colour = as.numeric(lag)))

merged_dist_clean <- filter(merged_dist, abs(lag) < 2)

lag_by_float <- data.frame("lovbio" = NA, "lag_min" = NA)
for (i in unique(merged_dist$lovbio)){
  t <- filter(merged_dist, lovbio == i)
  min_lag <- min(t$lag)
  ttable <- data.frame("lovbio" = i, "lag_min" = min_lag)
  lag_by_float <- bind_rows(lag_by_float, ttable)
}

merged_dist_clean <- left_join(merged_dist_clean, lag_by_float) %>% filter(lag == lag_min)

ggplot(merged_dist_clean)+
  geom_point(aes(x = lon.x, y = lat.x, colour = "hplc"), size = 3)+
  geom_point(aes(x = lon.y, y = lat.y, colour = "argo"))+
  geom_polygon(aes(x = long, y = lat, group = group), data = map_vec)+
  coord_quickmap()


merged_dist_clean <- merged_dist_clean[order(merged_dist_clean$id.x),]
merged_dist_clean$profile <- profile_numb(merged_dist_clean$depth.x, "downward")

NAT_LAS_list <- c("lovbio014b", "lovbio030b", "lovbio032b", "lovbio028b", "lovbio011b", "lovbio021c", "lovbio012b")
NAT_LAS <- filter(merged_dist_clean, lovbio %in% NAT_LAS_list)
LAS_profiles <- filter(first_profiles, lovbio %in% NAT_LAS$lovbio)
ggplot(LAS_profiles)+
  geom_path(aes(x = chl_smooth, y = -pres, colour = lovbio))+
  ylim(-250,0)

ggplot(NAT_LAS)+
  geom_point(aes(x = lon.x, y = lat.x, colour = "hplc"), size = 3)+
  geom_point(aes(x = lon.y, y = lat.y, colour = "argo"))+
  geom_polygon(aes(x = long, y = lat, group = group), data = map_vec)+
  coord_quickmap()
ggplot(NAT_LAS)+
  geom_point(aes(x = tchla, y = chla, colour = depth.x))


NAT_ICB_list <-  c("lovbio061c", "lovbio022b", "lovbio029b", "lovbio020b", "lovbio023b", "lovbio013b", "lovbio025c", "lovbio038b")
NAT_ICB <- filter(merged_dist_clean, lovbio %in% NAT_ICB_list)
ICB_profiles <- filter(first_profiles, lovbio %in% NAT_ICB$lovbio)
ggplot(ICB_profiles)+
  geom_path(aes(x = chla_adjusted, y = -pres, colour = lovbio))
ggplot(NAT_ICB)+
  geom_point(aes(x = tchla, y = chla))

NAT_IRS_list <- c("lovbio059c", "lovbio045b", "lovbio024c", "lovbio044b", "lovbio031c", "lovbio027b", "lovbio040b", "lovbio026c")
NAT_IRS <- filter(merged_dist_clean, lovbio %in% NAT_IRS_list)
IRS_profiles <- filter(first_profiles, lovbio %in% NAT_IRS$lovbio)

ggplot(IRS_profiles)+
  geom_path(aes(x = chla_adjusted, y = -pres, colour = lovbio))+
  ylim(-200,0)
ggplot(NAT_IRS)+
  geom_point(aes(x = tchla, y = chla))

ggplot(merged_dist_clean)+
  geom_point(aes(x = tchla, y = chla))+
  geom_point(aes(x = tchla, y = chla), data = NAT_IRS, colour = "Red")

g1 <- ggplot()+
  geom_path(aes(x = chla_adjusted, y = -pres, colour = lovbio), data = filter(ICB_profiles, lovbio != "lovbio061c"), size = 1.1)+
  geom_point(aes(x = tchla, y = -depth.x), data = filter(NAT_ICB, lovbio != "lovbio061c"))+
  ggtitle("Bassin Islandais")+
  ylim(-250,0)+
  xlab("Chlorophylle a")+
  ylab("Profondeur")+
  scale_color_brewer(name = "Flotteur", palette = "Set2")+
  theme_bw(base_size = 20)

g2 <- ggplot()+
  geom_path(aes(x = chla_adjusted, y = -pres, colour = lovbio), data = filter(IRS_profiles, lovbio != "lovbio059c"))+
  geom_point(aes(x = tchla, y = -depth.x), data = filter(NAT_IRS, lovbio != "lovbio059c"))+
  ggtitle("NAT_IRS")+
  ylim(-250,0)+
  theme_bw()

g3 <- ggplot()+
  geom_path(aes(x = chla_adjusted, y = -pres, colour = lovbio), data = LAS_profiles, size = 1.1)+
  geom_point(aes(x = tchla, y = -depth.x), data = NAT_LAS)+
  ggtitle("Mer du Labrador")+
  ylim(-250,0)+
  xlab("Chlorophylle a")+
  ylab("Profondeur")+
  scale_color_brewer(name = "Flotteur", palette = "Set2")+
  theme_bw(base_size = 20)

grid.arrange(g1,g3, ncol = 2)

#png(file = "argo/Plots/profiles_comp.png", width = 12.25, height = 7.73, units = "in", res = 100)
#grid.arrange(g1,g3, ncol = 2)
#dev.off()


#Takuvik data ####

hplc_tak <- read_csv("Data/hplc_tak")
hplc_tak$lon_round <- round(hplc_tak$lon)
hplc_tak$lat_round <- round(hplc_tak$lat)
hplc_tak <- filter(hplc_tak, depth > 5)
position_tak <- cbind(hplc_tak$lon, hplc_tak$lat)
tak_profiles <- first_profiles[grep("takapm", first_profiles$lovbio),]


merged_dist2 <- data.frame(matrix(ncol = 29, nrow = 0))

for (i in unique(tak_profiles$lovbio)){
  t <- filter(tak_profiles, lovbio == i)
  t_position <- cbind(rep(unique(t$lon), length(position_tak[,1])),  rep(unique(t$lat), length(position_tak[,1])))
  distance <- distHaversine(position_tak, t_position)
  tak_lonround <- unique(hplc_tak$lon_round[which(distance == min(distance))])
  tak_latround <- unique(hplc_tak$lat_round[which(distance == min(distance))])
  t2 <- filter(hplc_tak, lon_round == tak_lonround & lat_round == tak_latround)
  merged_t <- left_join(t2,t, by = c("depth"))
  names(merged_dist2) <- colnames(merged_t)
  merged_dist2 <- bind_rows(merged_dist2, merged_t)
}

merged_dist2 <- merged_dist2 %>% group_by(depth, lovbio) %>% summarise_all(mean) %>% ungroup()
merged_dist2 <- filter(merged_dist2, is.na(chla) == FALSE)
merged_dist2 <- filter(merged_dist2, chla_adjusted_qc != 4)

ggplot(merged_dist2)+
  geom_point(aes(x = tchla, y = chl_smooth))+
  ylim(0,2.5)

#Soclim####

hplc_soclim <- read_csv("Data/hplc_soclim")
hplc_soclim$date <- date(hplc_soclim$date)
hplc_soclim$depth <- round(hplc_soclim$depth)

profiles_soclim <- filter(first_profiles, lovbio %in% ref_bis$lovbio)
profiles_soclim$depth <- round(profiles_soclim$pres)

merged_soclim <- left_join(hplc_soclim, profiles_soclim, by = c("date", "depth"))

merged_soclim <- merged_soclim %>% filter(is.na(chla) == FALSE) %>% group_by(lovbio, depth, id.x) %>% summarise_all(mean) %>% ungroup()
merged_soclim <- filter(merged_soclim, lovbio %in% ref_bis$lovbio & lovbio != "lovbio103c" | lovbio == "lovbio103c" & lon.x > 72) 

ggplot(merged_soclim)+
  geom_point(aes(x = chla, y = chla))

#mobydick####

hplc_mduf <- read_csv("Data/argo/hplc_mduf")
hplc_mduf$date <- date(hplc_mduf$date) + 1

merged_mduf <- left_join(hplc_mduf, filter(first_profiles, lovbio == "lovbio111b"), by = c("date", "depth")) %>% filter(is.na(chla) == FALSE)

merged_mduf <- merged_mduf %>% group_by(depth, date, lovbio) %>% summarise_all(mean)

ggplot(merged_mduf)+
  geom_point(aes(x = tchla, y = chla))

#merge all database####


merged_data <- select(merged_dist_clean, date.y, date.x, id.x, depth.x, lon.x, lat.x, lon.y, lat.y, pigments, tchla, chla, chla_qc, chla_adjusted, chla_adjusted_qc, chl_smooth)
merged_data2 <- select(merged_dist2, date.y, date.x, lovbio, depth, lon.x, lat.x, lon.y, lat.y, pigments, tchla, chla, chla_qc, chla_adjusted, chla_adjusted_qc, chl_smooth)

merged_data <- rename_(merged_data, "lovbio" = "id.x", "depth" = "depth.x")
merged_data2 <- rename_(merged_data2)

merged_soclim <- select(merged_soclim, date, lovbio, depth, lon.x, lat.x, lon.y, lat.y, pigments, tchla, chla, chla_qc, chla_adjusted, chla_adjusted_qc, chl_smooth)

merged_mduf <- select(merged_mduf, date, lovbio, depth, lon.x, lat.x, lon.y, lat.y, pigments, tchla, chla, chla_qc, chla_adjusted, chla_adjusted_qc, chl_smooth)

merged_data2$'date.x' <- date(merged_data2$'date.x')
merged_data$'date.x' <- as_date(merged_data$'date.x')

merged_full <- bind_rows(merged_data, merged_data2)
table(merged_full$date.y - merged_full$date.x)
merged_full <- merged_full %>% mutate(lag = abs(date.x - date.y)) %>% 
  filter(lag < 2) %>% 
  select(-date.x) %>%
  rename_('date' = 'date.y')

merged_full <- bind_rows(merged_full, merged_soclim)
merged_full <- bind_rows(merged_full, merged_mduf)




merged_full <- merged_full[order(merged_full$lovbio),]
merged_full$profile <- profile_numb(merged_full$depth, "downward")

ggplot(merged_full)+
  geom_point(aes(x = lon.x, y = lat.x, colour = "hplc"), size = 5)+
  geom_point(aes(x = lon.y, y = lat.y, colour = "argo"), size = 3)+
  geom_polygon(aes(x = long, y = lat, group = group), data = map_vec)+
  coord_quickmap()+
  scale_color_brewer(name = "", palette = "Set1")+
  theme_bw(base_size = 20)+
  xlab("lon")+ ylab("lat")

#ggsave("argo/Plots/merge_map.png", scale = 2)

ggplot(merged_full)+
  geom_point(aes(x = tchla, y = chl_smooth))

table(merged_full$profile)


momadata <- merged_full %>% filter(lovbio != "lovapm002a" & profile != 27) %>% select(depth, tchla, profile, lovbio)
momadata <- na.omit(momadata)


momadata$ze <- NA
for (i in momadata$profile){
  dat <- filter(momadata, profile == i)
  momadata$ze[which(momadata$profile == i)] <- Zeu_moma(dat$tchla, dat$depth)
}

momadata <- momadata %>% ungroup() %>%  select(lovbio, ze) %>% group_by(lovbio) %>%  summarise_all(mean) %>% ungroup()
merged_full <- left_join(merged_full, momadata)

merged_argo <- merged_full %>% mutate(  zze = depth /ze,
                                        optical_layer = floor(zze/0.500001)+1,
                                        wdp = 1.56 * fuco + 0.92 * peri + 4.8 * allo + 1.02 * but + 1.12 * hex + 1.51 * zea + 0.69 * tchlb,
                                        micro = (1.56 * fuco + 0.92 * peri)/wdp,
                                        nano = (4.8 * allo + 1.02 * but + 1.51 * hex)/wdp,
                                        pico = (1.51 * zea + 0.69 * tchlb)/wdp,
                                        amicro = 0.0151 * exp(0.71*(depth/ze)),
                                        anano = 0.0822 * exp(-0.51*(depth/ze)),
                                        apico = 0.1116 * exp(-0.37*(depth/ze)),
                                        microfluo = amicro * micro * tchla,
                                        nanofluo = anano * nano * tchla,
                                        picofluo = apico * pico * tchla,
                                        a = amicro * micro + anano * nano + apico * pico,
                                        zeafluo = a * zea,
                                        fucofluo = a * fuco,
                                        perifluo = a * peri,
                                        allofluo = a * allo,
                                        tchlbfluo = a * tchlb,
                                        hexfluo = a * hex,
                                        butfluo = a * but,
                                        microquanti = micro * tchla,
                                        nanoquanti = nano * tchla,
                                        picoquanti = pico * tchla) %>% ungroup()


table(merged_argo$optical_layer)
merged_argo$optical_layer <- replace_na(merged_argo$optical_layer, 0)
merged_argo <- filter(merged_argo, optical_layer < 4)
test <- filter(merged_argo, lovbio == 'lovbio103c')

ggplot(merged_argo)+
  geom_point(aes(x = tchla, y = chl_smooth, colour = allo))+
  coord_trans(x = 'log', y = 'log')

database <- left_join(first_profiles, merged_argo)

ggplot(database)+
  geom_path(aes(x = chl_smooth, y = -depth))+
  geom_point(aes(x = tchla , y = - depth), colour = 'olivedrab3')+
  geom_point(aes(x = allo, y = -depth, colour = 'allo'))+
  geom_point(aes(x = zea, y = -depth, colour = 'zea'))+
  geom_point(aes(x = fuco, y = -depth, colour = 'fuco'))+
  facet_wrap(.~lovbio, scales = 'free_x')+
  theme_light()+
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )+
  scale_color_brewer(palette = 'Set1', name = '')+
  ylim(-200, 0)

cruise <- tibble('lovbio' = unique(database$lovbio), 'cruise' = c('REY1', 'REY1', 'Labrador', 'Labrador', 'OISO', 'OUTPACE', 'Bioargomed', 'NA', 'Takuvik', 'Bioargomed', 'Bioargomed', 'Bioargomed', 'Bioargomed', 'Bioargomed','Moose', 'Boussole', 'Bioargomed', 'Bioargomed', 'REY1', 'REY1', 'REY1', 'REY1', 'REY1', 'GEOVIDE', 'Islande2013', 'Islande2013', 'Islande2013', 'GEOVIDE', 'Islande2013', 'Islande2013', 'Islande2013', 'Islande2013', 'Labrador', 'Labrador', 'OISO', 'OISO', 'OISO', 'OISO', 'OUTPACE', 'Dewex', 'Moose', 'OUTPACE', 'Takuvik', 'Dewex', 'Takuvik', 'Takuvik', 'Bioargomed', 'NA'))
database <- left_join(database, cruise)

ggplot(database)+
  geom_path(aes(x = chl_smooth, y = -depth, group = lovbio))+
  geom_point(aes(x = tchla , y = - depth), colour = 'olivedrab3')+
  geom_point(aes(x = allo, y = -depth, colour = 'allo'))+
  geom_point(aes(x = zea, y = -depth, colour = 'zea'))+
  geom_point(aes(x = fuco, y = -depth, colour = 'fuco'))+
  facet_wrap(.~cruise, scales = 'free_x')+
  theme_light()+
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )+
  scale_color_brewer(palette = 'Set1', name = '')+
  ylim(-200, 0)

#write_csv(database, 'Data/database_argo')

#write_csv(merged_argo, "Data/merged_argo")



