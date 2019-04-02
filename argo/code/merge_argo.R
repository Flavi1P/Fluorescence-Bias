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
source("Scripts/profile_numb.r")
source("Scripts/zeu_moma.r")
source("Scripts/phi_lm.R")
source("Scripts/phi_stat.R")
source("Scripts/phi_boot.R")




pigments <- c("fuco", "peri", "hex", "but", "allo", "tchlb", "zea")

#global lovbio ####

hplc <- read_csv("Scripts/Data/argo/hplc_argo")
ref <- read_csv("Scripts/Data/argo/ref.csv")
ref_bis <- read_csv("Scripts/Data/argo/ref_bis")
map_vec <- read_csv("Scripts/Data/map_vec")
first_profiles <- read_csv("Scripts/Data/argo/first_profiles")
first_profiles <- first_profiles[-1,]
first_profiles$date <- date(first_profiles$date)
ref <- bind_rows(ref, ref_bis)

hplc <- left_join(hplc, ref, by = c("id" = "lovbio"))
first_profiles$number <- as.numeric(first_profiles$id)
first_profiles <- left_join(first_profiles, ref)
first_profiles <- left_join(first_profiles, ref_bis)


first_profiles$depth <- round(first_profiles$pres)
hplc$depth <- round(hplc$depth)


hplc$lat <- ifelse(hplc$id == "lovbio075b", - hplc$lat, hplc$lat)
hplc$lat <- ifelse(hplc$id %in% c("lovbio077b"), - hplc$lat, hplc$lat)
hplc$lon <- ifelse(hplc$id %in% c("lovbio077b"), - hplc$lon, hplc$lon)
first_profiles$date[first_profiles$lovbio == "lovbio057b"] <- date("2015-01-20")
first_profiles$date[first_profiles$lovbio == "lovbio043b"] <- date("2015-01-20")
#those two lines correct the date because the hplc sampling was done on this exact position at this day


a <- 1
for (i in 2:length(hplc$depth)){
  hplc[i, "profile"] <- ifelse(hplc[i,"depth"]>=hplc[i-1, "depth"], a+1, a)
  a <- ifelse(hplc[i,"depth"]>=hplc[i-1, "depth"], a+1, a)
}

hplc[1, "profile"] <- 1



lovbio <- unique(hplc$id)
lovbio <- lovbio[lovbio %in% unique(first_profiles$lovbio)]

hplc$lon_round <- round(hplc$lon)
hplc$lat_round <- round(hplc$lat)
merged_dist <- merged_t
merged_dist <- merged_dist[-c(1:14),]
for (i in lovbio){
  t1 <- filter(first_profiles, lovbio == i)
  t2 <- filter(hplc)
  hplc_point <- cbind(t2$lon, t2$lat)
  argo_point <- cbind(rep(unique(t1$lon), length(hplc_point[,1])), rep(unique(t1$lat), length(hplc_point[,1])))
  distance <- distHaversine(hplc_point, argo_point)
  t2_lonround <- t2$lon_round[which(distance == min(distance))]
  t2_latround <- t2$lat_round[which(distance == min(distance))]
  t2 <- filter(t2, lon_round == unique(t2_lonround) & lat_round == unique(t2_latround))
  merged_t <- left_join(t2,t1, by = c("depth", "id" = "lovbio"))
  merged_dist <- bind_rows(merged_dist, merged_t)
}
merged_dist <- filter(merged_dist, is.na(merged_dist$chla) == FALSE)
merged_dist <- filter(merged_dist, chla_qc != 4)


merged_dist <- merged_dist %>% group_by(depth, id, date.y) %>% summarise_all(mean)




merged_dist$lag <- merged_dist$date.y - merged_dist$new_date

ggplot(filter(merged_dist, chla_qc < 4))+
  geom_point(aes(x = tchla, y = chla, colour = as.numeric(lag)))

merged_dist_clean <- filter(merged_dist, lag < 3 & lag > -3)

ggplot(merged_dist_clean)+
  geom_point(aes(x = lon.x, y = lat.x, colour = "hplc"), size = 3)+
  geom_point(aes(x = lon.y, y = lat.y, colour = "argo"))+
  geom_polygon(aes(x = long, y = lat, group = group), data = map_vec)+
  coord_quickmap()


merged_dist_clean <- merged_dist_clean[order(merged_dist_clean$id),]
merged_dist_clean$profile <- profile_numb(merged_dist_clean$depth, "downward")
merged_dist_clean <- filter(merged_dist_clean, profile != 36)
#Takuvik data ####

hplc_tak <- read_csv("Scripts/Data/hplc_tak")
hplc_tak$lon_round <- round(hplc_tak$lon)
hplc_tak$lat_round <- round(hplc_tak$lat)
position_tak <- cbind(hplc_tak$lon, hplc_tak$lat)
tak_profiles <- first_profiles[grep("takapm", first_profiles$lovbio),]

merged_dist2 <- merged_t
merged_dist2 <- merged_dist2[-c(1:22),]

for (i in unique(tak_profiles$lovbio)){
  t <- filter(tak_profiles, lovbio == i)
  t_position <- cbind(rep(unique(t$lon), length(position_tak[,1])),  rep(unique(t$lat), length(position_tak[,1])))
  distance <- distHaversine(position_tak, t_position)
  tak_lonround <- hplc_tak$lon_round[which(distance == min(distance))]
  tak_latround <- hplc_tak$lat_round[which(distance == min(distance))]
  t2 <- filter(hplc_tak, lon_round == tak_lonround & lat_round == tak_latround)
  merged_t <- left_join(t2,t, by = c("depth"))
  merged_dist2 <- bind_rows(merged_dist2, merged_t)
}

merged_dist2 <- merged_dist2 %>% group_by(depth, lovbio) %>% summarise_all(mean) %>% ungroup()
merged_dist2 <- filter(merged_dist2, is.na(chla) == FALSE)
merged_dist2 <- filter(merged_dist2, chla_adjusted_qc != 4)

ggplot(merged_dist2)+
  geom_point(aes(x = tchla, y = chla))+
  ylim(0,2.5)

#Soclim####

hplc_soclim <- read_csv("Scripts/Data/hplc_soclim")
hplc_soclim$date <- date(hplc_soclim$date)
hplc_soclim$depth <- round(hplc_soclim$depth)

profiles_soclim <- filter(first_profiles, lovbio %in% ref_bis$lovbio)
profiles_soclim$depth <- round(profiles_soclim$pres)

merged_soclim <- left_join(hplc_soclim, profiles_soclim, by = c("date", "depth"))

merged_soclim <- merged_soclim %>% filter(is.na(chla) == FALSE) %>% group_by(lovbio, depth, id.x) %>% summarise_all(mean) %>% ungroup()
merged_soclim <- filter(merged_soclim, lovbio %in% ref_bis$lovbio & lovbio != "lovbio103c" | lovbio == "lovbio103c" & lon.x > 72) 

ggplot(merged_soclim)+
  geom_point(aes(x = tchla, y = chla))

#mobydick####

hplc_mduf <- read_csv("Scripts/Data/argo/hplc_mduf")
hplc_mduf$date <- date(hplc_mduf$date) +1

merged_mduf <- left_join(hplc_mduf, filter(first_profiles, lovbio == "lovbio111b"), by = c("date", "depth")) %>% filter(is.na(chla) == FALSE)

merged_mduf <- merged_mduf %>% group_by(depth, date, lovbio) %>% summarise_all(mean)

ggplot(merged_mduf)+
  geom_point(aes(x = tchla, y = chla))

#merge all database####


merged_data <- select(merged_dist_clean, new_date, id, depth, lon.x, lat.x, lon.y, lat.y, pigments, tchla, chla, chla_qc, chla_adjusted, chla_adjusted_qc)
merged_data2 <- select(merged_dist2, date.y, lovbio, depth, lon.x, lat.x, lon.y, lat.y, pigments, tchla, chla, chla_qc, chla_adjusted, chla_adjusted_qc)

merged_data <- rename_(merged_data, "date" = "new_date", "lovbio" = "id")
merged_data2 <- rename_(merged_data2, "date" = "date.y")

merged_soclim <- select(merged_soclim, date, lovbio, depth, lon.x, lat.x, lon.y, lat.y, pigments, tchla, chla, chla_qc, chla_adjusted, chla_adjusted_qc)

merged_mduf <- select(merged_mduf, date, lovbio, depth, lon.x, lat.x, lon.y, lat.y, pigments, tchla, chla, chla_qc, chla_adjusted, chla_adjusted_qc)

merged_full <- bind_rows(merged_data, merged_data2)
merged_full <- bind_rows(merged_full, merged_soclim)
merged_full <- bind_rows(merged_full, merged_mduf)




merged_full <- merged_full[order(merged_full$lovbio),]
merged_full$profile <- profile_numb(merged_full$depth, "downward")

ggplot(merged_full)+
  geom_point(aes(x = lon.x, y = lat.x, colour = "hplc"), size = 3)+
  geom_point(aes(x = lon.y, y = lat.y, colour = "argo"))+
  geom_polygon(aes(x = long, y = lat, group = group), data = map_vec)+
  coord_quickmap()




momadata <- merged_full %>% filter(lovbio != "lovapm002a" & lovbio != "lovbio067c" & profile != 36) %>% select(depth, tchla, profile, lovbio)
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
                                        amicro = 0.0164 * exp(0.79*(depth/ze)),
                                        anano = 0.0876 * exp(-0.45*(depth/ze)),
                                        apico = 0.1393 * exp(-0.69*(depth/ze)),
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

write_csv(merged_argo, "Scripts/Data/merged_argo")



