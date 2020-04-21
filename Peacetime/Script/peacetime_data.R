library('tidyverse')
library('readxl')
library('janitor')
source('functions/zeu_moma.R')


#open data ####

peacetime_pig <- read_excel("Peacetime/PEACETIME_pigmentsLOT1et2_170218_CTD.xlsx", 
                            na = "NaN") %>% clean_names()

peacetime_aphy <- read_delim("Peacetime/database_peacetime_aphy.csv", 
                             ";", escape_double = FALSE, col_names = FALSE, 
                             trim_ws = TRUE)

wl <- read_csv("Peacetime/wl_abs_spectra.csv", 
                      col_names = FALSE)

types <- c('c', 'c', rep('n', 186))
lov <- read_csv('Data/Absorbtion/lov_soclim.csv', col_types = as.list(types))

#attribute convenient colnames
colnames(peacetime_aphy) <- c('station', 'ctd', 'lon', 'lat', 'depth', 'tchla', paste('x', wl$X1, sep = ''))
plot_aphy <- peacetime_aphy %>% pivot_longer(7:507, names_to = 'wl', values_to = 'a')

#plot to check the good attribution of wavelength
plot_aphy$wl <- as.numeric(substr(plot_aphy$wl, 2, 4))
ggplot(filter(plot_aphy, station == 'S02'))+
  geom_path(aes(x = wl, y = a, group = depth))


#check if we have the same stations
which(!unique(peacetime_aphy$station)%in% unique(peacetime_pig$station))
which(!unique(peacetime_pig$station) %in% unique(peacetime_aphy$station))

#Change FAS to FAST
peacetime_aphy$station[peacetime_aphy$station == 'FAS'] <- 'FAST'

#arrange names of pigment table to match lov ref I predefined
peacetime_pig_sel <- peacetime_pig %>% 
  mutate(t_chlc = chlorophyll_c3 + chlorophyll_c1_c2,
         t_chlb = chlorophyll_b + divinyl_chlorophyll_b,
         a_car = NA,
         b_car = NA,
         dp = fucoxanthin + alloxanthin + x19_hexanoyloxyfucoxanthin + x19_butanoyloxyfucoxanthin + zeaxanthin + chlorophyll_b + peridinin) %>% 
  select(
                           'campagne' = cruise,
                           station,
                           ctd,
                           'lon' = longitude,
                           'lat' = latitude,
                           depth,
                           date,
                           'chla'= chlorophyll_a,
                           'dv_chla' = divinyl_chlorophyll_a,
                           't_chla' = total_chlorophyll_a,
                           'phaeo' = phaeophytin_a,
                           'chlc3' = chlorophyll_c3,
                           'chlc1c2' = chlorophyll_c1_c2,
                           t_chlc,
                           'chlb' = chlorophyll_b,
                           'dv_chlb' = divinyl_chlorophyll_b,
                           t_chlb,
                           'peri' = peridinin,
                           'fuco' = fucoxanthin,
                           'x19hf' = x19_hexanoyloxyfucoxanthin,
                           'x19bf' = x19_butanoyloxyfucoxanthin,
                           'prasino' = prasinoxanthin,
                           'viola_neo' = violaxanthin,
                           'allo' = alloxanthin,
                           'zea' = zeaxanthin,
                           'lut' = lutein,
                           'diato' = diatoxanthin,
                           'diadino' = diadinoxanthin,
                           a_car,
                           b_car,
                           'tot_car' = sum_carotenes,
                           dp)

#NaN, in the raw file means LOD not NA
peacetime_pig_sel[is.na(peacetime_pig_sel)] <- 0

#round depth of pigment all 5m, to match abs data
peacetime_pig_sel$depth <- round_any(peacetime_pig_sel$depth, 5)



#error in attribution of station names of peacetime aphy
peacetime_aphy$station[peacetime_aphy$station == 'S08'] <- 'S10'
peacetime_aphy$station[peacetime_aphy$lon == 16.631 & peacetime_aphy$lat == 36.2103] <- 'S08'

#create a unique id for station, ctd, depth
peacetime_aphy$ctd <- as.numeric(peacetime_aphy$ctd)
peacetime_pig_sel <- mutate(peacetime_pig_sel, id = paste(station, ctd, depth, sep = '_'))
peacetime_aphy <- mutate(peacetime_aphy, id = paste(station, ctd, depth, sep = '_'))

unique(peacetime_aphy$id)[which(!unique(peacetime_aphy$id) %in% unique(peacetime_pig_sel$id))]


peacetime_join <- left_join(peacetime_pig_sel, peacetime_aphy, by = c('lon', 'lat', 'depth', 'ctd', 'station'))
peacetime_join <- filter(peacetime_join, x400 != 'NA')

wl_list <- colnames(select(lov, x400:x700))

peacetime_join <- select(peacetime_join, - date, -'id.x', -'id.y', wl_list)

peacetime_join$ze <- NA
for(i in unique(peacetime_join$ctd)){
  df <- filter(peacetime_join, ctd  == i)
  zeu <- Zeu_moma(df$chla, df$depth)
  peacetime_join[peacetime_join$ctd == i, 'ze'] <- zeu
}
unique(peacetime_join$ze)

peacetime_ready <- peacetime_join %>% mutate(z_zeu = depth/ze,
                                             campagne = 'peacetime',
                                             wdp = 1.56 * fuco + 0.92 * peri + 4.8 * allo + 1.02 * x19bf + 1.12 * x19hf + 1.51 * zea + 0.69 * t_chlb,
                                             p_micro = (1.56 * fuco + 0.92 * peri)/wdp,
                                             p_nano = (4.8 * allo + 1.02 * x19bf + 1.51 * x19hf)/wdp,
                                             p_pico = (1.51 * zea + 0.69 * t_chlb)/wdp,
                                             pico_t_chla = p_pico * t_chla,
                                             nano_t_chla =  p_nano * t_chla,
                                             micro_t_chla = p_micro * t_chla) %>% 
  select(campagne,
         station,
         lon,
         lat, 
         depth,
         z_zeu,
         chla,
         dv_chla,
         t_chla,
         phaeo,
         chlc3,
         chlc1c2,
         t_chlc,
         chlb,
         dv_chlb,
         t_chlb,
         peri,
         fuco,
         x19hf,
         x19bf,
         prasino,
         viola_neo,
         allo,
         zea,
         lut,
         diato,
         diadino,
         a_car,
         b_car,
         tot_car,
         dp,
         p_pico,
         p_nano,
         p_micro,
         pico_t_chla,
         nano_t_chla,
         micro_t_chla,
         wl_list)

lov_full <- bind_rows(lov, peacetime_ready)

#write_csv(lov_full, 'Data/absorbtion/lov_soclim_peacetime.csv')
