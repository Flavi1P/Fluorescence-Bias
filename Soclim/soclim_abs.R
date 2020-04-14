library(tidyverse)
library(readxl)
library(janitor)

load('Soclim/Soclim_absorption/database_soclim.RData')

aph <- abs$aph
lambda <- abs$wl
z <- abs$z
z_vect <- as.vector(z)
st <- abs$st

plot(lambda, aph[1,1,])

df <- matrix(nrow = 13*12, ncol = 503)
df <- as.data.frame(df)

colnames(df) <- c('station', 'depth', paste('x', lambda, sep = ''))
df$station <- rep(st, 12)
df$depth <- z_vect

t <- 0
for(i in 1:13){
  t <- t + 1
  f <- t
  for(j in 1:12){
    df[f,3:503] <- aph[i,j,]
    f <- f + 13
  }
}

df_long <- df %>% 
  pivot_longer(x_800:x_300, names_to = 'lambda', values_to = 'aph')

df_long$lambda <- as.numeric(substr(df_long$lambda, 3,5))

ggplot(filter(df_long, station == 'O25' & depth == 40))+
  geom_path(aes(x = lambda, y = aph))

btl <- read_delim("Data/Raw/soclim/soclim_btl_V9.txt", 
                            "\t", escape_double = FALSE, trim_ws = TRUE, na = c("NA", "")) %>% clean_names()

pigment <- btl %>% select(station, latitude_deg_north, longitude_deg_east, depth_m, temperature_deg_celcius, salinity, chlorophyll_c3:total_chlorophyll_a, - contains('qa'))
pigment[pigment == 'LOD'] <- 0

st2 <- pigment$station
pigment <- pigment %>%
  select(-station) %>% 
  mutate_if(is.character, as.numeric) %>% 
  mutate('station' = st2) %>% 
  select(station, depth = depth_m, everything())

pigment_nona <- na.omit(pigment)
pigment_nona$depth <- round(pigment_nona$depth)
df$station <- gsub('-', '', df$station)

pigment_nona$z <- round(pigment_nona$depth / 5) * 5
pigment_nona <- select(pigment_nona, -depth)
soclim_full <- left_join(pigment_nona, df, by = c('station', 'z' = 'depth')) %>% 
  select(station, 'depth' = z, everything())

which(is.na(soclim_full$x_470))


lov <-read_excel("Dataset_LOV.xls", col_types = c("text", "text", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric"), 
                                 na = "NA") %>% clean_names()

lambda_sel <- colnames(select(lov, x400:x700))

colnames(soclim_full)
colnames(lov)

source('functions/zeu_moma.R')
soclim_full$ze <- NA
for(i in unique(soclim_full$station)){
  df <- filter(soclim_full, station == i)
  zeu <- Zeu_moma(df$chlorophyll_a, df$depth)
  soclim_full[soclim_full$station == i, 'ze'] <- zeu
}

size <- as.numeric(length(soclim_full$station))
soclim_full <- mutate(soclim_full, z_zeu = depth/ze,
                      campagne = 'soclim',
                      phaeo = as.numeric(rep(NA, size)),
                      t_chlc = chlorophyll_c3 + chlorophyll_c1_c2_mg_dvp,
                      t_chlb = chlorophyll_b + divinyl_chlorophyll_b,
                      a_car = as.numeric(rep(NA, size)),
                      b_car = as.numeric(rep(NA, size)),
                      tot_car = as.numeric(rep(NA, size)),
                      wdp = 1.56 * fucoxanthin + 0.92 * peridinin + 4.8 * alloxanthin + 1.02 * x19_butanoyloxyfucoxanthin + 1.12 * x19_hexanoyloxyfucoxanthin + 1.51 * zeaxanthin + 0.69 * t_chlb,
                      dp = fucoxanthin + peridinin + alloxanthin + x19_butanoyloxyfucoxanthin + x19_hexanoyloxyfucoxanthin + zeaxanthin + t_chlb,
                      p_micro = (1.56 * fucoxanthin + 0.92 * peridinin)/wdp,
                      p_nano = (4.8 * alloxanthin + 1.02 * x19_butanoyloxyfucoxanthin + 1.51 * x19_hexanoyloxyfucoxanthin)/wdp,
                      p_pico = (1.51 * zeaxanthin + 0.69 * t_chlb)/wdp,
                      pico_t_chla = p_pico * total_chlorophyll_a,
                      nano_t_chla =  p_nano * total_chlorophyll_a,
                      micro_t_chla = p_micro * total_chlorophyll_a)

soclim_sel <- soclim_full %>% select(campagne,
                                     station,
                                     lon = longitude_deg_east,
                                     lat = latitude_deg_north, 
                                     depth,
                                     z_zeu,
                                     chla = chlorophyll_a,
                                     dv_chla = divinyl_chlorophyll_a,
                                     t_chla = total_chlorophyll_a,
                                     phaeo,
                                     chlc3 = chlorophyll_c3,
                                     chlc1c2 = chlorophyll_c1_c2_mg_dvp,
                                     t_chlc,
                                     chlb = chlorophyll_b,
                                     dv_chlb = divinyl_chlorophyll_b,
                                     t_chlb,
                                     peri = peridinin,
                                     fuco = fucoxanthin,
                                     x19hf = x19_hexanoyloxyfucoxanthin,
                                     x19bf = x19_butanoyloxyfucoxanthin,
                                     prasino = prasinoxanthin,
                                     viola_neo = violaxanthin,
                                     allo = alloxanthin,
                                     zea = zeaxanthin,
                                     lut = lutein,
                                     diato = diatoxanthin,
                                     diadino = diadinoxanthin,
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
                                     lambda_sel)

lov_sel <- lov %>% select(campagne,
                          station = ctd_station,
                          lon,
                          lat,
                          depth = profondeur,
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
                          lambda_sel)



combine_df <- bind_rows(lov_sel, soclim_sel)

#write_csv(combine_df, 'Data/absorbtion/lov_soclim.csv')
