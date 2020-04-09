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

colnames(df) <- c('station', 'depth', paste('x', lambda, sep = '_'))
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
soclim_full[7,c(1,2)]




lov <- read_excel("Dataset_LOV.xls", na = "NA") %>% clean_names()
