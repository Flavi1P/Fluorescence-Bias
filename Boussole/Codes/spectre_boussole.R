library(tidyverse)
library(ggfortify)
library(castr)
library(stlplus)
library(janitor)
library(gridExtra)
library(lubridate)

#open boussole data
boussole <- read_csv("Boussole/Data/boussole.csv")

#open pigments absorbtion
spectre <- read_excel("Biosope/Data/Spectres_annick.xlsx")
spectre <- clean_names(spectre)

#select the absorbtion at 440 and 470nm
spectre440<- filter(spectre, lambda == 440)
spectre470 <- filter(spectre, lambda == 470)

#create columns that correspond to the total photosynthetic absorbance and non photosynthetic absorbance at 440 and 470. Create also a ratio between the two photosynthetic absorbtion
boussole <- boussole %>% mutate(photo_440 = peri * spectre440$peri + but * spectre440$x19_bf + hex * spectre440$x19_hf + fuco * spectre440$fuco + allo * spectre440$allox + dvchla * spectre440$dv_chla + chla * spectre440$chl_a,
                              protect_440 = zea * spectre440$zea,
                              photo_470 = peri * spectre470$peri + but * spectre470$x19_bf + hex * spectre470$x19_hf + fuco * spectre470$fuco + allo * spectre470$allox + dvchla * spectre470$dv_chla + chla * spectre470$chl_a,
                              protect_470 = zea * spectre470$zea,
                              ratio = photo_440/photo_470)

#plot the ratio as function of tchla

ggplot(filter(boussole, depth < 20))+
  geom_point(aes(x = tchla, y = ratio))+
  theme_bw()

#plot the ratio variation in time

ggplot(filter(boussole, depth < 20))+
  geom_point(aes(x = date, y = ratio))+
  theme_bw()

#ease the visualistation by doing the mean of each date

boussole_mean <- boussole %>% filter(depth <20) %>%  group_by(date) %>% summarise_all(mean)

ggplot(filter(boussole_mean, depth < 20))+
  geom_point(aes(x = date, y = ratio))+
  theme_bw()

#regularisation of the time serie

#see the difference between dates
difference <- as.numeric(diff(boussole_mean$date))
hist(difference)

#create a new date scale
date_new <- seq(min(boussole_mean$date), max(boussole_mean$date), by = 20)

#interpolate values on this new date scale
interp_values <- interpolate(boussole_mean$date, boussole_mean$ratio, date_new, method = "linear")
interp_values_tchla <- interpolate(boussole_mean$date, boussole_mean$tchla, date_new, method = "linear")


#create a df on the new scale with new values
boussole_ts <- data.frame("date" = date_new, "ratio" = interp_values, "tchla" = interp_values_tchla)

boussole_ts$year <- year(boussole_ts$date)

ratio_time <- ggplot(boussole_ts)+
  geom_path(aes(x = date, y = ratio), colour = "#fdae6b")+
  facet_wrap(.~year, ncol = 1, scales = "free")+
  theme_bw()+
  ylab("aPS440/aPS470")

tchla_time <- ggplot(boussole_ts)+
  geom_path(aes(x = date, y = tchla), colour = "#a1d99b")+
  facet_wrap(.~year, ncol = 1, scales = "free")+
  theme_bw()+
  ylab("[Chla]")



grid.arrange(tchla_time, ratio_time, ncol = 2)
#create a ts object
boussole_ts <- ts(boussole_ts$ratio)


#stl on ts (loess decomposition)

boussole_stl <- stlplus(boussole_ts, n.p = 4, s.window = "periodic", t.windows = 60)

#plot the object

plot(boussole_stl, scales = list(y = "free"))


#This waas not very good
#Now we try to see the variability of pigments composition

#Create a df with the date and pigment concentration per month at the surface
pigment_ts <- boussole_mean %>% select(date, fuco, peri, hex, but, allo, dvchla, chla)%>% 
  mutate(month = month(date), year = year(date)) %>%
  group_by(month, year) %>% 
  select(-date) %>% 
  summarise_all(mean) %>% 
  ungroup() 

#change the format to allowed a visualisatio with geom_bar
pigment_ts <- gather(pigment_ts, key = "pigment", value = "concentration", 3:9)

#plot the variation of different pigments along the years
pigment_time <- ggplot(filter(pigment_ts, pigment != "tchla"))+
  geom_bar(aes(x = month, y = concentration, fill = pigment), stat = "identity", position = "fill")+
  scale_fill_brewer(palette = "Set3")+
  theme_dark()+
  facet_wrap(. ~ year, ncol = 1)

#compare the variation of pigments with the variation of the ratio
grid.arrange(ratio_time, pigment_time, ncol = 2)

#plot the ratio a440/a470 vs fluo/chla

boussole <- boussole %>% mutate(ratio_abs = photo_440/photo_470) %>% filter(ratio_abs < 40)

ggplot(boussole)+
  geom_point(aes(x = tchla, y = ratio_abs))

             