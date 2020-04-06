library(tidyverse)
library(readxl)
library(lubridate)
source("functions/profile_numb.R")
source("functions/zeu_moma.R")

pigments <- c("fuco", "peri", "hex", "but", "allo", "tchlb", "zea")
hplc <- read_excel("Data/Donnees-pigments-campagnes.xlsx", na = "NA")

#Boussole 2014####

boussole_2014 <- read_excel("Data/Raw/Full/BDD_HPLC_fluo_cruises/Boussole-btl.xlsx", 
                            sheet = "2014", na = "NA")

boussole_2014 <- boussole_2014 %>% mutate(date = as_date(paste(Date, ...16, ...17, sep = "-")),
                                          depth = round(Press)) %>% 
  select(btle, date, Press, Fluo, depth, "F[volts]", "T[ITS90]")

names(boussole_2014) <- c("btle", "date", "press", "fluo", "depth", "fluo_volt", "temp")



hplc_boussole_2014 <- filter(hplc, Project == "Boussole" & Year == "2014")
names(hplc_boussole_2014) <- tolower(colnames(hplc_boussole_2014))
hplc_boussole_2014 <- hplc_boussole_2014 %>% mutate(date = as_date(paste(year, month, day, sep = "-"))) %>% select(date, depth, bottle_number, pigments, tchla, phytinea, dvchlb, dvchla, chla)

hplc_boussole_2014$bottle_number <- as.numeric(hplc_boussole_2014$bottle_number)
boussole_2014 <- left_join(boussole_2014, hplc_boussole_2014, by = c("date", "btle" = "bottle_number", "depth"))


#Boussole 2015 ####

boussole_2015 <- read_excel("Data/Raw/Full/BDD_HPLC_fluo_cruises/Boussole-btl.xlsx", 
                            sheet = "2015", na = "NA")

boussole_2015 <- boussole_2015 %>% mutate(date = as_date(paste(Date, ...16, ...17, sep = "-")),
                                          depth = round(Press)) %>% 
  select(btle, date, Press, Fluo, depth, "F[volts]", "T[ITS90]")

names(boussole_2015) <- c("btle", "date", "press", "fluo", "depth", "fluo_volt", "temp")



hplc_boussole_2015 <- filter(hplc, Project == "Boussole" & Year == "2015")
names(hplc_boussole_2015) <- tolower(colnames(hplc_boussole_2015))
hplc_boussole_2015 <- hplc_boussole_2015 %>% mutate(date = as_date(paste(year, month, day, sep = "-"))) %>% select(date, depth, bottle_number, pigments, tchla, phytinea, dvchlb, dvchla, chla)

hplc_boussole_2015$bottle_number <- as.numeric(hplc_boussole_2015$bottle_number)
boussole_2015 <- left_join(boussole_2015, hplc_boussole_2015, by = c("date", "btle" = "bottle_number", "depth"))

#boussole 2016 ####

boussole_2013 <- read_excel("Data/Raw/Full/BDD_HPLC_fluo_cruises/Boussole-btl.xlsx", 
                            sheet = "2013", na = "NA")

boussole_2013 <- boussole_2013 %>% mutate(date = as_date(paste(Date, ...18, ...19, sep = "-")),
                                          depth = round(Press)) %>% 
  select(btle, date, Press, Fluo, depth, "F[volts]", "T[ITS90]")

names(boussole_2013) <- c("btle", "date", "press", "fluo", "depth", "fluo_volt", "temp")



hplc_boussole_2013 <- filter(hplc, Project == "Boussole" & Year == "2013")
names(hplc_boussole_2013) <- tolower(colnames(hplc_boussole_2013))
hplc_boussole_2013 <- hplc_boussole_2013 %>% mutate(date = as_date(paste(year, month, day, sep = "-"))) %>% select(date, depth, bottle_number, pigments, tchla, phytinea, dvchlb, dvchla, chla)

hplc_boussole_2013$bottle_number <- as.numeric(hplc_boussole_2013$bottle_number)
boussole_2013 <- left_join(boussole_2013, hplc_boussole_2013, by = c("date", "btle" = "bottle_number", "depth"))

#binds 3 years ####

boussole <- bind_rows(boussole_2013, boussole_2014, boussole_2015)
rm(hplc_boussole_2013, hplc_boussole_2014, hplc_boussole_2015, hplc, boussole_2013, boussole_2014, boussole_2015)

#add ze####
boussole$profile <- NA
a <- 1
for (i in c(2 : length(boussole$date))){
  boussole$profile[i] <- ifelse(boussole$date[i] == boussole$date[i-1], a, a+1)
  a <- ifelse(boussole$date[i] == boussole$date[i-1], a, a+1)
}
boussole$profile[1] <- 1

boussole <- filter(boussole, tchla != "NA")
momadata <- boussole

momadata$ze <- NA

momadata <- filter(boussole, profile %in% which(table(momadata$profile) > 1)) %>%  select(depth, tchla, profile)

for (i in momadata$profile){
  dat <- filter(momadata, profile == i)
  momadata$ze[which(momadata$profile == i)] <- Zeu_moma(dat$tchla, dat$depth)
}

momadata <- momadata %>% group_by(profile) %>% summarize_all(mean)
boussole <- left_join(boussole, select(momadata, ze, profile), by ="profile")


#add other variables####

boussole <- boussole %>% mutate(wdp = 1.56 * fuco + 0.92 * peri + 4.8 * allo + 1.02 * but + 1.12 * hex + 1.51 * zea + 0.69 * tchlb,
                              micro = (1.56 * fuco + 0.92 * peri)/wdp,
                              nano = (4.8 * allo + 1.02 * but + 1.51 * hex)/wdp,
                              pico = (1.51 * zea + 0.69 * tchlb)/wdp,
                              amicro = 0.0164 * exp(0.79*(depth/ze)),
                              anano = 0.0876 * exp(-0.45*(depth/ze)),
                              apico = 0.1393 * exp(-0.69*(depth/ze)),
                              microfluo = amicro * micro * tchla,
                              nanofluo = anano * nano * tchla,
                              picofluo = apico * pico * tchla,
                              optical_layer = (depth/ze)/0.50001 + 1,
                              ratio = fluo/ tchla)
#wrtite csv####
#write_csv(boussole, "Boussole/Data/boussole.csv")
