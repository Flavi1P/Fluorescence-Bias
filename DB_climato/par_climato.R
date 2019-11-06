library(ncdf4)
library(tidyverse)
library(lubridate)
library(maps)
library(fields)

hplc <- read_csv("DB_climato/Data/global_hplc")
map <- read_csv("Data/map_vec")

# PAR values ####

files <- list.files("DB_climato/Data/PAR")



files_list <- data.frame("month" = NA, "file" = NA) #create an empty df
for(i in files){
  month <- round(as.numeric(substr(i, 6, 8))/30)+1 #compute the month of the data from the climatology file
  files_list_temp <- data.frame("month" = month, "file"= i) #associate the name of the file with the month of data
  files_list <- bind_rows(files_list, files_list_temp)
}

files_list <- files_list[-1,] #erase the na needed to create the ddf (ugly method)

hplc$par <- NA

for(i in 1:nrow(hplc)){
  month_of_sample <- hplc$month[i] #save the month of the hplc sample
  lon_of_sample <- hplc$lon[i] #the lon
  lat_of_sample <- hplc$lat[i] #the lat
  par_file <- files_list[files_list$month == month_of_sample, 2] #take the nc data at the month of interest
  nc_sample <- nc_open(paste("DB_climato/Data/PAR/", par_file, sep = "")) #open it
  #extract values at the colsest location of our sample
  par_sample <- ncvar_get(nc_sample, varid = 'par',
                          start= c(which.min(abs(nc_sample$dim$lon$vals - lon_of_sample)), # look for closest long
                                   which.min(abs(nc_sample$dim$lat$vals - lat_of_sample))),  # look for closest lat)
                          count = c(1,1)) 
  hplc$par[i] <- par_sample #attribute the value in the ncdf
  nc_close(nc_sample)
}



# rrs 667 ####

files <- list.files("DB_climato/Data/rrs_667")



files_list <- data.frame("month" = NA, "file" = NA) #create an empty df
for(i in files){
  month <- round(as.numeric(substr(i, 6, 8))/30)+1 #compute the month of the data from the climatology file
  files_list_temp <- data.frame("month" = month, "file"= i) #associate the name of the file with the month of data
  files_list <- bind_rows(files_list, files_list_temp)
}

files_list <- files_list[-1,] #erase the na needed to create the ddf (ugly method)

hplc$rrs667 <- NA
pb <- txtProgressBar(min = 0, max = nrow(hplc), style = 3)

for(i in 1:nrow(hplc)){
  month_of_sample <- hplc$month[i] #save the month of the hplc sample
  lon_of_sample <- hplc$lon[i] #the lon
  lat_of_sample <- hplc$lat[i] #the lat
  rrs667_file <- files_list[files_list$month == month_of_sample, 2] #take the nc data at the month of interest
  nc_sample <- nc_open(paste("DB_climato/Data/rrs_667/", rrs667_file, sep = "")) #open it
  #extract values at the colsest location of our sample
  rrs667_sample <- ncvar_get(nc_sample, varid = 'Rrs_667',
                          start= c(which.min(abs(nc_sample$dim$lon$vals - lon_of_sample)), # look for closest long
                                   which.min(abs(nc_sample$dim$lat$vals - lat_of_sample))),  # look for closest lat)
                          count = c(1,1)) 
  hplc$rrs667[i] <- rrs667_sample #attribute the value in the ncdf
  nc_close(nc_sample)
  setTxtProgressBar(pb, i)
}

table(is.na(hplc$rrs667))

# rrs 555 ####

files <- list.files("DB_climato/Data/rrs_555")



files_list <- data.frame("month" = NA, "file" = NA) #create an empty df
for(i in files){
  month <- round(as.numeric(substr(i, 6, 8))/30)+1 #compute the month of the data from the climatology file
  files_list_temp <- data.frame("month" = month, "file"= i) #associate the name of the file with the month of data
  files_list <- bind_rows(files_list, files_list_temp)
}

files_list <- files_list[-1,] #erase the na needed to create the ddf (ugly method)

hplc$rrs555 <- NA
pb <- txtProgressBar(min = 0, max = nrow(hplc), style = 3)

for(i in 1:nrow(hplc)){
  month_of_sample <- hplc$month[i] #save the month of the hplc sample
  lon_of_sample <- hplc$lon[i] #the lon
  lat_of_sample <- hplc$lat[i] #the lat
  rrs555_file <- files_list[files_list$month == month_of_sample, 2] #take the nc data at the month of interest
  nc_sample <- nc_open(paste("DB_climato/Data/rrs_555/", rrs555_file, sep = "")) #open it
  #extract values at the colsest location of our sample
  rrs555_sample <- ncvar_get(nc_sample, varid = 'Rrs_555',
                             start= c(which.min(abs(nc_sample$dim$lon$vals - lon_of_sample)), # look for closest long
                                      which.min(abs(nc_sample$dim$lat$vals - lat_of_sample))),  # look for closest lat)
                             count = c(1,1)) 
  hplc$rrs555[i] <- rrs555_sample #attribute the value in the ncdf
  nc_close(nc_sample)
  setTxtProgressBar(pb, i)
}

table(is.na(hplc$rrs555))

#save the file 
#write_csv(hplc, "DB_climato/Data/db_save")

# rs 488 ####

files <- list.files("DB_climato/Data/rrs_488")



files_list <- data.frame("month" = NA, "file" = NA) #create an empty df
for(i in files){
  month <- round(as.numeric(substr(i, 6, 8))/30)+1 #compute the month of the data from the climatology file
  files_list_temp <- data.frame("month" = month, "file"= i) #associate the name of the file with the month of data
  files_list <- bind_rows(files_list, files_list_temp)
}

files_list <- files_list[-1,] #erase the na needed to create the ddf (ugly method)

hplc$rrs488 <- NA
pb <- txtProgressBar(min = 0, max = nrow(hplc), style = 3)

for(i in 1:nrow(hplc)){
  month_of_sample <- hplc$month[i] #save the month of the hplc sample
  lon_of_sample <- hplc$lon[i] #the lon
  lat_of_sample <- hplc$lat[i] #the lat
  rrs488_file <- files_list[files_list$month == month_of_sample, 2] #take the nc data at the month of interest
  nc_sample <- nc_open(paste("DB_climato/Data/rrs_488/", rrs488_file, sep = "")) #open it
  #extract values at the colsest location of our sample
  rrs488_sample <- ncvar_get(nc_sample, varid = 'Rrs_488',
                             start= c(which.min(abs(nc_sample$dim$lon$vals - lon_of_sample)), # look for closest long
                                      which.min(abs(nc_sample$dim$lat$vals - lat_of_sample))),  # look for closest lat)
                             count = c(1,1)) 
  hplc$rrs488[i] <- rrs488_sample #attribute the value in the ncdf
  nc_close(nc_sample)
  setTxtProgressBar(pb, i)
}  

table(is.na(hplc$rrs488))

# rs 443 ####

files <- list.files("DB_climato/Data/rrs_443")



files_list <- data.frame("month" = NA, "file" = NA) #create an empty df
for(i in files){
  month <- round(as.numeric(substr(i, 6, 8))/30)+1 #compute the month of the data from the climatology file
  files_list_temp <- data.frame("month" = month, "file"= i) #associate the name of the file with the month of data
  files_list <- bind_rows(files_list, files_list_temp)
}

files_list <- files_list[-1,] #erase the na needed to create the ddf (ugly method)

hplc$rrs443 <- NA
pb <- txtProgressBar(min = 0, max = nrow(hplc), style = 3)

for(i in 1:nrow(hplc)){
  month_of_sample <- hplc$month[i] #save the month of the hplc sample
  lon_of_sample <- hplc$lon[i] #the lon
  lat_of_sample <- hplc$lat[i] #the lat
  rrs443_file <- files_list[files_list$month == month_of_sample, 2] #take the nc data at the month of interest
  nc_sample <- nc_open(paste("DB_climato/Data/rrs_443/", rrs443_file, sep = "")) #open it
  #extract values at the colsest location of our sample
  rrs443_sample <- ncvar_get(nc_sample, varid = 'Rrs_443',
                             start= c(which.min(abs(nc_sample$dim$lon$vals - lon_of_sample)), # look for closest long
                                      which.min(abs(nc_sample$dim$lat$vals - lat_of_sample))),  # look for closest lat)
                             count = c(1,1)) 
  hplc$rrs443[i] <- rrs443_sample #attribute the value in the ncdf
  nc_close(nc_sample)
  setTxtProgressBar(pb, i)
}  

table(is.na(hplc$rrs443))

# rs 412 ####

files <- list.files("DB_climato/Data/rrs_412")



files_list <- data.frame("month" = NA, "file" = NA) #create an empty df
for(i in files){
  month <- round(as.numeric(substr(i, 6, 8))/30)+1 #compute the month of the data from the climatology file
  files_list_temp <- data.frame("month" = month, "file"= i) #associate the name of the file with the month of data
  files_list <- bind_rows(files_list, files_list_temp)
}

files_list <- files_list[-1,] #erase the na needed to create the ddf (ugly method)

hplc$rrs412 <- NA
pb <- txtProgressBar(min = 0, max = nrow(hplc), style = 3)

for(i in 1:nrow(hplc)){
  month_of_sample <- hplc$month[i] #save the month of the hplc sample
  lon_of_sample <- hplc$lon[i] #the lon
  lat_of_sample <- hplc$lat[i] #the lat
  rrs412_file <- files_list[files_list$month == month_of_sample, 2] #take the nc data at the month of interest
  nc_sample <- nc_open(paste("DB_climato/Data/rrs_412/", rrs412_file, sep = "")) #open it
  #extract values at the colsest location of our sample
  rrs412_sample <- ncvar_get(nc_sample, varid = 'Rrs_412',
                             start= c(which.min(abs(nc_sample$dim$lon$vals - lon_of_sample)), # look for closest long
                                      which.min(abs(nc_sample$dim$lat$vals - lat_of_sample))),  # look for closest lat)
                             count = c(1,1)) 
  hplc$rrs412[i] <- rrs412_sample #attribute the value in the ncdf
  nc_close(nc_sample)
  setTxtProgressBar(pb, i)
}  

table(is.na(hplc$rrs412))


#MLD ####

nc_mld <- nc_open("DB_climato/Data/mld/mld_DR003_c1m_reg2.0.nc")

nc_mld$var$mld$missval
mld <- ncvar_get(nc_mld, nc_mld$var$mld)
lon <- ncvar_get(nc_mld, nc_mld$dim$lon)
lat <- ncvar_get(nc_mld, nc_mld$dim$lat)

lon2 <- lon
for(i in 1:length(lon)){
  if(lon[i]>=0 & lon[i]<180){
    lon2[i] <- lon[i]
  }
  else{
    lon2[i] <- lon[i]-360
  }
}

order_lon <- order(lon2)

mld2 <- mld[order_lon,,]
lon2 <- lon2[order_lon]
lon <- lon2

mld <- mld2
# image.plot(lon,lat,mld[,,3],zlim=c(0,900))
# map(add=T)

hplc$mld <- NA
pb <- txtProgressBar(min = 0, max = nrow(hplc), style = 3)

for(i in 1:nrow(hplc)){
  month_of_sample <- hplc$month[i] #save the month of the hplc sample
  lon_of_sample <- hplc$lon[i] #the lon
  lat_of_sample <- hplc$lat[i] #the lat
  mld_of_sample <- mld[which.min(abs(lon - lon_of_sample)), which.min(abs(lat - lat_of_sample)), month_of_sample] 
  #extract values at the colsest location of our sample
  hplc$mld[i] <- mld_of_sample #attribute the value in the ncdf
  setTxtProgressBar(pb, i)
}  

nc_close(nc_mld)
table(is.na(hplc$mld))


#checkpoint####
# ggplot()+
#   geom_point(aes(x = lon, y = lat, colour = mld), data = filter(hplc, dataset == "lov"))+
#   geom_polygon(aes(x = long, y = lat, group = group), data = map)+
#   coord_quickmap()

usuable_hplc <- filter(hplc, dataset == "lov" & rrs667 != "NA" & mld < 1000) #define the points where we can compute all we need
#we exclude maredat information because it's mainly surface points


# ggplot()+
#   geom_point(aes(x = lon, y = lat, colour = rrs555), data = usuable_hplc)+
#   geom_polygon(aes(x = long, y = lat, group = group), data = map)

#write_csv(usuable_hplc, "DB_climato/Data/database_final")

hplc2 <- usuable_hplc %>% mutate(lon = round(lon, 2),
                lat = round(lat, 2)) #usuablee_hplc is deprecated from now

LON_LAT_MONTH <- unique(hplc2[,c("lon","lat","month")]) %>% mutate(nprof = seq(1, nrow(.), by = 1)) #create a num of profile

hplc2 <- left_join(hplc2, LON_LAT_MONTH) #add the num of the profile on hplc2
hplc2 <- filter(hplc2, nprof %in% which(table(hplc2$nprof)>4)) #reject profiles with less than 4 data

# ggplot()+
#   geom_point(aes(x = lon, y = lat, colour = mld), data = usuable_hplc)+
#   geom_polygon(aes(x = long, y = lat, group = group), data = map)+
#   coord_quickmap()


# bathymetric data ####

#we want to select only data that come from a profile which is in a zone deeper than 500m

nc_bath <- nc_open("DB_climato/Data/bathy/GEBCO_2014_6x6min_Global.nc")

bath <- ncvar_get(nc_bath, "Height")
lon <- ncvar_get(nc_bath, nc_bath$dim$lon)
lat <- ncvar_get(nc_bath, nc_bath$dim$lat)

hplc2$bath <- NA

for(i in 1:nrow(hplc2)){
  lon_of_sample <- hplc2$lon[i] #the lon
  lat_of_sample <- hplc2$lat[i] #the lat
  bath_of_sample <- bath[which.min(abs(lon - lon_of_sample)), which.min(abs(lat - lat_of_sample))] 
  #extract values at the colsest location of our sample
  hplc2$bath[i] <- bath_of_sample #attribute the value in the ncdf
}  

#ggplot()+
  # geom_point(aes(x = lon, y = lat, colour = bath), data = hplc2)+
  # geom_polygon(aes(x = long, y = lat, group = group), data = map)+
  # coord_quickmap()

hplc_open_ocean <- filter(hplc2, bath < - 500) # we filter on the bathymetric criteria

# ggplot()+
#   geom_point(aes(x = lon, y = lat, colour = bath), data = hplc_open_ocean)+
#   geom_polygon(aes(x = long, y = lat, group = group), data = map)+
#   coord_quickmap()

good_profile <- c()
for(i in unique(hplc_open_ocean$nprof)){
  t_prof <- filter(hplc_open_ocean, nprof == i)
  if(max(t_prof$depth) > 100 & min(t_prof$depth < 10)){
    good_profile <- c(good_profile, i)
  }
}




#write_csv(hplc_open_ocean, "DB_climato/Data/database_final")

#ze####

hplc <- read_csv("DB_climato/Data/database_final")
source("functions/fonction_calcul_Ze_Monte_Carlo.r")
source("functions/zeu_moma.R")



hplc$ze_monte_carlo <- NA
for(i in unique(hplc$nprof)){
  t_prof <- filter(hplc, nprof == i)
  ze <- calcul_Ze_Monte_Carlo(t_prof$depth, t_prof$chla)
  hplc[hplc$nprof == i,]$ze_monte_carlo <- ze$chlze
}

hplc$ze_morel <- NA
for(i in unique(hplc$nprof)){
  t_prof <- filter(hplc, nprof == i)
  ze <- Zeu_moma(t_prof$chla, t_prof$depth)
  hplc[hplc$nprof == i,]$ze_morel <- ze
}

ggplot()+
  geom_point(aes(x = lon, y = lat, colour = ze_morel), data = hplc)+
  geom_polygon(aes(x = long, y = lat, group = group), data = map)+
  coord_quickmap()
#the ze_morel looks pretty much better
good_profile <- c()
for(i in unique(hplc$nprof)){
  t_prof <- filter(hplc, nprof == i)
    if(is.na(t_prof$ze_morel)){}
  else{
  if(unique(t_prof$ze_morel) < max(t_prof$depth) & min(t_prof$depth) < 10){
    good_profile <- c(good_profile, i)
  }
  }
}

hplc_qc <- filter(hplc, nprof %in% good_profile) #remove 283 profiles

ggplot()+
  geom_point(aes(x = lon, y = lat, colour = ze_morel), data = hplc_qc)+
  geom_polygon(aes(x = long, y = lat, group = group), data = map)+
  coord_quickmap()
#now we have our filtered dataset 
# we add absorbtion variables

#aps variables ####

#open pigments absorbtion
spectre <- read_excel("Biosope/Data/Spectres_annick.xlsx")
spectre <- clean_names(spectre)

#select the absorbtion at 440 and 470nm
spectre440<- filter(spectre, lambda == 440)
spectre470 <- filter(spectre, lambda == 470)

hplc_qc <- hplc_qc %>% mutate(photo_440 = peri * spectre440$peri + but * spectre440$x19_bf + hex * spectre440$x19_hf + fuco * spectre440$fuco + allo * spectre440$allox + chla * spectre440$chl_a + dv_chla * spectre440$dv_chla,
                              photo_470 = peri * spectre470$peri + but * spectre470$x19_bf + hex * spectre470$x19_hf + fuco * spectre470$fuco + allo * spectre470$allox + chla * spectre470$chl_a + dv_chla * spectre470$dv_chla,
                              ratio = photo_440/photo_470)

hplc_qc$doy <- yday(hplc_qc$date)#compute day of year

hplc_qc <- select(hplc_qc, date, doy, month, everything(), -test, - profile_id, - dataset, - ze_monte_carlo, )
#we write our final csv

#write_csv(hplc_qc, "DB_climato/Data/database_final")



