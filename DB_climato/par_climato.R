library(ncdf4)
library(tidyverse)
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

lon <- lon2

mld <- mld2

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


#conclusion ####
ggplot()+
  geom_point(aes(x = lon, y = lat, colour = mld), data = filter(hplc, mld < 150 & dataset == "lov"))+
  geom_polygon(aes(x = long, y = lat, group = group), data = map)

usuable_hplc <- filter(hplc, dataset == "lov" & rrs667 != "NA" & mld < 200) #define the points where we can compute all we need
#we exclude maredat information because it's mainly surface points


ggplot()+
  geom_point(aes(x = lon, y = lat, colour = rrs555), data = usuable_hplc)+
  geom_polygon(aes(x = long, y = lat, group = group), data = map)

#write_csv(usuable_hplc, "DB_climato/Data/database_final")
