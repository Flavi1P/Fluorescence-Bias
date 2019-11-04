library(ncdf4)
library(tidyverse)
library(maps)
library(fields)

files <- list.files("DB_climato/Data/PAR")
hplc <- read_csv("DB_climato/Data/global_hplc")


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




