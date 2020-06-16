library(ncdf4)
library(tidyverse)
library(readxl)
library(lubridate)
library(plyr)


#download.file("ftp://ftp.ifremer.fr/ifremer/argo/argo_merge-profile_index.txt","argo/Data/index_merge.txt", quiet = FALSE, mode = "w",cacheOK = TRUE) 
#biofiles <- read_csv("argo/Data/index_bio.txt", skip = 7)
# 
# index_ifremer<-read.table("argo/Data/index_merge.txt", skip=9, sep = ",")
# files<-as.character(index_ifremer[,1])
# ident<-strsplit(files,"/")
# ident<-matrix(unlist(ident), ncol=4, byrow=TRUE)
# dac<-ident[,1]
# wod<-ident[,2]
# prof_id<-ident[,4]
# variables<-as.character(index_ifremer[,8])
# lat<-index_ifremer[,3]
# lon<-index_ifremer[,4]
# time_all<-as.character(paste(index_ifremer$V2))


profiles <- list.files("argo/Data/biofiles", full.names = TRUE)

ref <- read_csv("Data/argo/ref.csv")
ref_bis <- read_csv("Data/argo/ref_bis")
ref <- bind_rows(ref, ref_bis)


argo_data <- data.frame("date" = NA, "lat" = NA, "lon" = NA, "pres" = NA,
                        "chla" = NA, "chla_qc" = NA, "chla_adjusted" = NA, "chla_adjusted_qc" = NA, "id"= NA)

a <- 1

#suite####

pb <- txtProgressBar(min = 1, max = length(profiles), style = 3)
for(i in profiles){
  nc <- NA
  nc <- nc_open(i)
  variables <- names(nc$var) %>% data.frame %>% filter(. == "CHLA")
  juld <- ncvar_get(nc, "JULD")
  juldqc <- ncvar_get(nc, "JULD_QC")
  origin<-NA
  origin<-as.POSIXct("1950-01-01 00:00:00", order="ymdhms") #convert juld->time
  time<-NA
  time<-origin + juld*3600*24
  time <- date(time)
  jd_qc<-NA
  jd_qc<-substr(ncvar_get(nc,"JULD_QC"),1,1)
    
    
  lat <- ncvar_get(nc, "LATITUDE")
  lon <- ncvar_get(nc, "LONGITUDE")
    
  pres <- ncvar_get(nc, "PRES")
    
  chla <- ncvar_get(nc, "CHLA")
  chla_qc <- ncvar_get(nc,"CHLA_QC")
  chla_qctab       <- llply(chla_qc,function(qcstring){
      as.numeric(unlist(strsplit(qcstring,split="")))
    })
  chla_qctab <- do.call(cbind,chla_qctab)
  p <- which(!is.na(chla[1,]))
  chla <- chla[,p]
  chla_qc <- chla_qctab[,p]
    
  bbp <- ncvar_get(nc, "BBP700")
  bbp <- bbp[,p]
    
  chla_adjusted <- ncvar_get(nc, "CHLA_ADJUSTED")
  p <- which(!is.na(chla_adjusted[1,]))
  chla_adjusted <- chla_adjusted[,p]
    
  chla_adjusted_qc <- ncvar_get(nc,"CHLA_ADJUSTED_QC")
  chla_adjusted_qctab  <- llply(chla_adjusted_qc,function(qcstring){
      as.numeric(unlist(strsplit(qcstring,split="")))
    })
  chla_adjusted_qctab <- do.call(cbind,chla_adjusted_qctab)
  chla_adjusted_qc <- chla_adjusted_qctab[,p]
  
  pres <- ncvar_get(nc, 'PRES')
  pres <- pres[,p]
  time <- rep(time[p], length(chla))
  lon <- rep(lon[p], length(chla))
  lat <- rep(lat[p], length(chla))
    
  id <- rep(substr(i, 22,28), length(chla))
    
  nc_df <- data.frame("date" = time, "lat" = lat, "lon"=lon ,"pres" = pres, "chla" = chla, "chla_qc" = chla_qc, "chla_adjusted" = chla_adjusted, "chla_adjusted_qc" = chla_adjusted_qc, "id" = id)

  argo_data <- bind_rows(argo_data, nc_df)
  a <- a+1
  setTxtProgressBar(pb, a)
  nc_close(nc)
}

#quality check

table(is.na(argo_data$lon)) #404 NA for lon
table(is.na(argo_data$lat)) #404 NA for lat
table(is.na(argo_data$chla)) #62 NA for chla
table(is.na(argo_data$chla_adjusted)) #62 NA chla_adjusted
table(is.na(argo_data$date)) #1 NA for date
table(is.na(argo_data$pres)) #62 NA for pres

argo_data_nona <- na.omit(argo_data) # we lose 465 data point

#profiles check

ggplot(argo_data_nona)+
  geom_point(aes(x = chla_adjusted, y = - pres))+
  facet_wrap(.~id)

table(argo_data_nona$id) #some profiles have only 2 or 8 points
