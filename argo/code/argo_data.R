library(ncdf4)
library(tidyverse)
library(readxl)
library(lubridate)
library(plyr)


#download.file("ftp://ftp.ifremer.fr/ifremer/argo/argo_bio-profile_index.txt",
#              "argo/Data/index_bio.txt", quiet = FALSE, mode = "w",
#              cacheOK = TRUE) 
#biofiles <- read_csv("argo/Data/index_bio.txt", skip = 7)
profiles <- list.files("Data/argo", full.names = TRUE)
profile_names <- list.files("Data/argo")
profiles <- profiles[grep(pattern = "^[M][R]", profile_names)]
profile_names <- profile_names[grep(pattern = "^[M][R]", profile_names)]
profile_code <- as.numeric(substr(profile_names, 3, 9))

ref <- read_csv("Data/argo/ref.csv")
refbis <- data.frame("number" = c(6902737, 6902739, 6902735, 6902742, 6902743, 6902880), "lovbio" = c("lovbio103c", "lovbio107c", "lovbio100c", "lovapm002a", "lovapm004a", "lovbio111b"))
refbis$lovbio <- as.character(refbis$lovbio)
ref <- bind_rows(ref, refbis)
ref <- filter(ref, number != 6901526)

argo_data <- data.frame("date" = NA, "lat" = NA, "lon" = NA, "pres" = NA, "temp" = NA,
                     "chla" = NA, "chla_qc" = NA, "chla_adjusted" = NA, "chla_adjusted_qc" = NA, "id"= NA)

first_profiles <- paste("Data/argo", "/MR", unique(ref$number), "_001.nc", sep = "")
first_profiles[which(first_profiles == "Data/argo/MR6901521_001.nc")] <- "Data/argo/MR6901521_001D.nc"
first_profiles[which(first_profiles == "Data/argo/MR6901524_001.nc")] <- "Data/argo/MR6901524_001D.nc"

a <- 1

#suite####

pb <- txtProgressBar(min = 1, max = length(first_profiles), style = 3)
for(i in first_profiles){
  nc <- NA
  nc <- nc_open(i)
  variables <- names(nc$var) %>% data.frame %>% filter(. == "CHLA")
  if (nrow(variables)==1){
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
  temp <- ncvar_get(nc, "TEMP")
  temp <- temp[,1]
  
  
  
  chla <- ncvar_get(nc, "CHLA")
  chla_qc <- ncvar_get(nc,"CHLA_QC")
  chla_qctab       <- llply(chla_qc,function(qcstring){
    as.numeric(unlist(strsplit(qcstring,split="")))
  })
  chla_qctab <- do.call(cbind,chla_qctab)
  p <- ifelse(length(chla[1,]) < 5, which(!is.na(chla[1,])), 5)
  chla <- chla[,p]
  chla_qc <- chla_qctab[,p]
  
  bbp <- ncvar_get(nc, "BBP700")
  bbp <- bbp[,p]
  
  chla_adjusted <- ncvar_get(nc, "CHLA_ADJUSTED")
  chla_adjusted <- chla_adjusted[,p]
  
  chla_adjusted_qc <- ncvar_get(nc,"CHLA_ADJUSTED_QC")
  chla_adjusted_qctab  <- llply(chla_adjusted_qc,function(qcstring){
    as.numeric(unlist(strsplit(qcstring,split="")))
  })
  chla_adjusted_qctab <- do.call(cbind,chla_adjusted_qctab)
  chla_adjusted_qc <- chla_adjusted_qctab[,p]
  
  pres <- pres[,p]
  time <- rep(time[p], length(chla))
  lon <- rep(lon[p], length(chla))
  lat <- rep(lat[p], length(chla))
  
  id <- rep(substr(first_profiles[a], 13,19), length(chla))
  
  nc_df <- data.frame("date" = time, "lat" = lat, "lon"=lon ,"pres" = pres, "temp" = temp, "chla" = chla, "chla_qc" = chla_qc, "chla_adjusted" = chla_adjusted, "chla_adjusted_qc" = chla_adjusted_qc, "id" = id, "down" = (rep("non_evaluate", length(pres))))
  nc_df$down <- as.character(nc_df$down)
  nc
  
  if(unique(nc_df$id) %in% refbis$number){
    for(j in c(1:length(na.omit(nc_df$pres)))-1){
  nc_df$down[j] <- ifelse(nc_df$pres[j] < nc_df$pres[j+1], "down", "up")
    }
    nc_df <- filter(nc_df, down == "down")
  }
  
nc_df <- nc_df %>% select(-down)
  argo_data <- bind_rows(argo_data, nc_df)
  }
  a <- a+1
  #setTxtProgressBar(pb, a)
  nc_close(nc)
}

write_csv(argo_data, "Data/argo/first_profiles")
write_csv(refbis, "Scripts/Data/argo/ref_bis")

ggplot(profile_shit)+
  geom_path(aes(x = chla, y = -pres, colour = "chla"))+
  ylim(-45,0)

which(duplicated(merged$tchla))

subref <- subref[grep("lovbio", subref$lovbio),]
ref <- ref[grep("lovbio", ref$lovbio),]
no_date <-"first"

hplc_data <- data.frame("depth" = NA, "lon" = NA, "lat" = NA, "date" = NA, "tchla" = NA, "fuco" = NA, "zea" = NA, "allo" = NA, "peri" = NA, "tchlb" = NA, "hex" = NA, "but" = NA, "id" = NA)

for (i in ref$lovbio){
 hplc_name <- list.files(paste("Scripts/Data/IN_SITU", i, "PIGMENTS", sep = "/"))
 hplc_name <- hplc_name[grep(".xls", hplc_name)]
 temp_pig <- read_excel(paste("Scripts/Data/IN_SITU", i, "PIGMENTS", hplc_name, sep = "/"))
 
 colnames(temp_pig) <- tolower(colnames(temp_pig))
 temp_pig <- data.frame(lapply(temp_pig, function(v) {
   if (is.character(v)) return(tolower(v))
   else return(v)
 }))
 if (any(grep("sampling.date|date$", colnames(temp_pig))) == TRUE){
 pig_new <- data.frame("depth" = temp_pig[,grep("depth", colnames(temp_pig))],
                       "lon" = as.numeric(as.character(temp_pig[,grep("lon", colnames(temp_pig))])),
                       "lat" = as.numeric(as.character(temp_pig[,grep("lat", colnames(temp_pig))])),
                       "date" = as.character(temp_pig[,grep("sampling.date|date$", colnames(temp_pig))]),
                       "tchla" = as.numeric(as.character(temp_pig[,grep("total.chlorophyll.a|tchla$", colnames(temp_pig))])),
                       "fuco" = as.numeric(as.character(temp_pig[,grep("^fucoxanthin",colnames(temp_pig))])),
                       "zea" = as.numeric(as.character(temp_pig[,grep("^zeaxanthin",colnames(temp_pig))])),
                       "allo" = as.numeric(as.character(temp_pig[,grep("^alloxanthin",colnames(temp_pig))])),
                       "peri" = as.numeric(as.character(temp_pig[,grep("^peridinin",colnames(temp_pig))])),
                       "tchlb" = as.numeric(as.character(temp_pig[,grep("total.chlorophyll.b|tchlb$",colnames(temp_pig))])),
                       "hex" = as.numeric(as.character(temp_pig[,grep("hexanoyloxyfucoxanthin$",colnames(temp_pig))])),
                       "but" = as.numeric(as.character(temp_pig[,grep("butanoyloxyfucoxanthin$",colnames(temp_pig))])))
 pig_new$id <- rep(i, length(pig_new$depth))
 hplc_data <- rbind(hplc_data, pig_new)
 }
 if (any(grep("sampling.date|date$", colnames(temp_pig))) == FALSE){
  no_date <- append(no_date, i)
}
 
}

hplc_data <- hplc_data[-1,]
hplc_data[is.na(hplc_data)] <- 0



hplc_positiv <- hplc_data
hplc_positiv[hplc_positiv < 0] <- NA
hplc_data <- bind_cols(select(hplc_data, -pigments, - tchla),select(hplc_positiv, pigments, tchla))

hplc_mduf <- read_excel("Scripts/Data/IN_SITU/MOBYDICK_pigments_130718_.xlsx", 
                        col_types = c("numeric", "text", "numeric", 
                                      "numeric", "numeric", "date", "numeric", 
                                      "numeric", "numeric", "numeric", "numeric", 
                                      "numeric", "numeric", "numeric", 
                                      "numeric", "numeric", "numeric", 
                                      "numeric", "numeric", "numeric", "numeric", 
                                      "numeric", "numeric", "numeric", "numeric", 
                                      "numeric", "numeric", "numeric", "numeric", 
                                      "numeric", "numeric", "numeric", "numeric", 
                                      "numeric", "numeric", "numeric", "numeric", 
                                      "numeric", "numeric", "numeric", "numeric", 
                                      "numeric", "numeric", "numeric", "numeric", 
                                      "numeric", "numeric", "numeric", "numeric", 
                                      "numeric", "numeric", "numeric", "numeric", 
                                      "numeric", "numeric", "numeric", "numeric", 
                                      "numeric", "numeric", "numeric", 
                                      "numeric", "numeric", "numeric", "numeric", 
                                      "numeric", "numeric"))


hplc_mduf <- select(hplc_mduf, "Depth_(m)", Longitude, Latitude, "Sampling_date_(UTC)", "Ship", Fucoxanthin, Peridinin, "19-Hexanoyloxyfucoxanthin",
                    "19-Butanoyloxyfucoxanthin", Alloxanthin, "Total_Chlb", Zeaxanthin, "Total_Chlorophyll_a")
names(hplc_mduf) <- c("depth", "lon", "lat", "date", "id", "fuco", "peri", "hex", "but", "allo", "tchlb", "zea", "tchla")
hplc_mduf[is.na(hplc_mduf)] <- 0
#write_csv(hplc_mduf, "Scripts/Data/argo/hplc_mduf")

hplc_tak <- read_excel("Scripts/Data/IN_SITU/GreenEdge-Amundsen-pigments_flotteurs-300117.xlsx", 
                       col_types = c("numeric", "text", "text", 
                                     "text", "numeric", "numeric", "date", 
                                     "numeric", "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", "numeric"))

hplc_tak <- select(hplc_tak, "Depth (m)", Longitude, Latitude, "Sampling date (UTC)", "Leg", Fucoxanthin, Peridinin, "19'-Hexanoyloxyfucoxanthin",
                   "19'-Butanoyloxyfucoxanthin", Alloxanthin, "Chlorophyll b", Zeaxanthin, "Total Chlorophyll a")
names(hplc_tak) <- c("depth", "lon", "lat", "date", "id", "fuco", "peri", "hex", "but", "allo", "tchlb", "zea", "tchla")
#write_csv(hplc_tak, "Scripts/Data/hplc_tak")


hplc_soclim <- read_excel("Scripts/Data/IN_SITU/SOCLIM2016_241017.xlsx", 
                          col_types = c("numeric", "numeric", "numeric", 
                                        "numeric", "date", "numeric", "text", 
                                        "text", "text", "text", "numeric", 
                                        "numeric", "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric"))

hplc_soclim <- select(hplc_soclim, "Depth (m)", Longitude, Latitude, "Sampling date (UTC)", "Station", Fucoxanthin, Peridinin, "19'-Hexanoyloxyfucoxanthin",
                   "19'-Butanoyloxyfucoxanthin", Alloxanthin, "Chlorophyll b", Zeaxanthin, "Total Chlorophyll a")

names(hplc_soclim) <- c("depth", "lon", "lat", "date", "id", "fuco", "peri", "hex", "but", "allo", "tchlb", "zea", "tchla")

#write_csv(hplc_soclim, "Scripts/Data/hplc_soclim")


ggplot(hplc_data)+
  geom_density(aes(x = as.numeric(date)))

origin<-as.POSIXct("1950-01-01 00:00:00", order="ymdhms") #convert juld->time

date_test <- substr(hplc_data$date, 1, 10)
date_test<- as.POSIXct("2015-05-12")
date_vec <- c()

for (i in c(1:length(hplc_data$date))){
  if (nchar(hplc_data[i, "date"]) == 10){
    temp_date <- gsub("/", "-", hplc_data[i, "date"])
    temp_date <- date(temp_date)
  }
  if (nchar(hplc_data[i, "date"]) < 10){
    temp_date <- as.numeric(hplc_data[i, "date"])
    temp_date <- origin + temp_date
  }
  date_vec <- append(date_vec, temp_date)
}


filter(hplc_data, new_date == "2013-05-10")
unique(date_vec)

origin<-as.POSIXct("1899-12-30 00:00:00", order="ymdhms")
hplc_data$new_date <- date_vec



hplc_data <- left_join(hplc_data, ref, by = c("id" = "lovbio"))





map_vec <- read_csv("Scripts/Data/map_vec")
ggplot(hplc_data)+
  geom_point(aes(x = lon, y = lat, colour = id))+
  geom_polygon(aes(x = long, y = lat, group = group), data = map_vec)+
  coord_quickmap()

write_csv(hplc_data, path = "Scripts/Data/argo/hplc_argo")
ncfile <- nc
Var <- "CHLA"
ExtractVar <- function(Var,FloatInfo){
  with(FloatInfo,{
    # This function should return a dataframe for variable with value, qc, iprofile and ilevel
    lvar             <- ncvar_get(ncfile,Var)
    lvar_qc          <- ncvar_get(ncfile,paste0(Var,"_QC"))
    lvar_qctab       <- llply(lvar_qc,function(qcstring){
      as.numeric(unlist(strsplit(qcstring,split="")))
    })
    lvar_qctab <- do.call(cbind,lvar_qctab)
    
    lvar_dir       <- ncvar_get(ncfile,"DIRECTION")
    lvar_direction <- llply(lvar_dir,function(dirstring){
      strsplit(dirstring,split="")
    })
    lvar_direction <- unlist(lvar_direction)
    # making dataframes, removing the NANs  
    alevels <- 1:6
    d <- ldply(as.list(1:6),function(iprof){
      indexes <- !(is.na(lvar[,iprof])|is.na(lvar[,iprof]))
      if(sum(indexes) == 0){
        return (data.frame())
      }
      
      data.frame(value    = lvar[indexes,iprof],
                 qc       = as.integer(lvar_qctab[indexes,iprof]),
                 alevel   = alevels[indexes],
                 depth    = pres[indexes,iprof],
                 dir      = lvar_direction[iprof],
                 aprofile = iprof,
                 variable = Var)
    })
    
    d$juld <- juld[d$aprofile]
    d$lon  <- lon[d$aprofile]
    d$lat  <- lat[d$aprofile]
    
    return(d=d)
  })
} 
    
ExtractVar("CHLA", nc)


p <- which(!is.na(chla[1,]))


profiles_bis <- filter(argo_data, id %in% refbis$number)

