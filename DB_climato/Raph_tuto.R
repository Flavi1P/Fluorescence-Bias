library(ncdf4)
library(maps)
library(fields)

nc <- nc_open("Data/climato_mld/mld_DR003_c1m_reg2.0.nc")

mld <- ncvar_get(nc,nc$var$mld)

time <- ncvar_get(nc,nc$dim$time)
lat <- ncvar_get(nc,nc$dim$lat)
lon <- ncvar_get(nc,nc$dim$lon)

#MLD for june
dim(mld)
length(lon)
length(lat)
length(time)

mld_june <- mld[,,6]

image.plot(lon,lat,mld_june,zlim=c(0,500))
map(add=T,fill=T,col="black")

nc_close(nc)

nc <- nc_open("Data/climato_oc/A20021822019212.L3m_MC_RRS_Rrs_412_9km.nc")

rrs_412 <- ncvar_get(nc,nc$var$Rrs_412)
lat <- ncvar_get(nc,nc$dim$lat)
lon <- ncvar_get(nc,nc$dim$lon)

rrs_412 <- rrs_412[,order(lat)]
lat <- lat[order(lat)]
image.plot(lon,lat,rrs_412)
map(add=T,fill=T,col="black")

