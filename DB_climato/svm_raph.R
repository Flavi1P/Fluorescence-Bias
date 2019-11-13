library(maps)
library(e1071)
library(tidyverse)
library(dplyr)
library(hydroGOF)
library(ggplot2)
library(janitor)
library(readxl)

#set the path to directories for data and plots
dir_data <- "DB_climato/Data/"
dir_plot <- "DB_climato"
# dir_lush <- "/home/sauzede/Documents/LUSH/FLAVIEN_NN/"

#load database Flavien
database <- read.table(paste(dir_data,"database_final",sep=""),h=T,sep=",")

#remove ratio = NA and PAR = NA
database <- database[!is.na(database$ratio),]
database <- database[!is.na(database$par),]

#remove outliers ratio --> Flavien check pourquoi ces ratios sont là?
#database <- filter(database, ratio<10 & ratio>1)

#Transform ratio in lo because lognormal distribution
database$logratio <- log10(database$ratio)

#add a N to each porfile with unique lon/lat/date
LON_LAT_DATE <- unique(database[,1:5])
dim(LON_LAT_DATE) #1728 profiles
LON_LAT_DATE$N <- 1:length(LON_LAT_DATE[,1])
#add this N profile to the database
database <- merge(database,LON_LAT_DATE)


#function to transform in radians the longitude and doy (circular variables)
rad.lon <- function(lon){
  radians <- (lon*pi)/180
  radians
}
rad.doy <- function(dat){
  radians <- (dat*pi)/182.5
  radians
}

#use sin and cos of lon and doy radians as inputs (see Sauzede et al. 2015 for details)
database$sin_lon <- sin(rad.lon(round(database$lon,2)))
database$cos_lon <- cos(rad.lon(round(database$lon,2)))

database$sin_doy <- sin(rad.doy(round(database$doy)))
database$cos_doy <- cos(rad.doy(round(database$doy)))

#seperate randomly the database in 80% of the data used for the training and 20% for the validation
set.seed(123)
iii_data_train <- sample(1:length(LON_LAT_DATE[,1]), round(0.8 * length(LON_LAT_DATE[,1])))
iii_data_train <- sort(iii_data_train)

length(iii_data_train)
# 1382
length(LON_LAT_DATE[,1])
# 1728

# 1728-1382 #346 stations pour validation

# #for exemple check the dfference between the distribution of rrs and log(rrs)
# hist(database$rrs412)
# hist(log(database$rrs412))

#add to the database the logtransformation of rrs, mld and depth
database <- mutate(.data = database, log_rrs667 = log(rrs667),
                   log_rrs555 = log(rrs555), log_rrs488 = log(rrs488),
                   log_rrs443 = log(rrs443), log_rrs412 = log(rrs412),
                   log_mld = log(mld), log_depth = log(depth))

#plot the histograms of each variable from the database
# pdf(paste(dir_plot,"Hist_data.pdf"))
# par(mfrow=c(2,2))
# for(i in 14:dim(database)[2]){
#   hist(database[,i],main=paste("Histogramme",colnames(database)[i]))
# }
# graphics.off()

# database2 <- database[,c("N","depth","sin_lon","cos_lon","sin_doy","cos_doy","lat","par", "log_rrs667", "log_rrs555", 
#                          "log_rrs488", "log_rrs443", "log_rrs412", "log_mld", "logratio")]


#create the database with the inputs that we want to have in our model (save the N for the separation between training and validation profiles)
database2 <- database[,c("N","depth","sin_lon","cos_lon","sin_doy","cos_doy","lat","par", "log_rrs667", "log_rrs555", 

                                                  "log_rrs488", "log_rrs443", "log_rrs412", "log_mld", "logratio")]

database2 <- database2[!is.na(database2$log_rrs667),]


#seperate the database for training and validation and remove the N input
DATABASE_TRAIN <- database2[database2$N %in% iii_data_train,-1]
DATABASE_VALID <- database2[!database2$N %in% iii_data_train,-1]

#plot the geographical distribution of training and validation profiles
# png(paste(dir_plot,"Database.png",sep=""),res=300, width=400*7,height=400*4)
# par(las=1)
# plot(database$lon[database$N %in% iii_data_train],database$lat[database$N %in% iii_data_train], pch=19, col="turquoise4", xlab="Longitude", ylab= "Latitude", 
#      main = paste(length(iii_data_train),"training stations and", length(LON_LAT_DATE[,1])-length(iii_data_train), "validation stations",sep=" "))
# points(database$lon[!database$N %in% iii_data_train],database$lat[!database$N %in% iii_data_train], pch=19, col="violetred3")
# map(add=T,col="black", fill=T)
# legend("bottomright", col= c("turquoise4","violetred3"), legend = c("training","validation"), pch=19, box.lty=0)
# graphics.off()

#funcion to transform the pressure input as in Sauzede et al. 2017 and Bittig et al. 2018 (CANYON method)
P_trans <- function(P){
  P_trans <- P/2e4 + 1 / (1+exp(-P/300))^3
  return(P_trans)
}

#transform the depth input (we should transform before the depth in pressure with oce library... J'ai pas eu le temps, SORRY)
DATABASE_TRAIN$P <- P_trans(DATABASE_TRAIN$depth)
DATABASE_VALID$P <- P_trans(DATABASE_VALID$depth)

DATABASE_TRAIN <- select(DATABASE_TRAIN, - depth)
DATABASE_VALID <- select(DATABASE_VALID, - depth)

#Training subsets
#input training subset
x_train <- subset(DATABASE_TRAIN, select = -logratio)
#output training subset
y_train <- DATABASE_TRAIN$logratio

#Validation subsets
#input validation subset
x_valid <- subset(DATABASE_VALID, select = -logratio)
#output validation subset
y_valid <- DATABASE_VALID$logratio

#compute mean and standard deviation to centre and reduce the inputs and outputs using the TRAINING data only
MEAN_DATA <- c(apply(x_train,2, mean), mean(y_train)) 
SD_DATA <- c(apply(x_train,2, sd), mean(y_train)) 

#centre and reduce each input
for(c in 1:dim(x_train)[2]){
  x_train[,c] <- 2/3 * (x_train[,c]-MEAN_DATA[c])/SD_DATA[c]
  x_valid[,c] <- 2/3 * (x_valid[,c]-MEAN_DATA[c])/SD_DATA[c]
}
#centre and reduce each output
y_train <- 2/3 * (y_train-MEAN_DATA[14])/SD_DATA[14]
y_valid <- 2/3 * (y_valid-MEAN_DATA[14])/SD_DATA[14]

#x_train$ratio <- y_train
#select best parameters
#tuned_parameters <- tune.svm(ratio~., data = x_train, gamma = 10^(-5:-1), cost = 10^(-3:1)) # this take a while
#summary(tuned_parameters) best parameters gamma = 0.1 & cost = 1, but if we take default ettings we have better rmse and R²

#x_train <- select(x_train, - ratio)

#Train the model of support vector machine
model <- svm(y_train~., data= x_train, gamma = 0.001, cost = 1)


#predict the validation outputs ftom the validation input subset
pred_valid <- predict(model, x_valid)

#decentre annd denormalize the data and delogtransform the ratio for each the estimation and the observation
Estimated_ratio <- 10.^(1.5*pred_valid*SD_DATA[14]+MEAN_DATA[14])
Obs_ratio <- 10.^(1.5*y_valid*SD_DATA[14]+MEAN_DATA[14])

#compute statistics
RMSE <- rmse(sim = Estimated_ratio, obs = Obs_ratio)
LM2 <- lm(log(Estimated_ratio)~log(Obs_ratio))
SM2 <- summary(LM2)
slope2 <- SM2$coefficients[2,1]
intercept2 <- SM2$coefficients[1,1]
R_squared2 <- SM2$adj.r.squared

#plot the scatterplot between the obs and pred values + statistics
# png(paste(dir_plot,"SVM_Ratio_P_trans.png",sep=""),res=300,height=300*8,width=300*8)
# plot(Estimated_ratio,Obs_ratio,log="xy",xlim=range(Obs_ratio,Estimated_ratio),ylim=range(Obs_ratio,Estimated_ratio),
#      xlab = "Ratio from SVM",ylab = "Ratio from pigments")
# abline(a=0,b=1,col="red")
# mtext(paste("RMSE =", round(RMSE,2), ",", "R²", "=", round(R_squared2,2), ", y = a *", round(slope2,2), "+", round(intercept2,2), sep=" "))
# graphics.off()

#ggplot scatterplot to have the colour for depth --> to improve!
DATABASE_VALID2 <- cbind(DATABASE_VALID,Estimated_ratio,Obs_ratio)
ggplot(data = DATABASE_VALID2) +
  geom_point(aes(x = log(Estimated_ratio), y = log(Obs_ratio)))+
  geom_line(aes(x = log(Estimated_ratio), y = log(Estimated_ratio)), col = "Red")+
  theme_minimal()


# predict the ratio on the argo/HPLC database ####

argo <- read_csv("DB_climato/Data/argo_climato")

spectre <- read_excel("Biosope/Data/Spectres_annick.xlsx")
spectre <- clean_names(spectre)

#select the absorbtion at 440 and 470nm
spectre440<- filter(spectre, lambda == 440)
spectre470 <- filter(spectre, lambda == 470)

#create columns that correspond to the total photosynthetic absorbance and non photosynthetic absorbance at 440 and 470. Create also a ratio between the two photosynthetic absorbtion
argo <- argo %>% mutate(photo_440 = peri * spectre440$peri + but * spectre440$x19_bf + hex * spectre440$x19_hf + fuco * spectre440$fuco + allo * spectre440$allox + tchla * spectre440$chl_a,
                        protect_440 = zea * spectre440$zea,
                        photo_470 = peri * spectre470$peri + but * spectre470$x19_bf + hex * spectre470$x19_hf + fuco * spectre470$fuco + allo * spectre470$allox + tchla * spectre470$chl_a,
                        protect_470 = zea * spectre470$zea, 
                        ratio_abs = photo_440/photo_470)

Obs_ratio <- argo$ratio_abs

argo$sin_lon <- sin(rad.lon(round(argo$lon,2)))
argo$cos_lon <- cos(rad.lon(round(argo$lon,2)))

argo$sin_doy <- sin(rad.doy(round(argo$doy)))
argo$cos_doy <- cos(rad.doy(round(argo$doy)))


#add to the database the logtransformation of rrs, mld and depth
argo_predict <- argo %>%  mutate( log_rrs667 = log(rrs667),
                   log_rrs555 = log(rrs555), log_rrs488 = log(rrs488),
                   log_rrs443 = log(rrs443), log_rrs412 = log(rrs412),
                   log_mld = log(mld), log_depth = log(depth),
                   P = P_trans(depth)) %>% 
  select(colnames(x_train))


#centre and reduce each input
for(c in 1:dim(argo_predict)[2]){
  argo_predict[,c] <- 2/3 * (argo_predict[,c]-MEAN_DATA[c])/SD_DATA[c]
}
#centre and reduce each output
# Obs_ratio <- 2/3 * (Obs_ratio-MEAN_DATA[14])/SD_DATA[14]


pred_argo <- predict(model, argo_predict)

#decentre annd denormalize the data and delogtransform the ratio for each the estimation and the observation
Estimated_ratio <- 10.^(1.5*pred_argo*SD_DATA[14]+MEAN_DATA[14])




RMSE <- rmse(sim = Estimated_ratio, obs = Obs_ratio)
LM2 <- lm(log(Estimated_ratio)~log(Obs_ratio))
SM2 <- summary(LM2)
slope2 <- SM2$coefficients[2,1]
intercept2 <- SM2$coefficients[1,1]
R_squared2 <- SM2$adj.r.squared


argo_predict <- cbind(argo_predict, Estimated_ratio,Obs_ratio)
ggplot(data = argo_predict) +
  geom_point(aes(x = log(Estimated_ratio), y = log(Obs_ratio)))+
  geom_line(aes(x = log(Estimated_ratio), y = log(Estimated_ratio)), col = "Red")+
  theme_minimal()

#Correction of the fluo ####

argo$ratio_estimated <- Estimated_ratio

argo <- argo %>% mutate(fitted_math = 13 * ratio_estimated^-1.83) #recompute the estimated values

argo$fluo <- argo$chla_adjusted * 2
#correct the fluorescence signal from there
argo <- argo %>% mutate(corrected_fluo = fluo/fitted_math)

#plot the corrected fluorescence

ggplot(argo)+
  geom_point(aes(x = log(tchla), y = log(fluo)))+
  geom_point(aes(x = log(tchla), y = log(corrected_fluo)), colour = "Red")+
  geom_line(aes(x = log(tchla) ,y = log(tchla)))

rmse(argo$corrected_fluo, argo$tchla)
rmse(argo$chla_adjusted, argo$tchla)


model_chla_adjusted <- lm(argo$chla_adjusted ~ argo$tchla)
model_fluo_corrected <- lm(argo$corrected_fluo~ argo$tchla)

summary(model_chla_adjusted)
summary(model_fluo_corrected)


ggplot(argo)+
  geom_point(aes(x = tchla, y = chla_adjusted))+
  geom_point(aes(x = tchla, y = corrected_fluo), colour = "Red")+
  geom_line(aes(x = tchla ,y = tchla))

#write_csv(argo, "DB_climato/argo_fluo_corrected")

#create and save DATABASE TRAIN and VALID for NN training on lush####
# 
# P_trans <- function(P){
#   P_trans <- P/2e4 + 1 / (1+exp(-P/300))^3
#   return(P_trans)
# }
# 
# database2 <- database[,c("N","depth","sin_doy","cos_doy","sin_lon","cos_lon","lat","par", "rrs667", "rrs555", 
#                          "rrs488", "rrs443", "rrs412", "mld", "logratio")]
# 
# database2$depth <- P_trans(database2$depth)
# 
# DATABASE_TRAIN <- database2[database2$N %in% iii_data_train,-1]
# DATABASE_VALID <- database2[!database2$N %in% iii_data_train,-1]
# 
# MEAN_DATA <- apply(DATABASE_TRAIN,2, mean) 
# SD_DATA <- apply(DATABASE_TRAIN,2, sd)
# 
# for(c in 1:dim(DATABASE_TRAIN)[2]){
#   DATABASE_TRAIN[,c] <- 2/3 * (DATABASE_TRAIN[,c]-MEAN_DATA[c])/SD_DATA[c]
#   DATABASE_VALID[,c] <- 2/3 * (DATABASE_VALID[,c]-MEAN_DATA[c])/SD_DATA[c]
# }
# 
# ne <- 13
# 
# write.table(MEAN_DATA[ne+1],paste(dir_lush,"/MOY_LOGRATIO_NN.dat",sep=""),row.names=F,col.names=F)
# write.table(SD_DATA[ne+1],paste(dir_lush,"/SD_LOGRATIO_NN.dat",sep=""),row.names=F,col.names=F)
# 
# write.table(DATABASE_TRAIN,file=paste(dir_lush,"/DATABASE_TRAIN.dat",sep=""),row.names=F,col.names = F)
# write.table(DATABASE_VALID,file=paste(dir_lush,"/DATABASE_VALID.dat",sep=""),row.names=F,col.names = F)
# 
# dim(DATABASE_TRAIN)
