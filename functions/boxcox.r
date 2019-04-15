# Transormation with boxcox :

# Switch the values below 1 to above 1
min.cols <- sapply(zt, function(x) min(x, na.rm=T))
idx <- which(min.cols < 1)
for (ii in idx) {
  zt[[ii]] <- zt[[ii]] + abs(min.cols[ii]) + 1
}
summary(zt)
#-> no more values under 1 (problematic for powerTransform function of car package)


#Find the coefficients for the boxcox transformation

# Calcultate lambda :
lambda <- c()
for (ii in 1:36){
  tic()
  x <- car::powerTransform(zt[[ii]])
  lambda[ii] <- x$lambda
  cat("Columns",ii,"done, lambda is:",x$lambda,". ")
  toc()
}
names(lambda) <- colnames(zt)

# Transform with boxcox :
#no_trans <- c("range", "nb1", "nb2", "nb3", "mode", "min", "fcons", "%area")
for (ii in colnames(zt)){
  zt[[ii]] <- bcPower(zt[[ii]], lambda[[ii]])
}
