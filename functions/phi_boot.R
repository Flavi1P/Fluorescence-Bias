phi_boot <- function(.data, variable = "fluo"){
  temp <- select(.data, "microfluo", "nanofluo", "picofluo", variable, optical_layer)
  coeff_boot <- data.frame("phi" = NA, "se" = NA, "optical_layer" = NA)
  for (i in unique(temp$optical_layer)){
    t <- filter(temp, optical_layer ==  i)
    t <- na.omit(t)
    A <- as.matrix(select(t, microfluo, nanofluo, picofluo))
    B <- as.matrix(select(t, variable))
    coeff_t <- data.frame(coef(nnls(A,B)))
    
    vY <- as.matrix(t[,-5])[,4]
    mX <- as.matrix(t[,-5])[,-4]
    
    vBeta <- as.matrix(coeff_t[,1])
    dSigmaSq <- sum((vY - mX%*%vBeta)^2)/(nrow(mX)-ncol(mX))  # estimate of sigma-squared
    mVarCovar <- dSigmaSq*chol2inv(chol(t(mX)%*%mX))          # variance covariance matrix
    vStdErr <- sqrt(diag(mVarCovar))                          # coeff. est. standard errors
    
    
    coeff_t$se <- vStdErr
    coeff_t$optical_layer <- i
    names(coeff_t) <- c("phi", "se", "optical_layer")
    coeff_boot <- rbind(coeff_boot, coeff_t)
  }
  
  coeff_boot <- na.omit(coeff_boot)
  coeff_boot$size <- c("micro", "nano", "pico")
  return(coeff_boot)
}