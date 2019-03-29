phi_lm<- function(.data, variable = "corrected"){
  temp <- select(.data, "microfluo", "nanofluo", "picofluo", variable, optical_layer)
  names(temp) <- c("microfluo", "nanofluo", "picofluo", "corrected", "optical_layer")
  coeff_lm <- data.frame("phi" = NA, "se" = NA, "optical_layer" = NA)
  for (i in unique(temp$optical_layer)){
    model <- lm(corrected~microfluo + nanofluo + picofluo + 0, data = filter(temp, optical_layer == i))
    coeff_t <- data.frame(coef(summary(model)))
    coeff_t <- coeff_t[,c(1,2)]
    coeff_t$optical_layer <- i
    names(coeff_t) <- c("phi","se", "optical_layer")
    coeff_lm <- rbind(coeff_lm, coeff_t)
  }
  
  coeff_lm <- na.omit(coeff_lm)
  coeff_lm$size <- c("micro", "nano", "pico")
  return(coeff_lm)
}