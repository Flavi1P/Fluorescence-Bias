phi_stat <- function(.data, n , repetition = 100, weighted = TRUE, method = "nnls", variable = "fluo"){
  temp <- select(.data, "microfluo", "nanofluo", "picofluo", variable, optical_layer, id)
  names(temp) <- c("microfluo", "nanofluo", "picofluo", "corrected", "optical_layer", "id")
  temp$weight2 <- NA
  if (method == "nnls"){
  if ( weighted == TRUE){
    for (i in unique(temp$id)){
      for (j in temp$optical_layer){
        temp$weight2[temp$id == i & temp$optical_layer == j] <- 1/(as.numeric(table(temp$id[temp$id == i & temp$optical_layer == j])) * as.numeric(length(temp$id)))
      }
    }
    
    pb <- txtProgressBar(min = 0, max = repetition, style = 3)
    phi_weighted <- data.frame("phi" = NA, "optical_layer" = NA, "size" = NA, "boot"  = NA)
    for (i in c(1:repetition)){
      t <- sample_n(temp, n, weight = weight2)
      phi_t <- phi_boot(t)
      phi_t$boot <- i
      phi_weighted <- rbind(phi_weighted, phi_t)
      setTxtProgressBar(pb, i)
    }
    
    phi_weighted <- na.omit(phi_weighted)
    phi_stat_weighted <- phi_weighted %>% group_by(size, optical_layer) %>% summarise("phi_mean" = mean(phi), "phi_sd" = sd(phi))
    
    phi_stat_weighted$phi_sd <- ifelse(phi_stat_weighted$phi_sd>phi_stat_weighted$phi_mean, phi_stat_weighted$phi_mean, phi_stat_weighted$phi_sd)
    
    return(phi_stat_weighted)
    }
  
  if (weighted == FALSE){

    pb <- txtProgressBar(min = 0, max = repetition, style = 3)
    phi_raw <- data.frame("phi" = NA, "optical_layer" = NA, "size" = NA, "boot"  = NA)
    for (i in c(1:repetition)){
      t <- sample_n(temp, n)
      phi_t <- phi_boot(t)
      phi_t$boot <- i
      phi_raw <- rbind(phi_raw, phi_t)
      setTxtProgressBar(pb, i)
    }
    
    phi_raw <- na.omit(phi_raw)
    phi_stat_raw <- phi_raw %>% group_by(size, optical_layer) %>% summarise("phi_mean" = mean(phi), "phi_sd" = sd(phi))
    
    phi_stat_raw$phi_sd <- ifelse(phi_stat_raw$phi_sd>phi_stat_raw$phi_mean, phi_stat_raw$phi_mean, phi_stat_raw$phi_sd)
    
    return(phi_stat_raw)
    
  }
  
  }
  if (method == "lm"){
    if (weighted == TRUE){
      for (i in unique(temp$id)){
        for (j in temp$optical_layer){
          temp$weight2[temp$id == i & temp$optical_layer == j] <- 1/(as.numeric(table(temp$id[temp$id == i & temp$optical_layer == j])) * as.numeric(length(temp$id)))
        }
      }
      pb <- txtProgressBar(min = 0, max = repetition, style = 3)
      phi_raw <- data.frame("phi" = NA, "se" = NA, "optical_layer" = NA, "size" = NA, "boot"  = NA)
      for (i in c(1:repetition)){
        t <- sample_n(temp, n, weight = weight2)
        phi_t <- phi_lm(t)
        phi_t$boot <- i
        phi_raw <- rbind(phi_raw, phi_t)
        setTxtProgressBar(pb, i)
      }
      phi_raw <- na.omit(phi_raw)
      phi_lm <- phi_raw %>% group_by(size, optical_layer) %>% summarise("phi_mean" = mean(phi), "phi_sd" = mean(se))
      
      return(phi_lm)
    }
  if (weighted == FALSE){
    pb <- txtProgressBar(min = 0, max = repetition, style = 3)
    phi_raw <- data.frame("phi" = NA, "se" = NA, "optical_layer" = NA, "size" = NA, "boot"  = NA)
    for (i in c(1:repetition)){
      t <- sample_n(temp, n, weight = weight2)
      phi_t <- phi_lm(t)
      phi_t$boot <- i
      phi_raw <- rbind(phi_raw, phi_t)
      setTxtProgressBar(pb, i)
      }
    phi_raw <- na.omit(phi_raw)
    phi_lm <- phi_raw %>% group_by(size, optical_layer) %>% summarise("phi_mean" = mean(phi), "phi_sd" = mean(se))
    
    return(phi_lm)
    
    }
  }
}