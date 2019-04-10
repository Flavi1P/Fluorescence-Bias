outliers <- function(.data, X = .data$fluo){
  mod <- lm(X~., data = .data)
  
  cooksd <- cooks.distance(mod)
  
  influential <- as.numeric(names(cooksd)[(cooksd > 4*mean(cooksd, na.rm=T))]) 
  return(influential)
}
