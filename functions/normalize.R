# Transormation with boxcox :

normalize <- function(.data){
  
  # Switch the values below 1 to above 1
  d <- .data

  for(i in 1:ncol(d)){
    var_trans <- jtrans::jtrans(d[[i]])
    d[,i] <- var_trans$transformed
  }

  return(d)
}
