profile_numb <- function(depth, direction){
  a <- 1
  profile <- rep(1, length(depth))
  if (direction == "upward"){
  for (i in 2:length(depth)){
  profile[i] <- ifelse(depth[i]>=depth[i-1], a+1, a)
  a <- ifelse(depth[i]>=depth[i-1], a+1, a)
}}
  if (direction == "downward"){
    for (i in 2:length(depth)){
      profile[i] <- ifelse(depth[i]<=depth[i-1], a+1, a)
      a <- ifelse(depth[i]<=depth[i-1], a+1, a)
    }}
  return(profile)
}