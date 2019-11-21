library(tidyverse)
library(lubridate)

ref <- read_csv("Data/argo/ref.csv")
ref_bis <- read_csv("Data/argo/ref_bis")

first_profiles <- read_csv("Data/argo/first_profiles")
first_profiles <- first_profiles[-1,]
first_profiles$date <- date(first_profiles$date)

first_profiles$number <- as.numeric(first_profiles$id)
first_profiles <- left_join(first_profiles, ref)
first_profiles <- left_join(first_profiles, ref_bis)


ggplot(first_profiles)+
  geom_path(aes(x = chla_adjusted, y = - pres))+
  facet_wrap(.~id, scales = "free_x")+
  ylim(-500,0)

