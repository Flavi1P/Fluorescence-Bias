library(tidyverse)

merged <- read_csv("Data/merged_argo")
profiles <- read_csv("Data/argo/first_profiles")
ref <- read_csv("Data/argo/ref.csv")

pigments <- c("fuco", "peri", "hex", "but", "allo", "tchlb", "zea")

merge_easy <- merged %>% select(date, lovbio, depth, lon.x = lon.y, lat.x = lat.y, pigments, chla, chla_qc, chla_adjusted, chla_adjusted_qc) %>% 
  rename_( "lon" = "lon.x", "lat" = "lat.x")

profiles <- left_join(profiles, ref, by = c("id" = "number"))

for (i in unique(profiles$lovbio)){
t <- filter(profiles, lovbio == i)
t2 <- filter(merge_easy, lovbio == i)
prof_merge <- difference_left_join(t, t2, by = c("pres" = "depth"), max_dist = 1) %>% 
  mutate(dist = abs(depth - pres)) %>% group_by(depth) %>% mutate(mindist = min(dist))
}

table(is.na(prof_merge$tchla))


