library(tidyverse)
library(Metrics)
library(vegan)
library(nnls)
library(sf)
library(FactoMineR)
source("functions/phi_lm.R")
source("functions/outliers.R")
source("functions/phi_boot.R")
path = "Data/Longhurst"
pigments <- c("fuco", "peri", "hex", "but", "allo", "tchlb", "zea")
NAT_IRS_list <- c("lovbio059c", "lovbio045b", "lovbio024c", "lovbio044b", "lovbio031c", "lovbio027b", "lovbio040b", "lovbio026c")

map_vec <- read_csv("Data/map_vec")
merged_argo <- read_csv("Data/merged_argo")
merged_argo <- filter(merged_argo, lovbio != "lovbio067c" & lovbio != "lovbio083d" & lovbio != "lovbio085d" & lovbio != "lovbio090d") #filter profile with only 1 match
merged_argo <- filter(merged_argo, !(lovbio %in% NAT_IRS_list))

argo_data <- merged_argo %>%
  select(date:zze) %>% 
  rename("lat_float" = "lat.x", "lon_float" = "lon.x", "lat_hplc" = "lat.y", "lon_hplc" = "lon.y")

write_csv(argo_data, "argo/Data/merged_data_argo_hplc")
