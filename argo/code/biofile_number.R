library(tidyverse)

ref <- read_csv("Data/argo/ref.csv")
ref_bis <- read_csv("Data/argo/ref_bis")
ref <- bind_rows(ref, ref_bis)

datalist <- list.files('argo/Data/biofiles')
datalist <- substr(datalist, 3,9)

missing_bio <- which(!ref$number %in% datalist)
ref$lovbio[missing_bio]
