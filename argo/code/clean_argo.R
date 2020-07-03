library(tidyverse)
library(vegan)
library(zoo)
library(janitor)
mergefile <- read_csv('Data/argo/first_profiles') %>% select(-temp)
biofile <- read_csv('Data/argo/Biodata')


mergefile$type <- rep('MR', length(mergefile$date))
biofile$type <- rep('BD', length(biofile$date))

data_argo <- bind_rows(biofile, mergefile)

ggplot(data_argo)+
  geom_path(aes(x = chla_adjusted, y = - pres, colour = type))+
  ylim(-300, 0)+
  facet_wrap(.~id, scales = 'free_x')

# No data 6901032 6901521 6901524 6901649 6902739 6902743 6902880 6902742  6902735
# Blurry signal 6901511 6901580 6901689 6902735 6902737 6902739
# Choes BD instead of MR 6901514 6901515 6901525 6901527 6901581 6901770 691658 6901656

data_argo <- data_argo %>% group_by(id, type) %>% 
  mutate(chl_smooth = rollmean(chla_adjusted, 20, fill = c(0,0,0)))

# ggplot(data_argo)+
#   geom_path(aes(x = chl_smooth, y = - pres, colour = type))+
#   ylim(-300, 0)+
#   facet_wrap(.~id, scales = 'free_x')
# Blurry signal 6901580 6902735 6902739

exclude_id <- c(6901580, 6902735, 6902739, 6901032, 6901524, 6901649, 6902739, 6902743, 6902880, 6902742, 6902735, 6901521)
BD_choice <-  c(6901514, 6901515, 6901525, 6901527, 6901581, 6901770, 691658, 6901656)
data_argo_filter <- filter(data_argo, !id %in% exclude_id) %>% filter((id %in% BD_choice & type == 'BD') | (!id %in% BD_choice & type == 'MR'))

test <- data_argo_filter %>% filter(id == 6901574)
unique(test$type)

# ggplot(data_argo_filter)+
#   geom_path(aes(x = chl_smooth, y = - pres, colour = type))+
#   ylim(-200, 0)+
#   facet_wrap(.~id, scales = 'free_x')


argo_afc <- data_argo_filter %>% 
  ungroup() %>% 
  select(id, pres, chl_smooth) %>% 
  filter(pres < 300) %>% 
  mutate(depth = round(pres)) %>%
  group_by(depth, id) %>% 
  summarise_at('chl_smooth', list(mean), na.rm = TRUE)

new_table <- data.frame('depth' = rep(1:200, 47), 'id' = rep(unique(argo_afc$id), each = 200))
new_table <- left_join(new_table, argo_afc)

new_table <- arrange(new_table, id, depth)
for (i in unique(new_table$id)){
  t <- filter(new_table, id == i) 
  for (j in c(1:length(t$depth))) {
    if (is.na(t$chl_smooth[j])) {
      t$chl_smooth[j] <- t$chl_smooth[j+1]
    }
    if (is.na(t$chl_smooth[j])) {
      jj <- j
       while (is.na(t$chl_smooth[j]) & jj %in% c(1:length(t$depth))) {
         t$chl_smooth[j] <- t$chl_smooth[jj]
        jj <- jj+1
       }
     }
      if (is.na(t$chl_smooth[j])){
        t$chl_smooth[j] <- 0
     }  
  t <- arrange(t, depth)  
  new_table[new_table$id == i,] <- t  
  }
}

table(is.na(new_table$chl_smooth))

ggplot(new_table)+
  geom_path(aes(x = chl_smooth, y = - depth))+
  facet_wrap(.~id, scales = 'free_x')
 
argo_afc <- new_table %>%   pivot_wider(names_from = depth, values_from = chl_smooth) %>% 
  clean_names()


afc_chl <- cca(select(argo_afc, x50:x150), scale = TRUE)
#extract scores of eache point
scores <- data.frame(scores(afc_chl, choices = c(1,2,3), display = "site"))
argo_afc <- bind_cols(argo_afc, scores) #addthem to lov df

#plot the all
ggplot(argo_afc)+
  geom_point(aes(x = CA1, y = CA2))

#create tree cluster
distargo <- dist(select(argo_afc, CA1, CA2))
argo_afc$group <- as.factor(cutree(hclust(distargo, method = "ward.D"), k = 2))

ggplot(argo_afc)+
  geom_point(aes(x = CA1, y = CA2, colour = group))

argo_tall <- argo_afc %>% pivot_longer(2:201, names_to = 'depth') %>% 
  mutate(depth2 = as.numeric(gsub('x', '', .$depth))) %>% 
  select(-depth)

ggplot(argo_tall)+
  geom_path(aes(x = value, y = -depth2, colour = group))+
  facet_wrap(.~ id, scale = 'free_x')
#write_csv(data_argo_filter, 'Data/argo/clean_argo')
