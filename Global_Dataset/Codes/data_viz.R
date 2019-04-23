library(tidyverse)

pigments <- c("fuco", "peri", "hex", "but", "allo", "tchlb", "zea")

argo <- read_csv("argo/Data/merged_argo") %>%
  select(pigments, tchla, micro, nano, pico, ze, chla_adjusted) %>%
  mutate(campagne = "argo",
         fluo = chla_adjusted * 2,
         ratio = fluo/tchla) %>% 
  select(- chla_adjusted)
biosope <-  read_csv("Biosope/Data/biosope") %>%
  select(pigments, tchla, micro, nano, pico, ze, fluo_urel) %>%
  mutate(campagne = "biosope",
         fluo = fluo_urel,
         ratio = fluo/tchla) %>% 
  select(- fluo_urel)
boussole <- read_csv("Boussole/Data/boussole.csv") %>%
  select(pigments, tchla, micro, nano, pico, ze, fluo) %>%
  mutate(campagne = "boussole",
         ratio = fluo/tchla)


dataset_tall <- bind_rows(argo, biosope, boussole) %>% 
  select(micro, nano, pico, campagne, ratio) %>% 
  group_by(campagne) %>% 
  gather(key = "size_class", value = "frequence", 1:3) %>% 
  ungroup()

ggplot(dataset)+
  geom_violin(aes(x = campagne, y = frequence, fill = frequence))+
  facet_grid(vars(size_class))+
  theme_bw()

ggplot(dataset)+
  geom_boxplot(aes(x = size_class, y = frequence, fill = campagne))+
  scale_fill_brewer(palette = "Set1")+
  ylim(0,1)+ylab("fraction")+
  theme_bw()

ggplot(dataset)+
  geom_boxplot(aes(x = campagne, y = ratio))+
  ylim(0,15)+
  theme_bw()




