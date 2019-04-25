library(tidyverse)
library(FactoMineR)

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

ggplot(dataset_tall)+
  geom_boxplot(aes(x = size_class, y = frequence, fill = campagne))+
  scale_fill_brewer(palette = "Set1")+
  ylim(0,1)+ylab("proportion de la classe de taille")+ xlab("classe de taille")+
  theme_bw(base_size = 14)

ggsave("Global_Dataset/Plots/size_class_boxplot.png", scale = 2)

ggplot(dataset_tall)+
  geom_violin(aes(x = campagne, y = ratio, fill = campagne))+
  scale_fill_brewer(palette = "Set1")+
  guides(fill = FALSE)+
  ylim(0,15)+
  theme_bw(base_size = 14)
ggsave("Global_Dataset/Plots/ratio_violin.png", scale = 2)

dataset <- bind_rows(argo, biosope, boussole) %>%
  mutate(pigsum = rowSums(.[, 1:7])) %>%
  filter(pigsum != 0)


afc <- CA(select(dataset, pigments), row.sup = as.numeric(rownames(dataset[dataset$campagne != "argo",])))

scores <- data.frame(afc$row$coord)
scores_sup <- data.frame(afc$row.sup$coord)
scores <- bind_rows(scores, scores_sup)

dataset <- bind_cols(dataset,scores)

pigscore <- data.frame(afc$col$coord)

ggplot()+
  geom_point(aes(x = Dim.1, y = Dim.2, colour = campagne), size = 0.9, alpha = 0.8, data = dataset)+
  geom_segment(aes(x = 0, xend = Dim.1, y = 0, yend = Dim.2), data = pigscore)+
  geom_text(aes(x = Dim.1, y = Dim.2, label = rownames(pigscore)), data = pigscore, size = 5)+
  scale_color_brewer(palette = "Set1")+
  ylab("CA2 (23%)") + xlab("CA1 (43%)")+
  coord_equal()+
  theme_bw(base_size = 14)


ggsave("Global_Dataset/Plots/CA_comp.png", scale = 2)
biosope <- filter(dataset, campagne == "biosope") %>% 
  mutate(coord = paste(round(Dim.1,1), round(Dim.2,1), sep = ";"),
         dup = duplicated(coord)) %>% 
  filter(dup == "FALSE")

ggplot()+
  geom_point(aes(x = Dim.1, y = Dim.2, colour = campagne), data = dataset)+
  geom_point(aes(x = Dim.1, y = Dim.2), data = biosope, colour = "blue")+
  geom_segment(aes(x = 0, xend = Dim.1, y = 0, yend = Dim.2), data = pigscore)+
  geom_text(aes(x = Dim.1, y = Dim.2, label = rownames(pigscore)), data = pigscore)+
  scale_color_brewer(palette = "Set1")+
  theme_bw()
