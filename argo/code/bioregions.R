library(tidyverse)
library(sf)
library(readxl)
library(vegan)
library(ggrepel)
library(gridExtra)
library(grid)
path = "Data/Longhurst"
argo <- read_csv("Data/merged_argo")
map_vec <- read_csv("Data/map_vec")
pigments <- c("fuco", "peri", "hex", "but", "allo", "tchlb", "zea")

longhurst_sf <- read_sf(dsn = path.expand(path), quiet = TRUE)

names(longhurst_sf) <- c("code", "region", "geometry")

#longhurst_sf %>% ggplot() + geom_sf(aes(fill = code))

pnts_sf <- do.call("st_sfc",c(lapply(1:nrow(argo),
                                     function(i) {st_point(as.numeric(argo[i,c("lon.x", "lat.x") ]))}), list("crs" = 4326))) 
pnts_trans <- st_transform(pnts_sf, 4326)
longhurst_trans <- st_transform(longhurst_sf, 4326)  
argo$code <- apply(st_intersects(longhurst_trans, pnts_trans, sparse = FALSE), 2, 
                     function(col) { 
                       longhurst_trans[which(col), ]$code
                     })

ggplot(argo)+
  geom_point(aes(x = lon.x, y = lat.x, colour = code), size = 3)+
  geom_polygon(aes(x = long, y = lat, group = group), data = map_vec)+
  coord_quickmap()

argo <- filter(argo, optical_layer < 4)

nbr_match_region<- count(argo, "code")
nbr_float_region <- argo %>% filter(duplicated(lovbio) == FALSE) %>% count("code")
resume_region <- bind_cols(nbr_float_region, nbr_match_region) %>% select(1,2,4)

names(resume_region) <- c("code", "nbr_of_float", "nbr_of_match")


argo <- argo %>% mutate(fluo = chla_adjusted * 2,
                        ratio = fluo/tchla)

region_argo <- argo  %>% group_by(code) %>% summarise_at(vars(ratio), c(mean, sd), na.rm = TRUE) %>% ungroup()
names(region_argo) <- c("code", "mean", "sd")

region_argo <- left_join(region_argo, resume_region)

region_argo$sd <- ifelse(region_argo$sd > region_argo$mean, region_argo$mean, region_argo$sd)

codref <- read_excel("Data/Longhurst_Province_Summary.xls", 
                     range = "A17:B70", col_names = FALSE)
names(codref) <- c("code", "region")
region_argo <- left_join(region_argo, codref)

region_argo <- filter(region_argo, code != "ANTA")#delete this region because we only have 1 observation

g1 <- ggplot(region_argo)+
  geom_col(aes(x = reorder(code, mean), y = mean, fill = code))+
  geom_errorbar(aes(x = code, ymin = mean - sd, ymax = mean + sd))+
  xlab("code of longhurst oceanic bioregion")+
  geom_errorbar(aes(code, ymax = 2, ymin = 2),
                size=0.5, linetype = "longdash", inherit.aes = F, width = 1)+
  geom_errorbar(aes(code, ymax = 1, ymin = 1),
                size=0.5, linetype = "longdash", inherit.aes = F, width = 1)+
  theme(axis.text.x = element_text(angle = 45))+
  geom_text(aes(x = code, y = mean + sd + 0.5, label = nbr_of_float))+
  guides(fill = FALSE)+
  scale_fill_brewer(palette = "Set1")
g1  


#AFC####
afc_table <- na.omit(select(argo, pigments, code, micro, nano, pico,  ratio, lon.x, lat.x))
afc_table <- filter(afc_table, code != "ANTA")

afc_argo <- cca(na.omit(select(afc_table, pigments)))
test <- envfit(afc_argo, select(afc_table, micro, nano, pico, ratio))

env_arrow <- as.data.frame(test$vectors$arrows) #On récupère les donnes du modele sur les variables environnementales

argo_score <- as.data.frame(scores(afc_argo, choices = c(1,2,3,4,5), display = c("site"))) #this is a dataframe with the score on the 5 axes of each sample

afc_table <- bind_cols(afc_table, argo_score)

pig_score <- as.data.frame(scores(afc_argo, choices = c(1,2,3,4,5), display = c("species")))


g2 <- ggplot()+
  geom_point(aes(x = CA1, y = CA2, colour = code), size = 2, data = afc_table)+
  geom_segment(aes(x = 0, xend = CA1 *1.5, y = 0, yend = CA2*1.5), data = pig_score)+
  geom_text_repel(aes(x = CA1*1.5, y = CA2*1.5, label = rownames(pig_score)), data = pig_score)+
  geom_segment(aes(x = 0, y = 0, xend = CA1*1.7, yend = CA2*1.7), data = env_arrow, colour = "#33a02c")+
  geom_text(aes(x = CA1*1.7, y = CA2*1.7, label=rownames(env_arrow), fontface = 2), data = env_arrow)+
  scale_color_brewer(palette = "Set1") + coord_equal() +
  guides(colour = FALSE)

g2
grid.arrange(g1,g2, ncol = 2)

g3 <- ggplot(afc_table)+
  geom_point(aes(x = lon.x, y = lat.x, colour = code), size = 2)+
  geom_polygon(aes(x = long, y = lat, group = group), data = map_vec)+
  xlab("lon")+ylab("lat")+
  coord_quickmap()+
  scale_color_brewer(palette = "Set1")+
  theme_bw()

g3



grid.newpage()
pushViewport(viewport(layout = grid.layout(2,2)))
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
}
print(g1, vp = define_region(1,1))
print(g2, vp = define_region(1:2, 2))
print(g3, vp = define_region(2,1))


