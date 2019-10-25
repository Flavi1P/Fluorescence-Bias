library(tidyverse)
library(sf)
library(readxl)
library(vegan)
library(ggrepel)
library(gridExtra)
library(grid)
library(nnls)
library(FactoMineR)
library(janitor)
path = "Data/Longhurst"
argo <- read_csv("Data/merged_argo")
source("functions/outliers.R")
source("functions/phi_simple.R")
map_vec <- read_csv("Data/map_vec")
pigments <- c("fuco", "peri", "hex", "but", "allo", "tchlb", "zea")
NAT_IRS_list <- c("lovbio059c", "lovbio045b", "lovbio024c", "lovbio044b", "lovbio031c", "lovbio027b", "lovbio040b", "lovbio026c")


argo <- filter(argo, !(lovbio == "takapm005b" & depth == 20))#remove an hplc match in a spike
argo <- filter(argo, lovbio != "lovbio067c" & lovbio != "lovbio083d" & lovbio != "lovbio085d" & lovbio != "lovbio090d") #filter profile with only 1 match
argo <- filter(argo, !(lovbio %in% NAT_IRS_list)) #filter dubious hplc in arctique


longhurst_sf <- read_sf(dsn = path.expand(path), quiet = TRUE)

names(longhurst_sf) <- c("code", "region", "geometry")

#longhurst_sf %>% ggplot() + geom_sf(aes(fill = code))

pnts_sf <- do.call("st_sfc",c(lapply(1:nrow(argo),
                                     function(i) {st_point(as.numeric(argo[i,c("lon.y", "lat.y") ]))}), list("crs" = 4326))) 
pnts_trans <- st_transform(pnts_sf, 4326)
longhurst_trans <- st_transform(longhurst_sf, 4326)  
argo$code <- apply(st_intersects(longhurst_trans, pnts_trans, sparse = FALSE), 2, 
                     function(col) { 
                       longhurst_trans[which(col), ]$code
                     })

ggplot(argo)+
  geom_point(aes(x = lon.x, y = lat.x), size = 4, shape = 21, fill = "red2")+
  geom_polygon(aes(x = long, y = lat, group = group), data = map_vec, fill = "gray28")+
  theme_bw(base_size = 16)+
  xlab("Longitude (°E)")+
  ylab("Latitude (°N)")+
  coord_quickmap()
#ggsave("argo/Plots/map.png")




argo <- filter(argo, optical_layer < 4)



nbr_match_region<- count(argo, "code")
nbr_float_region <- argo %>% filter(duplicated(lovbio) == FALSE) %>% count("code")
resume_region <- bind_cols(nbr_float_region, nbr_match_region) %>% dplyr::select(1,2,4)

names(resume_region) <- c("code", "nbr_of_float", "nbr_of_match")


argo <- argo %>% mutate(fluo = chla_adjusted * 2,
                        ratio = fluo/tchla)

influential <- outliers(select(argo, fluo, tchla, micro, nano, pico))

#argo <- argo[-influential,]

#lets mean the 7 profile of sarc because they all correspond to the same launch

sarc <- argo %>% filter(code == "SARC") %>% group_by(depth) %>% summarise_all(mean) %>% ungroup()
sarc$code <- "SARC"
sarc$lovbio <- "lovbio_sarc"

argo <- argo %>% filter(code != "SARC")
argo <- bind_rows(argo, sarc)


phi_multiple <- argo %>% group_by(code) %>% phi_simple("fluo")

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
  xlab("Province océanique")+
  ylab("Rapport Fluo/Chla")+
  geom_errorbar(aes(code, ymax = 1, ymin = 1),
                size=1, linetype = "longdash", inherit.aes = F, width = 1)+
  guides(fill = FALSE)+
  scale_fill_brewer(palette = "Dark2")+
  theme_bw(base_size = 20)+
  ylim(0,9)
g1  

#compute the absorbance ratio

#open pigments absorbtion
spectre <- read_excel("Biosope/Data/Spectres_annick.xlsx")
spectre <- clean_names(spectre)

#select the absorbtion at 440 and 470nm
spectre440<- filter(spectre, lambda == 440)
spectre470 <- filter(spectre, lambda == 470)

#create columns that correspond to the total photosynthetic absorbance and non photosynthetic absorbance at 440 and 470. Create also a ratio between the two photosynthetic absorbtion
argo<- argo %>% mutate(photo_440 = peri * spectre440$peri + but * spectre440$x19_bf + hex * spectre440$x19_hf + fuco * spectre440$fuco + allo * spectre440$allox + tchla * spectre440$chl_a,
                              protect_440 = zea * spectre440$zea,
                              photo_470 = peri * spectre470$peri + but * spectre470$x19_bf + hex * spectre470$x19_hf + fuco * spectre470$fuco + allo * spectre470$allox + tchla * spectre470$chl_a,
                              protect_470 = zea * spectre470$zea, 
                              ratio_abs = photo_440/photo_470)

#resume it by region
region_argo_absorbance <- argo  %>% group_by(code) %>% summarise_at(vars(ratio_abs), c(mean, sd), na.rm = TRUE) %>% ungroup()
names(region_argo_absorbance) <- c("code", "mean", "sd")

region_argo_absorbance <- filter(region_argo_absorbance, code != "ANTA")#delete this region because we only have 1 observation

#plot

gabs <- ggplot(region_argo_absorbance)+
  geom_col(aes(x = reorder(region_argo$code, region_argo$mean), y = mean, fill = code))+
  geom_errorbar(aes(x = code, ymin = mean - sd, ymax = mean + sd))+
  xlab("Province océanique")+
  ylab("Rapport a440/a470")+
  guides(fill = FALSE)+
  scale_fill_brewer(palette = "Dark2")+
  theme_bw(base_size = 20)

grid.arrange(g1, gabs, ncol = 1)


ggplot(argo)+
  geom_point(aes(x = ratio, y = ratio_abs, colour = code))+
  xlab("rapport fluo/[tchla]")+
  ylab("rapport a440/a470")+
  ylim(1,5)+
  xlim(0,10)+
  theme_minimal()

#variance des deux absorbtions photosynthétiques

#create a df to be able to use the geom boxplot function
absorbtion <- data.frame("lambda" = as.factor(c(rep(440, 180), rep(470,180))), "aps" = c(argo$photo_440, argo$photo_470))

ggplot(absorbtion)+
  geom_boxplot(aes(x = lambda, y = aps))+
  ylab("aPS")+
  theme_classic()

#analyse des pigments

pigment_region <- argo %>% group_by(code) %>%
  summarise_at(vars(peri, but, hex, fuco, allo, tchla), mean) %>% ungroup() %>% 
  gather(2:7, key = "pigment", value = "concentration")

g_pig <- ggplot(filter(pigment_region, pigment != "tchla"))+
  geom_col(aes(x = code, y = concentration, fill = pigment), position = "Fill")+
  scale_fill_brewer(palette = "Dark2")+
  theme_bw(base_size = 20)

grid.arrange(g1, g_pig, ncol = 1)


#compute the size index 

argo <- argo %>% mutate(size_index = 1 * micro + 5 * nano + 50 * micro)

ggplot(argo)+
  geom_point(aes(x = size_index, y = photo_440, colour = code))+
  xlab("size index")+
  ylab("rapport a440/a470")+
  theme_minimal()

#resume it by region
region_argo_size_index <- argo  %>% group_by(code) %>% summarise_at(vars(size_index), c(mean, sd), na.rm = TRUE) %>% ungroup()
names(region_argo_size_index) <- c("code", "mean", "sd")

region_argo_size_index <- filter(region_argo_size_index, code != "ANTA")#delete this region because we only have 1 observation

#plot

g_size_index <- ggplot(region_argo_size_index)+
  geom_col(aes(x = reorder(region_argo$code, region_argo$mean), y = mean, fill = code))+
  geom_errorbar(aes(x = code, ymin = mean - sd, ymax = mean + sd))+
  xlab("Province océanique")+
  ylab("Size index")+
  guides(fill = FALSE)+
  scale_fill_brewer(palette = "Dark2")+
  theme_bw(base_size = 20)
g_size_index

grid.arrange(g1, gabs, g_size_index, ncol = 1 )

ggplot(argo)+
  geom_point(aes(x = size_index, y = ratio, colour = code))+
  scale_color_brewer(palette = "Dark2")+
  theme_classic()
# #AFC####
# afc_table <- na.omit(select(argo, pigments, code, micro, nano, pico,  ratio, lon.y, lat.y))
# afc_table <- filter(afc_table, code != "ANTA")
# 
# afc_argo <- cca(na.omit(select(afc_table, pigments)))
# test <- envfit(afc_argo, select(afc_table, micro, nano, pico, ratio))
# 
# env_arrow <- as.data.frame(test$vectors$arrows) #On récupère les donnes du modele sur les variables environnementales
# 
# argo_score <- as.data.frame(scores(afc_argo, choices = c(1,2,3,4,5), display = c("site"))) #this is a dataframe with the score on the 5 axes of each sample
# 
# afc_table <- bind_cols(afc_table, argo_score)
# 
# pig_score <- as.data.frame(scores(afc_argo, choices = c(1,2,3,4,5), display = c("species")))
# 
# 
# ggplot()+
#   geom_point(aes(x = CA1, y = CA2, colour = code), size = 2, data = afc_table)+
#   geom_segment(aes(x = 0, xend = CA1 *1.5, y = 0, yend = CA2*1.5), data = pig_score)+
#   geom_text_repel(aes(x = CA1*1.5, y = CA2*1.5, label = rownames(pig_score)), data = pig_score)+
#   geom_segment(aes(x = 0, y = 0, xend = CA1*1.7, yend = CA2*1.7), data = env_arrow, colour = "#33a02c")+
#   geom_text(aes(x = CA1*1.7, y = CA2*1.7, label=rownames(env_arrow), fontface = 2), data = env_arrow)+
#   scale_color_brewer(palette = "Dark2) + coord_equal() +
#   guides(colour = FALSE)+
#   theme_bw()
# 
# g2
# grid.arrange(g1,g2, ncol = 2)
# 
# g3 <- ggplot(afc_table)+
#   geom_point(aes(x = lon.y, y = lat.y, fill = code), pch = 21, colour = "black", size = 4)+
#   geom_polygon(aes(x = long, y = lat, group = group), data = map_vec)+
#   xlab("Longitude (°E)")+ylab("Latitude (°N)")+
#   coord_quickmap()+
#   scale_fill_brewer(palette = "Dark2, name = "Province océanique", labels = c("Archipelagos", "Arctique Atlantique", "Boréal Polaire", "Méditerranée", "Courant Circumpolaire Antarctique", "Atlantique Subarctique", "Gyre Subtropical Pacifique Sud"))+
#   guides("fill" = FALSE)+
#   theme_bw(base_size = 18)
# 
# g3
# 
# #ggsave("argo/Plots/map_simple.png")
# 
# grid.newpage()
# pushViewport(viewport(layout = grid.layout(2,2)))
# define_region <- function(row, col){
#   viewport(layout.pos.row = row, layout.pos.col = col)
# }
# print(g1, vp = define_region(1,1))
# print(g2, vp = define_region(1:2, 2))
# print(g3, vp = define_region(2,1))
# 
# png(file = "argo/Plots/barplot_map.png", width = 8.25, height = 7.73, units = "in", res = 100)
# grid.arrange(g3,g1)
# dev.off()
# 
# #DCA####
# 
# AFC_detrend <- decorana(select(afc_table, pigments))
# AFC_detrend
# 
# scores_detrend <- data.frame(scores(AFC_detrend, choices = c(1,2,3), display = "site"))
# afc_table <- bind_cols(afc_table, scores_detrend)
# 
# pigscore_detrend <- data.frame(scores(AFC_detrend, choices = c(1,2,3), display = "species"))
# 
# fitscore_detrend <- envfit(AFC_detrend, select(afc_table, micro, nano, pico, ratio))
# fitarrow_detrend <- as.data.frame(fitscore_detrend$vectors$arrows)
# 
# 
# 
# ggplot(afc_table)+
#   geom_point(aes(x = DCA1, y = DCA2, colour = code), size = 3)+
#   geom_segment(aes(x = 0, xend = DCA1 * 3/4, y = 0, yend = DCA2 *3/4), data = pigscore_detrend)+
#   geom_text(aes(x = DCA1 * 3/4, y = DCA2 * 3/4, label = rownames(pigscore_detrend)), data = pigscore_detrend, size = 6)+
#   geom_segment(aes(x = 0, y = 0, xend = DCA1, yend = DCA2), data = fitarrow_detrend, colour = "#33a02c", size = 1)+
#   geom_text(aes(x = DCA1, y = DCA2, label=rownames(fitarrow_detrend), fontface = 2), data = fitarrow_detrend, size = 7)+
#   scale_color_brewer(name = "Code", palette = "Dark2)+
#   coord_equal()+
#   xlab("DCA1 47%")+ylab("DCA2 7%")+
# #  ggtitle("Detrend Correspondance analysis on Argo HPLC data")+
#   theme_bw(base_size = 16)
# ggsave("argo/Plots/DCA.png")
# 
# fit1 <- lm(ratio~microfluo + nanofluo + picofluo, data =argo)
# plot(fit1)
# summary(fit1)

#### phi global part (optional), phi global = fluorescent yield of the community
# 
# source("functions/phi_boot.R")
# 
# phi_argo <- phi_boot(argo, variable = "fluo")
# phi_argo$se <- ifelse(phi_argo$se > phi_argo$phi, phi_argo$phi, phi_argo$se)
# 
# ggplot(phi_argo, aes(x=size, y = phi, fill = as.factor(optical_layer))) +
#   geom_bar(position=position_dodge(), stat="identity") +
#   geom_errorbar(aes(ymin = phi-se, ymax = phi+se),position=position_dodge())+
#   scale_fill_viridis_d( name = "optical layer")+
#   ylab("Phi")+ xlab("size classe")
# 
# 
# phi_argo <- phi_argo %>% select(phi, optical_layer, size) %>% spread(key = size, value = phi)
# names(phi_argo) <- c("optical_layer", "phi_micro", "phi_nano", "phi_pico") 
# 
# 
# argo <- left_join(argo, phi_argo)
# argo <- argo %>% mutate(phi_glob = micro * phi_micro + nano * phi_nano + pico * phi_pico)
# 
# 
# phi_glob_argo <- argo  %>% group_by(code) %>% summarise_at(vars(phi_glob), c(mean, sd), na.rm = TRUE) %>% ungroup()
# names(phi_glob_argo) <- c("code", "mean", "sd")
# 
# g4 <- ggplot(phi_glob_argo)+
#   geom_col(aes(x = reorder(code, mean), y = mean, fill = code))+
#   geom_errorbar(aes(x = code, ymin = mean - sd, ymax = mean + sd))+
#   xlab("code of longhurst oceanic bioregion")+
#   geom_errorbar(aes(code, ymax = 2, ymin = 2),
#                 size=0.5, linetype = "longdash", inherit.aes = F, width = 1)+
#   geom_errorbar(aes(code, ymax = 1, ymin = 1),
#                 size=0.5, linetype = "longdash", inherit.aes = F, width = 1)+
#   theme(axis.text.x = element_text(angle = 45))+
#   guides(fill = FALSE)+
#   scale_fill_brewer(palette = "Dark2)
# 
# grid.arrange(g1, g4, ncol = 1)
# 
# afc_table <- na.omit(select(argo, pigments, code, micro, nano, pico,  ratio, lon.y, lat.y, phi_glob))
# afc_table <- filter(afc_table, code != "ANTA")
# 
# afc_argo <- cca(na.omit(select(afc_table, pigments)))
# test <- envfit(afc_argo, select(afc_table, micro, nano, pico, ratio, phi_glob))
# 
# env_arrow <- as.data.frame(test$vectors$arrows) #On récupère les donnes du modele sur les variables environnementales
# 
# argo_score <- as.data.frame(scores(afc_argo, choices = c(1,2,3,4,5), display = c("site"))) #this is a dataframe with the score on the 5 axes of each sample
# 
# afc_table <- bind_cols(afc_table, argo_score)
# 
# pig_score <- as.data.frame(scores(afc_argo, choices = c(1,2,3,4,5), display = c("species")))
# 
# 
# ggplot()+
#   geom_point(aes(x = CA1, y = CA2, colour = code), size = 2, data = afc_table)+
#   geom_segment(aes(x = 0, xend = CA1 *1.5, y = 0, yend = CA2*1.5), data = pig_score)+
#   geom_text_repel(aes(x = CA1*1.5, y = CA2*1.5, label = rownames(pig_score)), data = pig_score)+
#   geom_segment(aes(x = 0, y = 0, xend = CA1*1.7, yend = CA2*1.7), data = env_arrow, colour = "#33a02c")+
#   geom_text(aes(x = CA1*1.7, y = CA2*1.7, label=rownames(env_arrow), fontface = 2), data = env_arrow)+
#   scale_color_brewer(palette = "Dark2) + coord_equal() +
#   guides(colour = FALSE)
# 
# summary(lm(ratio~phi_glob, data = argo))
