library(tidyverse)
library(janitor)
library(readxl)
library(patchwork)
library(vegan)
library(ggrepel)
library(treemap)

lov <- read_excel("Dataset_LOV.xls", na = "NA") %>% clean_names()



#lov_nested <- lov %>% nest(pigments = c(chla : tot_car), aph = c(x400 : x700)) %>% 
  # mutate(plot_aph = map(aph, function(.x){pivot_longer(.x, 1:151, names_to = "wavelength", values_to = "abs")})) %>% 
  # mutate(plot_aph = map(plot_aph, ~.x %>% mutate(wavelength = as.numeric(substr(wavelength, 2, 4))))) %>% 
  # mutate(plot_aph = map2(campagne, plot_aph, function(.x,.y){
  #   ggplot(data = .y, aes(x = wavelength, y = abs))+
  #            geom_path()+
  #            theme_bw()+
  #            ggtitle(label = .x)
  # }))


#lov_nested$plot_aph[[25]] + lov_nested$plot_aph[[26]]

spectre <- read_excel("Biosope/Data/Spectres_annick.xlsx")
spectre <- clean_names(spectre)
spectre440<- filter(spectre, lambda == 440)
spectre470 <- filter(spectre, lambda == 470)

lov_pig <- select(lov, chla:tot_car)

lov_pig$zea <- as.numeric(lov_pig$zea)
lov_pig$lut <- as.numeric(lov_pig$lut)

lov_pig <- select(lov_pig, chlc12 = chlc1c2, peri, x19bf, fuco, x19hf, diad = diadino, allox = allo, zea, dv_chlb, chl_b = chlb, dv_chla, chl_a = chla, ss_car = b_car, a_car)
lov_pig[is.na(lov_pig)] <- 0
lov_pur <- select(lov_pig, - zea, -diad, - ss_car, -a_car)

lov_mat <- as.matrix(lov_pig)
lov_mat <- t(lov_mat)

lov_pur_mat <- as.matrix(lov_pur)
lov_pur_mat <- t(lov_pur_mat)

spectre[is.na(spectre)] <- 0
spectre_pur <- select(spectre, -zea, -diad, -ss_car, -a_car)

spectre_mat <- as.matrix(select(spectre, - lambda))
spectre_pur_mat <- as.matrix(select(spectre_pur, -lambda))

result <- spectre_mat %*%  lov_mat
result <- t(result)

result_pur <- spectre_pur_mat %*% lov_pur_mat
result_pur <- t(result_pur)

result_df <- as.data.frame(result)
pur_df <- as.data.frame(result_pur)

colnames(result_df) = paste('a', spectre$lambda, sep = '')
colnames(pur_df) = paste('pur', spectre$lambda, sep = '')


lov_tot <- bind_cols(lov, result_df, pur_df)
lov_long <- lov_tot %>% mutate(round_depth = round(profondeur), ratio = x440/x470) %>%
  pivot_longer(c(x400:pur700), names_to = "wavelength", values_to = "abs")

lov_long$wavelength <- gsub('^([a-z]+)([0-9]{3})$', '\\1_\\2', lov_long$wavelength)
lov_long <- lov_long %>% separate(wavelength, into = c('type', 'lambda'), '_')

lov_long <- lov_long %>% pivot_wider(names_from = type, values_from = abs)

lov_long <- lov_long %>% mutate(factor = pur/a,
                                real = x * factor)

real_df <- select(lov_long, campagne, profondeur, lambda, real)
pur_df <- select(lov_long, campagne, profondeur, lambda, pur)
a_df <- select(lov_long, campagne, profondeur, lambda, a)
x_df <- select(lov_long, campagne, profondeur, lambda, x)

transform_wide <- function(data, value, name){
  df <- data %>% group_by(lambda) %>% 
    mutate(row = row_number()) %>% 
    pivot_wider(names_from = lambda, values_from = {{value}}) %>% 
    select(-row, - campagne, - profondeur)
  colnames(df) <- gsub(' ', '',paste(name, colnames(df), sep = ''))
  return(df)
}

real_df <- transform_wide(real_df, real, 'real')
pur_df <- transform_wide(pur_df, pur, 'pur')
a_df <- transform_wide(a_df, a, 'a')
x_df <- transform_wide(x_df, x, 'x')

lov_tot <- bind_cols(lov, real_df, a_df, pur_df)

lov$campagne <- sub("[1-9](.*)", "", lov$campagne)
lov_campagne <- lov %>%
  mutate(round_depth = round(profondeur), ratio = x440/x470) %>%
  group_by(campagne, round_depth) %>%
  summarise_at(vars(x400:x700, ratio, z_zeu), c(mean, sd)) %>%
  pivot_longer(c(x400_fn1:x700_fn1, x400_fn2:x700_fn2), names_to = "wavelength", values_to = "Abs")
lov_campagne$wavelength <- substr(lov_campagne$wavelength, 2,8)
lov_campagne <- separate(lov_campagne, wavelength, into = c('wavelength', 'operation'), '_')
lov_campagne$wavelength <- as.numeric(lov_campagne$wavelength)

lov_campagne_plot <- lov_campagne %>% pivot_wider(names_from = operation, values_from = Abs) %>%
  arrange(campagne, round_depth, wavelength)
names(lov_campagne_plot) <- c('campagne', 'depth', 'ratio_mean', 'z_zeu_mean', 'ratio_sd', 'z_zeu_sd', 'wavelength', 'mean', 'sd')


ggplot(filter(lov_campagne_plot, campagne != 'BENCAL'))+
  geom_path(aes(x = wavelength, y = mean, colour = - depth, group = depth), size = 1)+
  theme_dark()+
  scale_color_distiller(palette = 'YlGnBu', direction = 1, name = 'profondeur')+
  ylab('mean aph')+
  facet_wrap(.~campagne, scales = 'free_y')

ggplot(filter(lov_campagne_plot, campagne != 'BENCAL' & z_zeu_mean < 4))+
  geom_path(aes(x = wavelength, y = mean, colour = ratio_mean, group = depth), size = 1)+
  theme_dark()+
  scale_color_distiller(palette = 'RdYlBu', name = 'ratio a40/a470')+
  ylab('mean aph')+
  facet_wrap(.~campagne, scales = 'free_y')

ggplot(filter(lov_campagne_plot, campagne != 'BENCAL' & z_zeu_mean < 4))+
  geom_path(aes(x = wavelength, y = mean, colour = z_zeu_mean, group = depth), size = 1)+
  theme_dark()+
  scale_color_binned(name = 'Z/Zeu')+
  ylab('mean aph')+
  facet_wrap(.~campagne, scales = 'free_y')


ggplot(lov_campagne_plot)+
  geom_boxplot(aes(y = ratio_mean, x = campagne))

lov_ratio <- lov_campagne_plot %>% mutate(ratio_type = round(ratio_mean, 1)) %>% 
  arrange(wavelength, depth, campagne)

ggplot(filter(lov_ratio, campagne != 'BENCAL' & z_zeu_mean < 4))+
  geom_path(aes(x = wavelength, y = mean, colour = depth, group = depth), size = 1)+
  theme_dark()+
  scale_color_distiller(palette = 'RdYlBu', name = 'profondeur')+
  ylab('mean aph')+
  xlim(400,550)+
  facet_wrap(.~ratio_type, scales = 'free_y')

lov_afc <- lov_tot %>%
  mutate(ratio_440_470 = x440/x470, ratio_440_530 = x440/x530, real_440_470 = real440/real470) %>% 
  select(campagne, lat, lon, profondeur, z_zeu, p_pico, p_nano, p_micro, fuco, peri, x19hf, x19bf, allo, t_chlb, t_chla, zea, ratio_440_470, ratio_440_530, real_440_470, x400:x600, real400:real600, a400:a600) %>% 
  mutate(rowsum = rowSums(select(., x400:x550))) %>% 
  filter(rowsum > 0 & ratio_440_530 >= 0 & ratio_440_530 < 10 & ratio_440_470 <= 1.4 & real_440_470 <= 4)


AFC <- cca(select(lov_afc, x430:x480), scale = TRUE)

scores <- data.frame(scores(AFC, choices = c(1,2,3), display = "site"))
lov_afc <- bind_cols(lov_afc, scores)

fitscore <- envfit(AFC, select(lov_afc, lat:real_440_470), na.rm = TRUE)  
fitarrow <- as.data.frame(fitscore$vectors$arrows)

ggplot(lov_afc)+
  geom_point(aes(x = CA1, y = CA2, colour = ratio_440_470))+
  geom_segment(aes(x = 0, y = 0, xend = CA1*1.7, yend = CA2*1.7), data = fitarrow)+
  geom_text_repel(aes(x = CA1*1.7, y = CA2*1.7, label=rownames(fitarrow), fontface = 2), data = fitarrow)+
  scale_color_viridis_c()+
  xlim(-3,3)+
  ylim(-10,10)

distlov <- dist(select(lov_afc, CA1, CA2, CA3))
lov_afc$group <- as.factor(cutree(hclust(distlov, method = "ward.D"), k = 3))

ggplot(lov_afc)+
  geom_point(aes(x = CA1, y = CA2, colour = group))+
  geom_segment(aes(x = 0, y = 0, xend = CA1*1.7, yend = CA2*1.7), data = fitarrow, colour = "#33a02c")+
  geom_text_repel(aes(x = CA1*1.7, y = CA2*1.7, label=rownames(fitarrow), fontface = 2), data = fitarrow)+
  scale_color_viridis_d()+
  xlim(-3,3)+
  ylim(-10,10)

lov_clust <- lov_afc %>%
  group_by(group) %>% 
  summarise_at(vars(c(x400:x600, real400:real600, a400:a600, ratio_440_470)), c(mean, sd)) %>% 
  pivot_longer(c(x400_fn1:x600_fn1, x400_fn2:x600_fn2, a400_fn1:a600_fn1, a400_fn2:a600_fn2, real400_fn1:real600_fn1, real400_fn2:real600_fn2), names_to = 'wavelength', values_to = 'abs') %>% 
  separate(wavelength, into = c('wavelength', 'operation'), '_')

lov_clust$lambda <- as.numeric(str_sub(lov_clust$wavelength, -3,-1))
lov_clust$operation <- gsub('fn1', 'mean', lov_clust$operation)
lov_clust$operation <- gsub('fn2', 'sd', lov_clust$operation)
lov_clust$type <- str_sub(lov_clust$wavelength, end = -4)


lov_clust <- lov_clust %>% 
  pivot_wider(names_from = operation, values_from = abs)

fit <- aov(ratio_440_470~group , data = lov_afc)
hsd <- TukeyHSD(fit)

summary <- lov_afc %>% 
  group_by(group) %>% 
  summarise_at(vars(c(ratio_440_470, real_440_470)), c(mean, sd))

ggplot(filter(lov_clust, type == 'x'))+
  geom_path(aes(x = lambda, y = mean, colour = 'observed'))+
  geom_path(aes(x = lambda, y = mean, colour = 'real'), data = filter(lov_clust, type == 'real'))+
  geom_line(aes(x = lambda, y = mean + sd, colour = 'observed'), linetype = 'dotted')+
  geom_line(aes(x = lambda, y = mean - sd, colour = 'observed'), linetype = 'dotted')+
  geom_line(aes(x = lambda, y = mean + sd, colour = 'real'), linetype = 'dotted', data = filter(lov_clust, type == 'real'))+
  geom_line(aes(x = lambda, y = mean - sd, colour = 'real'), linetype = 'dotted', data = filter(lov_clust, type == 'real'))+
  geom_vline(xintercept = 440, colour = 'blue')+
  geom_vline(xintercept = 470, colour = 'green')+
  facet_wrap(.~ group, scales = 'free_y')

ggplot(lov_afc)+
  geom_boxplot(aes(x = group, y = ratio_440_470))

tplot <- lov_afc %>% 
  group_by(group) %>% 
  summarise_at(vars(c(p_pico, p_nano, p_micro, t_chlb, fuco, zea, peri, allo, x19hf, x19bf)), mean, na.rm = TRUE) %>% 
  ungroup() %>% 
  pivot_longer(t_chlb:x19bf, names_to = 'pigment', values_to = 'concentration') %>% 
  mutate(size = ifelse(pigment %in% c('zea', 't_chlb'), 'pico', ifelse(pigment %in% c('allo', 'x19hf', 'x19bf'), 'nano', ifelse(pigment %in% c('fuco', 'peri'), 'micro', 'error'))))

tplot1 <- filter(tplot, group == '1')
tplot2 <- filter(tplot, group == '2')
tplot3 <- filter(tplot, group == '3')


treemap(tplot1, index = c('size', 'pigment'), vSize = 'concentration', type = 'index', palette = 'Set1')
treemap(tplot2, index = c('size', 'pigment'), vSize = 'concentration', type = 'index', palette = 'Set1')
treemap(tplot3, index = c('size', 'pigment'), vSize = 'concentration', type = 'index', palette = 'Set1')

ggplot(tplot1, aes(area = concentration, fill = size, subgroup = pigment))+
  geom_treemap(position = )+
  geom_treemap_subgroup_text(size = 12)

par(mfrow = c(1,3))

ggplot(lov_afc)+
  geom_boxplot(aes(x = group, y = real_440_470))

fit <- aov(real_440_470~group , data = lov_afc)
hsd <- TukeyHSD(fit)


