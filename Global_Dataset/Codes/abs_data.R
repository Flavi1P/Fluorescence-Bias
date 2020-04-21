library(tidyverse)
library(janitor)
library(readxl)
library(patchwork)
library(vegan)
library(ggrepel)
library(treemap)
library(treemapify)

map_vec <- read_csv('Data/map_vec')

#lov <- read_excel("Dataset_LOV.xls", na = "NA") %>% clean_names()
types <- c('c', 'c', rep('n', 186))
lov <- read_csv('Data/Absorbtion/lov_soclim_peacetime.csv', col_types = as.list(types))


# lov_nested <- lov %>% nest(pigments = c(chla : tot_car), aph = c(x400 : x700)) %>% 
#   mutate(plot_aph = map(aph, function(.x){pivot_longer(.x, 1:151, names_to = "wavelength", values_to = "abs")})) %>% 
#   mutate(plot_aph = map(plot_aph, ~.x %>% mutate(wavelength = as.numeric(substr(wavelength, 2, 4))))) %>% 
#   mutate(plot_aph = map2(campagne, plot_aph, function(.x,.y){
#   ggplot(data = .y, aes(x = wavelength, y = abs))+
#               geom_path()+
#               theme_bw()+
#               ggtitle(label = .x)
#    }))


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
lov_long <- lov_tot %>% mutate(round_depth = round(depth), ratio = x440/x470) %>%
  pivot_longer(c(x400:pur700), names_to = "wavelength", values_to = "abs")

lov_long$wavelength <- gsub('^([a-z]+)([0-9]{3})$', '\\1_\\2', lov_long$wavelength)
lov_long <- lov_long %>% separate(wavelength, into = c('type', 'lambda'), '_')

lov_long <- lov_long %>% pivot_wider(names_from = type, values_from = abs)

lov_long <- lov_long %>% mutate(factor = pur/a,
                                real = x * factor)

real_df <- select(lov_long, campagne, depth, lambda, real)
pur_df <- select(lov_long, campagne, depth, lambda, pur)
a_df <- select(lov_long, campagne, depth, lambda, a)
x_df <- select(lov_long, campagne, depth, lambda, x)

transform_wide <- function(data, value, name){
  df <- data %>% group_by(lambda) %>% 
    mutate(row = row_number()) %>% 
    pivot_wider(names_from = lambda, values_from = {{value}}) %>% 
    select(-row, - campagne, - depth)
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
  mutate(round_depth = round(depth), ratio = x440/x470) %>%
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


ggplot(lov_campagne_plot)+
  geom_boxplot(aes(y = ratio_mean, x = campagne))

lov_afc <- lov_tot %>%
  mutate(ratio_440_470 = x440/x470, ratio_440_530 = x440/x530, real_440_470 = real440/real470) %>% 
  select(campagne, station, lat, lon, depth, z_zeu, p_pico, p_nano, p_micro, fuco, peri, x19hf, x19bf, allo, t_chlb, t_chla, zea, ratio_440_470, ratio_440_530, real_440_470, x400:x600, real400:real600, a400:a600) %>% 
  mutate(rowsum = rowSums(select(., x400:x550)),
         rowsum2 = rowSums(select(., real400:real600))) %>% 
  filter(rowsum > 0 & rowsum2 > 0 & ratio_440_530 >= 0 & ratio_440_530 < 10 & ratio_440_470 <= 1.4 & real_440_470 <= 4)


AFC <- cca(select(lov_afc, real400:real600), scale = TRUE)

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

distlov <- dist(select(lov_afc, CA1, CA2))
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

g1 <- ggplot(filter(lov_clust, type == 'x'))+
  geom_path(aes(x = lambda, y = mean, colour = 'observed'))+
  geom_path(aes(x = lambda, y = mean, colour = 'real'), data = filter(lov_clust, type == 'real'))+
  geom_line(aes(x = lambda, y = mean + sd, colour = 'observed'), linetype = 'dotted')+
  geom_line(aes(x = lambda, y = mean - sd, colour = 'observed'), linetype = 'dotted')+
  geom_line(aes(x = lambda, y = mean + sd, colour = 'real'), linetype = 'dotted', data = filter(lov_clust, type == 'real'))+
  geom_line(aes(x = lambda, y = mean - sd, colour = 'real'), linetype = 'dotted', data = filter(lov_clust, type == 'real'))+
  geom_vline(xintercept = 440, colour = 'blue')+
  geom_vline(xintercept = 470, colour = 'green')+
  facet_wrap(.~ group,scales = 'free_y')+
  theme_bw()


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





g2 <- ggplot(tplot1, aes(area = concentration, fill = size, subgroup = size, label = pigment))+
  geom_treemap(layout = 'fixed')+
  geom_treemap_subgroup_text(layout = 'fixed', place = 'middle', fontface = 'bold', size = 14)+
  geom_treemap_text(layout = 'fixed', place = 'bottomright', 'size' = 11, colour = 'white', fontface = 'italic')+
  guides(fill = FALSE)+
  scale_fill_brewer(palette = 'Dark2')

g3 <- ggplot(tplot2, aes(area = concentration, fill = size, subgroup = size, label = pigment))+
  geom_treemap(layout = 'fixed')+
  geom_treemap_subgroup_text(layout = 'fixed', place = 'middle', fontface = 'bold', size = 14)+
  geom_treemap_text(layout = 'fixed', place = 'bottomright', 'size' = 11, colour = 'white', fontface = 'italic')+
  guides(fill = FALSE)+
  scale_fill_brewer(palette = 'Dark2')

g4 <- ggplot(tplot3, aes(area = concentration, fill = size, subgroup = size, label = pigment))+
  geom_treemap(layout = 'fixed')+
  geom_treemap_subgroup_text(layout = 'fixed', place = 'middle', fontface = 'bold', size = 14)+
  geom_treemap_text(layout = 'fixed', place = 'bottomright', 'size' = 11, colour = 'white', fontface = 'italic')+
  guides(fill = FALSE)+
  scale_fill_brewer(palette = 'Dark2')

g1 /(g2 | g3 | g4)


lov_afc %>% select(ratio_440_470, real_440_470) %>% 
  pivot_longer(c(1,2), names_to = 'ratio', values_to = 'value') %>% 
  ggplot()+
  geom_violin(aes(y = value, x = ratio))

comp_raw <- lov_afc %>% select(real440, x440, real470, x470) %>% 
  pivot_longer(c(1:4), names_to = 'wl', values_to = 'abs')

model <- aov(abs~wl, data = comp_raw)
summary(model)

ggplot(lov_afc)+
  geom_bar(aes(x = campagne, fill = group), position = 'fill')

biosope <- filter(lov_afc, campagne == 'Biosope')

ggplot(biosope)+
  geom_point(aes(x = lon, y = - depth, colour = group), size = 2)

biosope %>% filter(station == "UPW1") %>% 
  ggplot()+
  geom_path(aes(x = x440, y = - depth, colour = 'x440'))+
  geom_path(aes(x = x470, y = - depth, colour = 'x470'))+
  geom_path(aes(x = real470, y = - depth, colour = 'real470'))+
  geom_path(aes(x = real440, y = - depth, colour = 'real440'))

biosope %>% filter(station == "STB14") %>% 
  ggplot()+
  geom_path(aes(x = ratio_440_470, y = - depth, colour = 'ratio'))+
  geom_path(aes(x = real_440_470, y = - depth, colour = 'real'))

biosope %>% filter(station == "EGY3") %>% 
  ggplot()+
  geom_point(aes(x = real_440_470, y = - depth, colour = group), size = 2)

bio_nested <- biosope %>% nest(aph = c(x400 : x600), real = c(real400:real600)) %>% 
mutate(plot_aph = map(aph, function(.x){pivot_longer(.x, 1:101, names_to = "wavelength", values_to = "abs")})) %>% 
mutate(plot_aph = map(plot_aph, ~.x %>% mutate(wavelength = as.numeric(substr(wavelength, 2, 4))))) %>% 
mutate(plot_aph = map2(depth, plot_aph, function(.x,.y){
ggplot(data = .y, aes(x = wavelength, y = abs))+
            geom_path()+
    geom_vline(aes(xintercept = 440), colour = 'Blue')+
    geom_vline(aes(xintercept = 470), colour = 'Green')+
          theme_bw()+
           ggtitle(label = .x)
}))


bio_nested %>% filter(group == 2 & lon < -130) %>% 
  .[['plot_aph']] %>% reduce(`%+%`)

bio_nested %>% filter(station == 'EGY3') %>% 
  .[['plot_aph']] %>% reduce(`%+%`)

soclim <- filter(lov_afc, campagne == 'soclim')

ggplot(soclim)+
  geom_point(aes(x = lat, y = - depth, colour = group), size = 2)+
  ylim(-210,0)

occur <- soclim %>% group_by(station) %>% count(group) %>% 
  mutate(max = max(n)) %>% filter(n == max) %>% 
  select(station, dominent = group) %>%
  left_join(soclim)

ggplot(occur)+
  geom_text(aes(x = lon, y = lat, label = station, colour = dominent))+
  geom_polygon(aes(x= long, y = lat, group = group), data = map_vec)+
  coord_quickmap(xlim = c((min(soclim$lon) - 5), max(soclim$lon) + 5), ylim = c(min(soclim$lat) - 5, max(soclim$lat) + 5))

ggplot(occur)+
  geom_path(aes(x = real_440_470, y = - depth, group = station))+
  facet_wrap(.~dominent)

treeso <- soclim %>% group_by(station) %>% 
  summarise_at(vars(c(p_pico, p_nano, p_micro, t_chlb, fuco, zea, peri, allo, x19hf, x19bf)), mean, na.rm = TRUE) %>% 
  ungroup() %>% 
  pivot_longer(t_chlb:x19bf, names_to = 'pigment', values_to = 'concentration') %>% 
  mutate(size = ifelse(pigment %in% c('zea', 't_chlb'), 'pico', ifelse(pigment %in% c('allo', 'x19hf', 'x19bf'), 'nano', ifelse(pigment %in% c('fuco', 'peri'), 'micro', 'error'))))


ggplot(treeso, aes(area = concentration, fill = size, subgroup = size, label = pigment))+
  geom_treemap(layout = 'fixed')+
  geom_treemap_subgroup_text(layout = 'fixed', place = 'middle', fontface = 'bold', size = 14)+
  geom_treemap_text(layout = 'fixed', place = 'bottomright', 'size' = 11, colour = 'white', fontface = 'italic')+
  guides(fill = FALSE)+
  scale_fill_brewer(palette = 'Dark2')+
  facet_wrap(.~station)

soclim_long <- lov_long %>% filter(campagne == 'soclim') %>% 
  arrange(lambda) %>% 
  arrange(station) %>% 
  mutate(id = paste(station, depth, sep = '_'))

ggplot(soclim_long)+
  geom_path(aes(x = lambda, y = factor, colour = station, group = id))

biosope_long <- lov_long %>% filter(campagne == 'Biosope') %>% 
  arrange(lambda) %>% 
  arrange(station) %>% 
  mutate(id = paste(station, depth, sep = '_'))

biosope_long$lambda <- as.numeric(biosope_long$lambda)
ggplot(biosope_long)+
  geom_path(aes(x = lambda, y = factor, colour = station, group = id))+
  scale_x_continuous(breaks = c(410, 440, 470, 550))

ggplot(biosope_long)+
  geom_path(aes(x = lambda, y = factor, colour = station, group = id))+
  scale_x_continuous(breaks = c(4300, 440, 470, 500, 550), limits = c(420, 600))

model_440 <- lm(log(t_chla)~log(a440), data = filter(lov_afc, group == '1'))
summary(model_440)                     
model_470 <- lm(t_chla~a470, data = filter(lov_afc, group == '1'))
summary(model_470)                     

ggplot(lov_afc)+
  geom_point(aes(x = t_chla, y = real440, colour = '440'))+
  geom_point(aes(x = t_chla, y = real470, colour = '470'))+
  coord_trans(x = 'log', y = 'log')+
  facet_wrap(.~group)

lov_afc <- lov_afc %>% mutate(chla_470 = real470/t_chla,
                              chla_440 = real440/t_chla)

rsq <- data.frame('chl' = seq(0.005, 5, 0.005),
                  'real_440' = NA,
                  'real_470' = NA,
                  'observed_440' = NA,
                  'observed_470' = NA,
                  'vitro_440' = NA,
                  'vitro_470' = NA)

for(i in c(1:length(rsq$chl))){
  threshold <- rsq$chl[i]
  df <- filter(lov_afc, t_chla < threshold)
  real_470 <- lm(t_chla~real470, data = df)
  real_440 <- lm(t_chla~real440, data = df)
  real_440 <- summary(real_440)$adj.r.squared
  real_470 <- summary(real_470)$adj.r.squared
  observed_470 <- summary(lm(t_chla~x470, data = df))$adj.r.squared
  observed_440 <- summary(lm(t_chla~x440, data = df))$adj.r.squared
  vitro_470 <- summary(lm(t_chla~a470, data = df))$adj.r.squared
  vitro_440 <- summary(lm(t_chla~a440, data = df))$adj.r.squared
  rsq$real_440[i] <- real_440
  rsq$real_470[i] <- real_470
  rsq$observed_440[i] <- observed_440
  rsq$observed_470[i] <- observed_470
  rsq$vitro_440[i] <- vitro_440
  rsq$vitro_470[i] <- vitro_470
}

ggplot(rsq)+
  geom_path(aes(x = chl, y = real_440, colour = 'Real 440'))+
  geom_path(aes(x = chl, y = real_470, colour = 'Real 470'))+
  geom_path(aes(x = chl, y = observed_440, colour = 'observed 440'))+
  geom_path(aes(x = chl, y = observed_470, colour = 'observed 470'))+
  geom_path(aes(x = chl, y = vitro_440, colour = 'vitro 440'))+
  geom_path(aes(x = chl, y = vitro_470, colour = 'vitro 470'))
