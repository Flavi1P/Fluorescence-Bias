library(tidyverse)
library(janitor)
library(readxl)
library(patchwork)
library(vegan)
library(ggrepel)
library(treemap)
library(treemapify)

types <- c('c', 'c', rep('n', 327))
lov_afc <- read_csv('Data/lov_afc.csv', col_types = as.list(types))

lov_afc <- lov_afc %>% mutate(chla_470 = real470/t_chla,
                              chla_440 = real440/t_chla)

rsq <- data.frame('chl' = seq(0.02, 1, 0.05),
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
  geom_path(aes(x = chl, y = vitro_470, colour = 'vitro 470'))+
  scale_color_brewer(palette = 'Paired')+
  ylab('R² value')+
  xlab('Chla concentration')+
  ggtitle('R² of the relation aphy~chla')+
  xlim(0,1)


model <- lm(real440~t_chla, data = lov_afc)
test <- lov_afc %>% arrange(t_chla)

test <- lov_afc %>%
  filter(t_chla < 1) %>% 
  mutate(t_chla = plyr::round_any(t_chla, 0.1))
table(test$t_chla)

test <- test %>% group_by(t_chla) %>% 
  sample_n(size = 25)

test$cv440 <- runsd(test$a440, 40)/runmean(test$a440, 40)
test$cv470 <- runsd(test$a470, 40)/runmean(test$a470, 40)

test_summarised <- test %>% group_by(t_chla) %>% summarise_at(c('cv440', 'cv470'), mean)
ggplot(test_summarised)+
  geom_point(aes(x = t_chla ,y = cv440, colour = '440'))+
  geom_point(aes(x = t_chla, y = cv470, colour = '470'))

lov_model <- lov_afc %>% filter(t_chla < 3 & z_zeu <= 3)

model_440 <- lm(t_chla~real440, data = lov_model)
model_470 <- lm(t_chla~real470, data = lov_model)

AIC(model_440)
AIC(model_470)
exp((AIC(model_440)-AIC(model_470))/2)

df_model <- data.frame('chla' = lov_model$t_chla , 'resid440' = model_440$residuals, 'resid470' = model_470$residuals)

ggplot(df_model)+
  geom_smooth(aes(x = chla, y = abs(resid440), colour = '440'), method = 'glm')+
  geom_smooth(aes(x = chla, y = abs(resid470), colour = '470'), method = 'glm')

ggplot(lov_model)+
  geom_point(aes(x = t_chla, y = real470/t_chla))+
  xlim(0,0.1)
