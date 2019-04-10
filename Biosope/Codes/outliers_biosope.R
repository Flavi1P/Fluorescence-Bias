library(tidyverse)

pigments <- c("fuco", "peri", "hex", "but", "allo", "tchlb", "zea")

biosope <- read_csv("Biosope/Data/biosope")

mod <- lm(fluo_urel~tchla + micro + nano + pico, data = biosope)
cooksd <- cooks.distance(mod)

plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")
abline(h = 4*mean(cooksd, na.rm=T), col="red") 
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")

influential <- as.numeric(names(cooksd)[(cooksd > 4*mean(cooksd, na.rm=T))]) 
biosope <- biosope[-influential,]

biosope <- select(biosope, pigments, fluo_urel, tchla, micro, nano, pico, ratio)

AFC <- cca(select(biosope, pigments))

scores <- data.frame(scores(AFC, choices = c(1,2,3), display = "site"))
biosope <- bind_cols(biosope, scores)

pigscore <- data.frame(scores(AFC, choices = c(1,2,3), display = "species"))



fitscore <- envfit(AFC, select(biosope, micro, nano, pico, ratio))
fitarrow <- as.data.frame(fitscore$vectors$arrows)

ggplot(biosope)+
  geom_point(aes(x = CA1, y = CA2, colour = pico))+
  geom_segment(aes(x = 0, xend = CA1, y = 0, yend = CA2), data = pigscore)+
  geom_text(aes(x = CA1, y = CA2, label = rownames(pigscore)), data = pigscore)+
  geom_segment(aes(x = 0, y = 0, xend = CA1*1.7, yend = CA2*1.7), data = fitarrow, colour = "#33a02c")+
  geom_text(aes(x = CA1*1.7, y = CA2*1.7, label=rownames(fitarrow), fontface = 2), data = fitarrow)+
  scale_color_viridis_c()