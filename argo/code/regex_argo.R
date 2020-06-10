#get numbers of float with hplc from url
url <- "http://www.oao.obs-vlfr.fr/DEPLOYMENT/Affiche_Ancillary_IS.php"
page <- readLines(con = url) %>% 
  as.data.frame()
names(page) <- "ligne_html"
page$ligne_html <- as.character(page$ligne_html)

page2 <- separate(page, ligne_html, sep = "<tr>", into = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17", "18","19","20","21","22","23"))
chain = c()
for (i in c(1:23)){
  t <- page2[,i]
  chain <- c(chain,t)
}
chain <- na.omit(chain)
chain_df <- data.frame("chain" = chain)
chain_df2 <- separate(chain_df, chain, sep = "<td>", into = c("nimp", "float", "ctd", "pigments", "nitrate", "doxy", "in_air", "residuals"), extra = "merge")

float_we_need <- chain_df2[grep("x", chain_df2$pigments),]

page <- data.frame("chain" = chain_df2$residuals)
page <- na.omit(page)
names(page) <- "ligne_html"

page2 <- separate(page, ligne_html, sep = "<tr>", into = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17", "18","19","20","21","22","23"))
chain = c()
for (i in c(1:23)){
  t <- page2[,i]
  chain <- c(chain,t)
}
chain <- na.omit(chain)
chain_df <- data.frame("chain" = chain)
chain_df2 <- separate(chain_df, chain, sep = "<td>", into = c("nimp", "float", "ctd", "pigments", "nitrate", "doxy", "in_air", "residuals"), extra = "merge")

float_we_need2 <- chain_df2[grep("x", chain_df2$ctd),]

float_we_need <- bind_rows(float_we_need,float_we_need2)


float_we_need$number <- ifelse(as.numeric(rownames(float_we_need))<40, sub("^.+,", "", float_we_need$float),sub("^.+,", "", float_we_need$nimp) )
float_we_need$number <- as.numeric(gsub("([0-9]+).*$", "\\1",float_we_need$number))

float_we_need$lovbio <- ifelse(as.numeric(rownames(float_we_need))<40, gsub("^.+[I][T][U]", "", float_we_need$float), gsub("^.+[I][T][U]", "", float_we_need$nimp))
float_we_need$lovbio <- substr(float_we_need$lovbio,2,11)
    
float_ref <- select(float_we_need, number, lovbio)
write_csv(float_ref, "Scripts/Data/argo/ref.csv")

unique(float_we_need$lovbio)
#download associate floats witj mobaxterm
#with the fiollowing code : "for i in number of my profiles 
# do wget -N —user=anonymous --password='your@adress' ftp://ftp.ifremer.fr/ifremer/argo/dac/coriolis/$i/profiles/MR*
#done

#new download for descent profiles with this line 
#do wget -N —user=anonymous --password='your@adress' ftp://ftp.ifremer.fr/ifremer/argo/dac/coriolis/$i/profiles/BD${i}_001D.nc
