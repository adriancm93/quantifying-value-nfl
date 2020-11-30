library(dplyr)
library(ggplot2)
library(ggthemes)
setwd("~/GitHub/mixed_effects_bootstrapping/qb")

files <- list.files(path = "results");lst = list()
for(fil in files){
  df_i <- readRDS(paste0("results/",fil))
  print(paste(fil,df_i$passer_player_id..Intercept. %>% mean()))
  lst[[fil]] <- df_i
}
results <- dplyr::bind_rows(lst)

lst <- list()
for (m in colnames(results)[3:4]){
  lst[[m]] <- data.frame(coef=results[,m],effect = m) 
  };plot<-dplyr::bind_rows(lst)
plot$effect <- factor(plot$effect,levels = c('passer_player_id.POSC_rank..Intercept.','passer_player_id.PDSC_rank..Intercept.'))

plot %>% 
  ggplot(aes(x=coef,fill=effect)) +
  labs(
    title= "Opposing Defenses Impact QB's Produciton more than Offensive Supporting Cast"
  )+
  geom_density(alpha=.5) + theme_bw() + 
  theme(
    legend.title = element_blank(),
    legend.position = 'top',
    plot.title = element_text(hjust=0.5)
  )  






