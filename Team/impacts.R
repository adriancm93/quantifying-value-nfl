library(dplyr)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
setwd("~/GitHub/mixed_effects_bootstrapping/Team")

files <- list.files(path = "results");lst = list()
for(fil in files){
  df_i <- readRDS(paste0("results/",fil))
  lst[[fil]] <- df_i
}
results <- dplyr::bind_rows(lst)

lst <- list()
for (m in colnames(results)[3:4]){
  lst[[m]] <- data.frame(coef=results[,m],effect = m) 
};plot<-dplyr::bind_rows(lst)
plot$effect <- factor(plot$effect,levels = c('POSC_rank.PDSC_rank..Intercept.','POSC_rank.passer_player_id..Intercept.'))

plot %>% 
  ggplot(aes(x=coef,fill=effect)) +
  geom_density(alpha=.5) +
  scale_fill_manual(values=RColorBrewer::brewer.pal(name='Set3',n=8)[5:6])+
  labs(
    title= "Offensive supporting casts are both helped by their QB \nand impacted by opposing defenses at a similar rate"
  ) + theme_bw() + 
  theme(
    legend.title = element_blank(),
    legend.position = 'top',
    plot.title = element_text(hjust=0.5)
  ) 
