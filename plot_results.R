library(dplyr)
library(ggplot2)
library(ggthemes)

# prepare -----------------------------------------------------------------
#QB
files <- list.files(path = "~/GitHub/quantifying-value-nfl/QB/results");lst = list()
for(fil in files){
  df_i <- readRDS(paste0("C:/Users/adrian-boss/Documents/GitHub/quantifying-value-nfl/QB/results/",fil))
  print(paste(fil,df_i$passer_player_id..Intercept. %>% mean()))
  lst[[fil]] <- df_i
}
qb <- dplyr::bind_rows(lst);qb<-data.frame(coef = qb$passer_player_id..Intercept.,effect = 'QB')
qb_lst = lst()
for (i in 1:300){
  m <- mean(sample(qb$coef,50,replace=FALSE))
  qb_lst[[i]] <- m
}; 
qb_mc<-data.frame(coef=qb_lst %>% unlist(),effect = 'QB')

#DefTeam
files <- list.files(path = "~/GitHub/quantifying-value-nfl/DefTeam/results");lst = list()
for(fil in files){
  df_i <- readRDS(paste0("~/GitHub/quantifying-value-nfl/DefTeam/results/",fil))
  print(paste(fil,df_i$PDSC_rank..Intercept. %>% mean()))
  lst[[fil]] <- df_i
}
defteam <- dplyr::bind_rows(lst);defteam<-data.frame(coef = defteam$PDSC_rank..Intercept.,effect = 'Defensive Roster')
defteam_lst = lst()
for (i in 1:300){
  m <- mean(sample(defteam$coef,2,replace=FALSE))
  defteam_lst[[i]] <- m
}; 
defteam_mc<-data.frame(coef=defteam_lst %>% unlist(),effect = 'Defensive Roster')

#Team
files <- list.files(path = "~/GitHub/quantifying-value-nfl/Team/results");lst = list()
for(fil in files){
  df_i <- readRDS(paste0("~/GitHub/quantifying-value-nfl/Team/results/",fil))
  print(paste(fil,df_i$POSC_rank..Intercept. %>% mean()))
  lst[[fil]] <- df_i
}
team <- dplyr::bind_rows(lst);team<-data.frame(coef = team$POSC_rank..Intercept.,effect = 'Offensive Supporting Cast')
team_lst = lst()
for (i in 1:300){
  m <- mean(sample(team$coef,2,replace=FALSE))
  team_lst[[i]] <- m
}; 
team_mc<-data.frame(coef=team_lst %>% unlist(),effect = 'Offensive Supporting Cast')
#Def Coach
files <- list.files(path = "~/GitHub/quantifying-value-nfl/DefCoach/results");lst = list()
for(fil in files){
  df_i <- readRDS(paste0("~/GitHub/quantifying-value-nfl/DefCoach/results/",fil))
  print(paste(fil,df_i$def_coach..Intercept. %>% mean()))
  lst[[fil]] <- df_i
}
defcoach <- dplyr::bind_rows(lst);defcoach<-data.frame(coef = defcoach$def_coach..Intercept.,effect = 'Defensive Coach')
defcoach_lst = lst()
for (i in 1:300){
  m <- mean(sample(defcoach$coef,8,replace=FALSE))
  defcoach_lst[[i]] <- m
}; 
defcoach_mc<-data.frame(coef=defcoach_lst %>% unlist(),effect = 'Defensive Coach')
#PosCoach
files <- list.files(path = "~/GitHub/quantifying-value-nfl/PosCoach/results");lst = list()
for(fil in files){
  df_i <- readRDS(paste0("~/GitHub/quantifying-value-nfl/PosCoach/results/",fil))
  print(paste(fil,df_i$pos_coach..Intercept. %>% mean()))
  lst[[fil]] <- df_i
}
poscoach <- dplyr::bind_rows(lst);poscoach<-data.frame(coef = poscoach$pos_coach..Intercept.,effect = 'Offensive Coach')
poscoach_lst = lst()
for (i in 1:300){
  m <- mean(sample(poscoach$coef,2,replace=FALSE))
  poscoach_lst[[i]] <- m
}; 
poscoach_mc<-data.frame(coef=poscoach_lst %>% unlist(),effect = 'Offensive Coach')



# Plot Results ------------------------------------------------------------

plot<-readRDS('~/GitHub/quantifying-value-nfl/plot_data.RDS')
plot$effect <- factor(plot$effect,levels = c('Offensive Coach',
                                             'Defensive Coach',
                                             'Defensive Roster',
                                             'Offensive Supporting Cast',
                                             'QB'))
setwd("~/GitHub/quantifying-value-nfl")

plot%>% ggplot(aes(x=coef,fill=effect))+
  geom_density(alpha=.4)+
  scale_fill_manual(
    values=c("red","#377EB8","#4DAF4A","#FF7F00","#984EA3")) +
  annotate(geom = "label", x = 0, y = 250, hjust = "left",fill="white",vjust = 1,
           label = "Extremely hard for a HC to have passing game successwithout an elite QB \nor a good QB + supporting cast combination",
           color = "black") +
  labs(
    title = 'Measuring Impact',
    subtitle = 'QB is still the most impactful; Defensive Coaches tend to matter more than Offensive',
    caption = 'By: Adrian Cadena @adrian_stats',
    y = 'Density',
    x = 'Impact on Passing Efficiency'
  )+
  theme_bw()+
  theme(
    plot.title = element_text(hjust=0.5,size=18),
    plot.subtitle = element_text(hjust=0.5,size=13),
    legend.position = 'top',
    legend.title = element_blank(),
    legend.text = element_text(size=10),
    panel.grid.minor = element_blank() )+
  coord_cartesian(ylim=c(0,250)) + ggsave('out.png',dpi=600, width = 8, height = 6,units = 'in')






