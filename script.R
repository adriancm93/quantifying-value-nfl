#Libraries 
library(lme4)
library(dplyr)
library(plyr)
library(purrr)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(extrafont)
loadfonts(device = "win")
library(parallel)
set1<-brewer.pal(n = 9, name = "Set1") 

seasons <- 2000:2019
pbp <- purrr::map_df(seasons, function(x) {
  readRDS(
    url(
      glue::glue("https://raw.githubusercontent.com/guga31bb/nflfastR-data/master/data/play_by_play_{x}.rds")
    )
  )
})

# data prep ---------------------------------------------------------------


npass <- pbp %>%
  filter(play_type == 'pass',
         season_type == 'REG') %>%
  group_by(passer_player_id) %>% 
  dplyr::summarise(
    num_plays = n()
  ); pbp<-merge(pbp,npass,by='passer_player_id',all.x = T,no.dups = T)

#nflfastR's EPA model already accounts for a lot of variables so we will keep it simple
pbp_mut<-pbp%>% 
  dplyr::filter(
    play_type == 'pass',
    season_type == 'REG',
    wp >= .20,
    wp <= .80,
    game_half != 'Overtime',
    !is.na(epa),
    num_plays >= 100,
    !is.na(down),
    penalty == 0
  )%>%
  dplyr::mutate(
    pos_coach = if_else(posteam==home_team,home_coach,away_coach) 
    ,
    def_coach = if_else(defteam==home_team,home_coach,away_coach) 
    ,
    wind = ifelse(is.na(wind)==T,min(wind,na.rm=T),wind)
    ,
    temp = ifelse(is.na(temp)==T,75,temp)
    ,
    team = paste(season,posteam,sep='_')
    ,
    def_team = paste(season,defteam,sep='_')
    ,
    outdoor = if_else(roof %in% c('outdoors','open'),1,0)
    ,
    era = if_else(season %in% 2000:2005,'era1',
                  if_else(season %in% 2006:2013,'era2',
                          if_else(season %in% 2014:2017,'era3','era4')))
  )   %>% 
  select(epa,temp,wind,passer_player_id,pos_coach,def_coach,team,wp,def_team,season,outdoor)

sample = pbp_mut[sample(nrow(pbp_mut), 50000, replace = FALSE),]

saveRDS(sample,'sample.rds')



# analysis ----------------------------------------------------------------

#Read Data
sample = readRDS('sample.rds')

#Mixed model 
mixed_model<-sample %>% 
  lmer(formula=
         epa ~
         wp + 
         (1|passer_player_id)+
         (1|pos_coach)+
         (1|team)+
         (1|def_team)+
         (1|def_coach)
       ,
       control=lmerControl(optimizer="nloptwrap", calc.derivs = FALSE))
#Summary
mixed_model %>% summary()

#This is what we are going to be getting for each bootstrap simulation
getME(mixed_model, "theta") #Research Cholesky factors | try ?getME

getME(mixed_model, "theta")['passer_player_id.(Intercept)']

#Bootstrap 
#Followed this tutorial: https://stats.idre.ucla.edu/r/dae/mixed-effects-logistic-regression/

#Create sampler function. 
#Source: http://biostat.mc.vanderbilt.edu/wiki/Main/HowToBootstrapCorrelatedData

sampler <- function(dat, clustervar, replace = TRUE, reps = 1) {
  cid <- unique(dat[, clustervar[1]])
  ncid <- length(cid)
  recid <- sample(cid, size = ncid * reps, replace = TRUE)
  if (replace) {
    rid <- lapply(seq_along(recid), function(i) {
      cbind(NewID = i, RowID = sample(which(dat[, clustervar] == recid[i]),
                                      size = length(which(dat[, clustervar] == recid[i])), replace = TRUE))
    })
  } else {
    rid <- lapply(seq_along(recid), function(i) {
      cbind(NewID = i, RowID = which(dat[, clustervar] == recid[i]))
    })
  }
  dat <- as.data.frame(do.call(rbind, rid))
  dat$Replicate <- factor(cut(dat$NewID, breaks = c(1, ncid * 1:reps), include.lowest = TRUE,
                              labels = FALSE))
  dat$NewID <- factor(dat$NewID)
  return(dat)
}

memory.limit()
memory.limit(size=30000)

#Set seed
set.seed(20)

#index
start_time <- Sys.time();indx <- sampler(sample, "passer_player_id", reps = 300);end_time <- Sys.time();end_time - start_time
#resampled
start_time <- Sys.time(); resampled_data <- cbind(indx, sample[indx$RowID, ]);end_time <- Sys.time();end_time - start_time

#Here we are getting coefficients from original model to use as starting point
f <- fixef(mixed_model)
r <- getME(mixed_model, "theta")

rm(indx)
rm(sample)
clus <- makeCluster(4)
start_time <-Sys.time();clusterExport(clus, c("resampled_data", "f", "r"));end_time <- Sys.time();end_time - start_time
start_time <-Sys.time();clusterEvalQ(clus, require(lme4));end_time <- Sys.time();end_time - start_time

#Create Bootstrapping  function
boot_function <- function(i) {
  simu <- try(
    lmer(
      formula=
        epa ~
        wp + 
        (1|passer_player_id)+
        (1|pos_coach)+
        (1|team)+
        (1|def_team)+
        (1|def_coach)
      ,
      data = resampled_data
      ,
      control=lmerControl(optimizer="nloptwrap", calc.derivs = FALSE)
      ,
      start = list(fixef = f, theta = r)
    )
    ,
    silent = TRUE)
  if (class(simu) == "try-error")
    return(simu)
  c(fixef(simu), getME(simu, "theta"))
}

#Run bootstrap (this will take a while) # 1 hours 21 minutes in my laptop
set.seed(20)
start <- proc.time(); output <- parLapplyLB(clus, X = levels(resampled_data$Replicate), fun = boot_function); end <- proc.time()
end-start

#Stop cluster
stopCluster(clus)

#I recommend to save resampled data
#saveRDS(resampled_data,'name.RDS')

#Delete heavy data (save first?)
#rm(resampled_data)

#Success rate
success <- sapply(output, is.numeric)
mean(success)

#Create dataframe
final <- do.call(cbind, output[success])
final_transposed<-t(final) %>% data.frame()
saveRDS(final_transposed,'final_transposed3.RDS')


#Prepare data for plot 
coach <- final_transposed %>% select(pos_coach..Intercept.) %>% mutate(Variable='Head Coach') %>% select(intercept = pos_coach..Intercept.,Variable)
qb <- final_transposed %>% select(passer_player_id..Intercept.) %>% mutate(Variable='Quarterback') %>% select(intercept = passer_player_id..Intercept.,Variable)
#def_team <- final_transposed %>% select(opponent..Intercept.) %>% mutate(Variable='Opposing Defensive Team')%>% select(intercept = opponent..Intercept.,Variable)
#def_coach <- final_transposed %>% select(def_coach..Intercept.) %>% mutate(Variable='Opposing Head Coach')%>% select(intercept = opponent..Intercept.,Variable)
team <- final_transposed %>% select(team..Intercept.) %>% mutate(Variable='Team (Proxy for Offensive Supporting Cast)')%>% select(intercept = team..Intercept.,Variable)

#(I decided not to include def_team. I feel more resampling is needed for an accurate estimate (around 1000): too much variance in results)
plot <- rbind(coach,qb,team)
plot$Variable <- factor(plot$Variable, levels = c('Head Coach','Team (Proxy for Offensive Supporting Cast)','Quarterback'))

#Straight forward plot
ggplot(plot,aes(x=intercept, fill=Variable))+
  geom_density(alpha=.4)+
  theme_fivethirtyeight()+ 
  scale_fill_manual(values=c("#4DAF4A","#FF7F00","#984EA3")) + 
  labs(x='Absolute Impact on Passing Offense Efficiency (EPA) | Further to the right means greater impact (positive or negative)', y='Density',
       title = "Quarterback-Ability Matters More Than Supporting Offensive Talent and Head Coach Competence",
       subtitle = 'When predicting passing efficiency, quarterback ability matters the most, followed by offensive personnel',
       caption = "Data and EPA models: nflfastR | Seasons 2006 - 2019 | Win Probability 20%-80% | Quarters 1-4 | Passing Plays Only
Mixed Effects Model to Find Intercept | Using 200 bootstrapped (re-sampled) datasets to find the distribution of Cholesky Factorized Coefficients"
  )





#This is the code I ised to add notes and prepare the plot to look nice for Twitter. 
#It will look pretty weird in Rstudio, but  will look good when you save it (using my specs) and open the file
ggplot(plot,aes(x=intercept, fill=Variable))+
  geom_density(alpha=.4,size=.15)+
  theme_fivethirtyeight()+ 
  scale_fill_manual(values=c("#377EB8","#4DAF4A","#FF7F00","#984EA3"))+
  labs(x='Absolute Impact on Passing Offense Efficiency (EPA) | Further to the right means greater impact (positive or negative)', y='Density',
       title = "Quarterback-Ability Matters More Than Supporting Offensive Talent and Head Coach Competence",
       subtitle = 'When predicting passing efficiency, quarterback ability matters the most, followed by offensive personnel',
       caption = "Data and EPA models: nflfastR | Seasons 2006 - 2019 | Win Probability 20%-80% | Quarters 1-4 | Passing Plays Only
Mixed Effects Model to Find Intercept | Using 200 bootstrapped (re-sampled) datasets to find the distribution of Cholesky Factorized Coefficients"
  )+
  theme(text = element_text(),
        plot.title = element_text(size = 7, family = "Trebuchet MS",color = "grey20",hjust = .5),
        plot.subtitle = element_text(size = 6, family = "Trebuchet MS",color = "grey20",hjust = .5),
        axis.title = element_text(size = 5, family = "Trebuchet MS",color = "grey20"),
        axis.text = element_text(size = 3, family = "Trebuchet MS",color = "grey20"),
        legend.text = element_text(size = 5, family = "Trebuchet MS",color = "grey20"),
        legend.title = element_blank(),
        legend.position="top",
        legend.key.size =unit(.3,"line"),
        panel.grid = element_line(size=.11),
        plot.caption = element_text(size=4, family = "Trebuchet MS",color = "grey20",hjust = 0),
        plot.caption.position = "panel",
  )  +
  annotate(geom = "label", x = .001, y = 135, hjust = "left",fill="#F0F0F0",vjust = 1,
           label = "It is extremely hard for a HC to have passing game success 
without an elite QB or a good QB + supporting cast combination",
           size=1.6, 
           family = "Trebuchet MS",
           color = "grey20")  + 
  annotate(
    geom = "curve", x = .11, y = 50,xend = .103, yend = 16, 
    curvature = -.2, size= unit(.12, "mm"),arrow = arrow(length = unit(1.2, "mm")),
  )+
  annotate(geom = "label", x = .095, y = 65, hjust = "left",size=1.45,family = "Trebuchet MS",color = "grey20",fill="#F0F0F0",vjust = 1,
           label = 'Quarterback ability has 
the biggest impact, variance, 
and upside') + 
  annotate(
    geom = "curve", x = .102, y = 101,xend = .078, yend = 67, 
    curvature = -.2, size= unit(.15, "mm"),arrow = arrow(length = unit(1.2, "mm")),
  )+
  annotate(geom = "text", x = .101, y = 115, hjust = "left",size=1.45,family = "Trebuchet MS",color = "grey20",
           label = 'A good supporting cast 
greatly improves passing 
efficiency')+ 
  annotate(
    geom = "curve", x = .0505, y = 90,xend = .043, yend = 44, 
    curvature = .1, size= unit(.15, "mm"),arrow = arrow(length = unit(1.2, "mm")),
  )+
  annotate(geom = "text", x = .051, y = 92, hjust = "left",size=1.45,family = "Trebuchet MS",color = "grey20",
           label = 'Very rarely a HC matters 
more than the QB')+ 
  annotate(
    geom = "curve", x = .014, y = 67,xend = .021, yend = 49, 
    curvature = .1, size= unit(.15, "mm"),arrow = arrow(length = unit(1.2, "mm")),
  )+
  annotate(geom = "text", x = .001, y = 77, hjust = "left",size=1.45,family = "Trebuchet MS",color = "grey20",
           label = 'Defending HC is the least 
impactful factor') 

ggsave('passing_factors.png', dpi=1000, width = 12.5, height = 8, units = "cm")