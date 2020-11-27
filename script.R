# Bootstraping Mixed Effects EPA impact --

setwd("~/GitHub/mixed_effects_bootstrapping")

# libraries ---------------------------------------------------------------

library(lme4)
library(dplyr)
library(plyr)
library(purrr)
library(ggplot2)
library(parallel)
#library(ggthemes)
#library(RColorBrewer)
#library(extrafont)
#loadfonts(device = "win")
#set1<-brewer.pal(n = 9, name = "Set1") 

# Read Raw Data -----------------------------------------------------------

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
    home = if_else(posteam==home_team,1,0)
    ,
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
  select(epa,passer_player_id,pos_coach,def_coach,team,wp,def_team,season,home,era)

rm(pbp)

#set.seed(20)

era1 = pbp_mut %>% filter(era == 'era1');era1=era1[sample(nrow(era1), 25000, replace = FALSE),]
era2 = pbp_mut %>% filter(era == 'era2');era2=era2[sample(nrow(era2), 25000, replace = FALSE),]
era3 = pbp_mut %>% filter(era == 'era3');era3=era3[sample(nrow(era3), 30000, replace = FALSE),]
era4 = pbp_mut %>% filter(era == 'era4');era4=era4[sample(nrow(era4), 20000, replace = FALSE),]

sample = rbind(era1,era2,era3,era4)

saveRDS(sample,'approach3/sample3.rds')

# Mixed model -------------------------------------------------------------

sample = readRDS('approach3/sample3.rds')

mixed_model<-sample %>% 
  lmer(formula=
         epa ~
         wp + 
         (1|team) + 
         (1|passer_player_id) +
         (1|team:passer_player_id)+
         (1|def_team)+
         (1|def_team:passer_player_id)+
         (1|def_team:team)
       ,
       control=lmerControl(optimizer="nloptwrap", calc.derivs = FALSE))
#Summary
mixed_model %>% summary()
getME(mixed_model, "theta")

# Sampler function --------------------------------------------------------

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


# Sampling and Clustering ----------------------------------------------------------------

memory.limit()
memory.limit(size=30000)

#set.seed(20)

#sample<-sample %>% select(epa,wind,def_coach) 

# index
start_time <- Sys.time();indx <- sampler(sample, "era", reps = 10);end_time <- Sys.time();end_time - start_time

# resample
start_time <- Sys.time(); resampled_data <- cbind(indx, sample[indx$RowID, ]);end_time <- Sys.time();end_time - start_time

#Here we are getting coefficients from original model to use as starting point
f <- fixef(mixed_model)
r <- getME(mixed_model, "theta")

rm(indx)
rm(sample)

clus <- makeCluster(4)
start_time <-Sys.time();clusterExport(clus, c("resampled_data", "f", "r"));end_time <- Sys.time();end_time - start_time
start_time <-Sys.time();clusterEvalQ(clus, require(lme4));end_time <- Sys.time();end_time - start_time

# Bootstrap function ------------------------------------------------------

boot_function <- function(i) {
  simu <- try(
    lmer(
      formula=
        epa ~
        wp + 
        (1|team) + 
        (1|passer_player_id) +
        (1|team:passer_player_id)+
        (1|def_team)+
        (1|def_team:passer_player_id)+
        (1|def_team:team)
      ,
      data = resampled_data
      ,
      subset = Replicate == i
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

# Run bootstrap -----------------------------------------------------------

rm(.Random.seed, envir=globalenv())
# min for 10 samples 
start <- proc.time(); output <- parLapplyLB(clus, X = levels(resampled_data$Replicate), fun = boot_function); end <- proc.time();end-start

#Success rate
success <- sapply(output, is.numeric)
mean(success)

#Stop cluster
stopCluster(clus)
rm(resampled_data)

# Get results -------------------------------------------------------------

final <- do.call(cbind, output[success])
final_transposed<-t(final) %>% data.frame()
coefficients <- final_transposed

saveRDS(coefficients,'approach3/bootstrap1_10samples.rds')

# Plot Result -------------------------------------------------------------

plot <- data.frame(coef = readRDS(''))

ggplot(plot,aes(x=coef))+
  geom_density(alpha=.4)+
  theme_bw()+ 
  labs(x='Coefficient', y='Density')

# cum_plot ----------------------------------------------------------------

qb <- data.frame(coef = readRDS('approach2/passer1.rds'),effect='QB')
team <- data.frame(coef = readRDS('approach2/team1.rds'),effect='Team')
coach <- data.frame(coef = readRDS('approach2/poscoach1.rds'),effect='Coach')
defteam <- data.frame(coef = readRDS('approach2/defteam1.rds'),effect='DefTeam')
defcoach <- data.frame(coef = readRDS('approach2/defcoach1.rds'),effect='DefCoach')

plot2 <- rbind(qb,team,coach,defteam,defcoach)

ggplot(plot2,aes(x=coef,fill=effect))+
  geom_density(alpha=.4)+
  theme_bw()+ 
  labs(x='Coefficient', y='Density')

# approach 3 plot ---------------------------------------------------------

plot3 <- readRDS('approach3/bootstrap1_10samples.rds')

qb<-data.frame(coef = plot$coef.passer_player_id..Intercept.,effect = 'QB')
team<- data.frame(coef = plot$coef.team..Intercept.,effect='Team')
defteam<- data.frame(coef = plot$coef.def_team..Intercept.,effect='DefTeam')

plot3 <- rbind(qb,team,defteam)

plot3$effect <- factor(plot3$effect,levels = c('DefTeam','Team','QB'))

ggplot(plot3,aes(x=coef,fill=effect))+
  geom_density(alpha=.4)+
  theme_bw()+ 
  labs(x='Coefficient', y='Density')


resampled_data
