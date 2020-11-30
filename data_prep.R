#Libraries
library(boot)
library(dplyr)
library(ggplot2)
library(ggthemes)
# Read and mutate data ----------------------------------------------------
seasons <- 2006:2019
pbp <- purrr::map_df(seasons, function(x) {
  readRDS(
    url(
      glue::glue("https://raw.githubusercontent.com/guga31bb/nflfastR-data/master/data/play_by_play_{x}.rds")
    )
  )
})
npass <- pbp %>%
  filter(play_type == 'pass',
         season_type == 'REG') %>%
  group_by(passer_player_id) %>% 
  dplyr::summarise(
    num_plays = n()
  ); pbp<-merge(pbp,npass,by='passer_player_id',all.x = T,no.dups = T)
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
    team = paste(season,posteam,sep='_')
    ,
    def_team = paste(season,defteam,sep='_')
  )   %>% 
  select(epa,passer_player_id,pos_coach,def_coach,team,wp,def_team,season)
saveRDS(pbp_mut,'~/GitHub/quantifying-value-nfl/pbp_mut.RDS')
# new Features -------------------------------------------------------------
#Offensive
pff<-read.csv("~/pff_grades_seas.csv")
pff['POSC'] <- pff['RECV']*.85 + pff['PBLK']*.15
pff['PDSC'] <- pff['Cov']*.60 + pff['PRSH']*.40
pff<-pff  %>%
  arrange(Year,desc(POSC))
rank <- pff %>% select(Year,Team_abr,POSC)%>%
  arrange(Year,desc(POSC)) %>% group_by(Year) %>%
  summarise(
    POSC_rank = 1:32) %>% ungroup()
pff<-cbind(pff,rank %>% select(POSC_rank))
sup_pass <- pff %>% select(POSC,POSC_rank,team = key)
#Defensive
pff<-pff  %>%
  arrange(Year,desc(PDSC))
rank <- pff %>% select(Year,Team_abr,PDSC)%>%
  arrange(Year,desc(PDSC)) %>% group_by(Year) %>%
  summarise(
    PDSC_rank = 1:32) %>% ungroup()
pff<-cbind(pff,rank %>% select(PDSC_rank))
def_pass <- pff %>% select(PDSC,PDSC_rank,def_team = key)
#Merge
pbp_mut<-readRDS('~/GitHub/quantifying-value-nfl/pbp_mut.RDS')
pbp_mut <- pbp_mut %>% filter(season >= 2006)
pbp_mut=merge(x=pbp_mut,y=sup_pass,by='team',how='left',no.dups = T)
pbp_mut=merge(x=pbp_mut,y=def_pass,by='def_team',how='left',no.dups = T)
saveRDS(pbp_mut,'~/GitHub/quantifying-value-nfl/pbp_mut.RDS')
