library(dplyr)
library(ggplot2)
library(ggthemes)

# Offense Analysis
pff<-read.csv("~/pff_grades_seas.csv")
sup_rec <- pff %>% select(PBLK,RECV,poskey = key)
sample<-readRDS('~/GitHub/mixed_effects_bootstrapping/pbp_mut.RDS')
sample$poskey <- sample$team
sample = sample %>% filter(season >= 2006)
merged <-merge(x=sample,y=sup_rec,by='poskey',how='left',no.dups = T)
#Function
return_coef <- function(formula,data,indices) {
  d <- data[indices,] 
  fit <- lm(formula, data=d)
  outp = data.frame(summary(fit)$coefficients)
  c <- outp$Estimate
  return(c)
}
#Model
fit<-lm(epa~PBLK+RECV,data=merged) %>% summary()
vars=fit$coefficients %>% rownames()
#Boot
results = boot(data=merged, statistic=return_coef,
               R=500, formula=epa~PBLK+RECV)
#Results
fr<-data.frame(results$t) 
colnames(fr) = vars
#PBLK
pblk_ci= round(quantile(fr$PBLK,probs=c(.025,.95)),4)
pblk_uci = pblk_ci[2]
pblk_lci = pblk_ci[1]
#RECV
recv_ci= round(quantile(fr$RECV,probs=c(.025,.95)),4)
recv_uci = recv_ci[2]
recv_lci = recv_ci[1]
# view results
plot <- rbind(data.frame(coef=fr$PBLK,variable='PBLK'),data.frame(coef=fr$RECV,variable='RECV'))
plot %>% ggplot(aes(x=coef,fill=variable)) +
  geom_density(alpha=.5) + 
  labs(
    title = paste('PBLK: mean=',round(mean(fr$PBLK),4),'CI=',pblk_lci,'-',pblk_uci,
                  '\nRECV: mean=',round(mean(fr$RECV),4),'CI=',recv_lci,'-',recv_uci)
  ) + theme_fivethirtyeight()

mean(fr$PBLK)/(mean(fr$PBLK) + mean(fr$RECV))
mean(fr$RECV)/(mean(fr$PBLK) + mean(fr$RECV))

# Use a .85 / .15 ratio (RECV/PBLK)

# Defensive Analysis
pff<-read.csv("~/pff_grades_seas.csv")
sup_rec <- pff %>% select(PRSH,Cov,defkey = key)
sample<-readRDS('~/GitHub/mixed_effects_bootstrapping/pbp_mut.RDS')
sample$defkey <- sample$def_team
sample = sample %>% filter(season >= 2006)
merged_def <-merge(x=sample,y=sup_rec,by='defkey',how='left',no.dups = T)
# Function
return_coef <- function(formula,data,indices) {
  d <- data[indices,] 
  fit <- lm(formula, data=d)
  outp = data.frame(summary(fit)$coefficients)
  c <- outp$Estimate
  return(c)
}
#Model
fit<-lm(epa~PRSH+Cov,data=merged_def) %>% summary()
vars=fit$coefficients %>% rownames()
#Voot
results = boot(data=merged_def, statistic=return_coef,
               R=500, formula=epa~PRSH+Cov)
#Results
fr<-data.frame(results$t) 
colnames(fr) = vars
#PBLK
prsh_ci= round(quantile(-fr$PRSH,probs=c(.025,.95)),4)
prsh_uci = -prsh_ci[2]
prsh_lci = -prsh_ci[1]
#RECV
cov_ci= round(quantile(-fr$Cov,probs=c(.025,.95)),4)
cov_uci = -cov_ci[2]
cov_lci = -cov_ci[1]
# view results
plot <- rbind(data.frame(coef=fr$PRSH,variable='PRSH'),data.frame(coef=fr$Cov,variable='Cov'))

plot %>% ggplot(aes(x=coef,fill=variable)) +
  geom_density(alpha=.5) + 
  labs(
    title = paste('PRSH: mean=',round(mean(fr$PRSH),4),'CI=',prsh_lci,'to',prsh_uci,
                  '\nCov: mean=',round(mean(fr$Cov),4),'CI=',cov_lci,'to',cov_uci)
  ) + theme_fivethirtyeight()+
  scale_x_reverse()

mean(fr$PRSH)/(mean(fr$PRSH) + mean(fr$Cov))
mean(fr$Cov)/(mean(fr$PRSH) + mean(fr$Cov))
# Use a .60 / .40 ratio (RECV/PBLK)