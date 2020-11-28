setwd("~/GitHub/mixed_effects_bootstrapping/approach1")

library(lme4)
library(dplyr)
library(plyr)
library(purrr)
library(ggplot2)
library(parallel)

sample = readRDS('sample1.rds')

mixed_model<-sample %>% 
  lmer(formula=
         epa ~
         wp + 
         (1|team)+
         (1|pos_coach)+
         (1|pos_coach:team)+
         (1|pos_coach:passer_player_id)+
         (1|pos_coach:def_team)
       ,
       control=lmerControl(optimizer="nloptwrap", calc.derivs = FALSE))


#Summary
mixed_model %>% summary()
getME(mixed_model, "theta")

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

sample<-sample %>% select(epa,wp,team,def_team,passer_player_id,era,pos_coach) 

rm(.Random.seed, envir=globalenv())

# index
start_time <- Sys.time();indx <- sampler(sample, "era", reps = 5);end_time <- Sys.time();end_time - start_time
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
        (1|team)+
        (1|pos_coach)+
        (1|pos_coach:team)+
        (1|pos_coach:passer_player_id)+
        (1|pos_coach:def_team)
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

ef1 <- data.frame('coef'=coefficients$pos_coach..Intercept.,'vari'='coach') 
ef2 <- data.frame('coef'=coefficients$team..Intercept.,'vari'='team') 
p <- rbind(ef1,ef2)

ggplot(p,aes(x=coef,fill=vari))+
  geom_density(alpha=.4)









