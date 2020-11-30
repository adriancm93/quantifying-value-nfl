library(lme4)
library(dplyr)
library(plyr)
library(ggplot2)
library(ggthemes)
library(parallel)

# Read data --------------------------------------------------------
sample<-readRDS('~/GitHub/mixed_effects_bootstrapping/pbp_mut.RDS')
sample<-sample %>% select(epa,POSC_rank,PDSC_rank,passer_player_id,wp)
# Function --------------------------------------------------------
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

# Analysis --------------------------------------------------------

start_time <- Sys.time(); mixed_model<-sample %>% 
  lmer(formula=
         epa ~
         wp + 
         (1|PDSC_rank) +
         (1|PDSC_rank:POSC_rank)+
         (1|PDSC_rank:passer_player_id)
       ,
       control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE,
                             optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE))
  );end_time <- Sys.time();end_time - start_time
f <- fixef(mixed_model)
r <- getME(mixed_model, "theta")

# Model -------------------------------------------------------------------

main_start <- Sys.time()
# index
setwd("~/GitHub/mixed_effects_bootstrapping/approach3/DefTeam")
rm(.Random.seed, envir=globalenv())
start_time <- Sys.time();indx <- sampler(sample, "PDSC_rank", reps = 150);end_time <- Sys.time();(end_time - start_time)
# resample
start_time <- Sys.time(); resampled_data <- cbind(indx, sample[indx$RowID, ]);end_time <- Sys.time();print(end_time - start_time)
rm(indx)
clus <- makeCluster(4)
start_time <-Sys.time();clusterExport(clus, c("resampled_data", "f", "r"));end_time <- Sys.time();(end_time - start_time)
start_time <-Sys.time();clusterEvalQ(clus, require(lme4));end_time <- Sys.time();(end_time - start_time)
# Bootstrap function ------------------------------------------------------
boot_function <- function(i) {
  simu <- try(
    lmer(
      formula=
        epa ~
        wp + 
        (1|PDSC_rank) +
        (1|PDSC_rank:POSC_rank)+
        (1|PDSC_rank:passer_player_id)
      ,
      data = resampled_data
      ,
      subset = Replicate == i
      ,
      control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE,
                            optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE))
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
start <- Sys.time(); output <- parLapplyLB(clus, X = levels(resampled_data$Replicate), fun = boot_function); end <- Sys.time();(end-start)
#Success rate○
success <- sapply(output, is.numeric)
mean(success)
#Stop cluster
stopCluster(clus)
rm(resampled_data)
# Get results -------------------------------------------------------------
final <- do.call(cbind, output[success])
final_transposed<-t(final) %>% data.frame()
coefficients <- final_transposed
time<-stringr::str_split(Sys.time()," ")[[1]][2]
filename <- paste0('results/',stringr::str_replace_all(time,":",""),'_150.RDS')
saveRDS(object=coefficients,file=as.character(filename))
# Plot Results -------------------------------------------------------------
files <- list.files(path = "results");lst = list()
for(fil in files){
  df_i <- readRDS(paste0("results/",fil))
  print(paste(fil,df_i$PDSC_rank..Intercept. %>% mean()))
  lst[[fil]] <- df_i
}
boot <- dplyr::bind_rows(lst);defteam<-data.frame(coef = boot$PDSC_rank..Intercept.,effect = 'DefensiveRoster')
ggplot(defteam,aes(x=coef,fill=effect))+
  geom_density(alpha=.4)+theme_bw()+theme(legend.position = 'top')
print(paste('N samples =',length(defteam$coef)))
rm(boot)
main_end <- Sys.time();print(main_end - main_start)
