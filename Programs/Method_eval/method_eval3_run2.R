### Evaluating the methods (except for Hypred) in the AGAIG and reverse scenarios using the 
### simulated data with SD=1.0

library(tidyverse)
library(dplyr)
library(here)

#### sourcing the main source file with all the needed functions
source("Programs/method_eval_source.R")


##############################################################################################


##### AGAIG scenarios using SD=1.0 and estimated betas
data <- readRDS("SimData/Uganda_SimData_SD1.0.rds")
#### these betas were estimated using the sim_data_betas program
data <- predictVL(data, b0star=0.4414945, b1star=0.1170119, b2star=1.9283680, b3star=0.1369999)


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.05, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file="Results/SimResults/Results_AGAIG_SD1_ME.05_rand.rds")


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.12, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file="Results/SimResults/Results_AGAIG_SD1_ME.12_rand.rds")

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.25, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file="Results/SimResults/Results_AGAIG_SD1_ME.25_rand.rds")

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.5, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file="Results/SimResults/Results_AGAIG_SD1_ME.5_rand.rds")



##### reverse scenarios using SD=1.0 and estimated betas from reverse training set
data <- readRDS("SimData/Uganda_SimData_SD1.0.rds")
#### these estimated betas also came from the sim_data_betas program using the 'reverse' training set
data <- predictVL(data, b0star=5.8155002, b1star=-0.2442230, b2star=-3.7203965, b3star=0.1260099)


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.05, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file="Results/SimResults/Results_reverse_SD1_ME.05_rand.rds")


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.12, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file="Results/SimResults/Results_reverse_SD1_ME.12_rand.rds")

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.25, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file="Results/SimResults/Results_reverse_SD1_ME.25_rand.rds")

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.5, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file="Results/SimResults/Results_reverse_SD1_ME.5_rand.rds")
