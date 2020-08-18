### Evaluating the HyPred method in the AGAIG and reverse scenarios using the 
### simulated data with SD=1.0 and SD=0

library(tidyverse)
library(dplyr)
library(here)

#### sourcing the main source file with all the needed functions
source("Programs/method_eval_source.R")

##############################################################################################


##### AGAIG scenarios using SD=1.0 and estimated betas
data <- readRDS("SimData/Records_used_in_500_matrices_hypred/SD1.0_data_500_run2.rds")
data <- predictVL(data, b0star=0.4414945, b1star=0.1170119, b2star=1.9283680, b3star=0.1369999)



set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=0, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_AGAIG_SD1_ME0_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.025, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_AGAIG_SD1_ME.025_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.05, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_AGAIG_SD1_ME.05_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.075, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_AGAIG_SD1_ME.075_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.1, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_AGAIG_SD1_ME.1_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.125, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_AGAIG_SD1_ME.125_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.15, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_AGAIG_SD1_ME.15_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.175, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_AGAIG_SD1_ME.175_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.2, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_AGAIG_SD1_ME.2_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.225, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_AGAIG_SD1_ME.225_rand.rds")



set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.25, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_AGAIG_SD1_ME.25_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.5, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_AGAIG_SD1_ME.5_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.75, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_AGAIG_SD1_ME.75_rand.rds")


##################################################################################################################################

##### AGAIG scenarios using SD=0 and estimated betas
data <- readRDS("SimData/Records_used_in_500_matrices_hypred/SD0_data_500_run2.rds")
data <- predictVL(data, b0star=0.4414945, b1star=0.1170119, b2star=1.9283680, b3star=0.1369999)


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=0, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_AGAIG_SD0_ME0_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.025, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_AGAIG_SD0_ME.025_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.05, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_AGAIG_SD0_ME.05_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.075, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_AGAIG_SD0_ME.075_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.1, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_AGAIG_SD0_ME.1_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.125, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_AGAIG_SD0_ME.125_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.15, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_AGAIG_SD0_ME.15_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.175, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_AGAIG_SD0_ME.175_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.2, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_AGAIG_SD0_ME.2_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.225, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_AGAIG_SD0_ME.225_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.25, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_AGAIG_SD0_ME.25_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.5, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_AGAIG_SD0_ME.5_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.75, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_AGAIG_SD0_ME.75_rand.rds")

######################################################################################################################################

##### Reverse scenarios using SD=1.0 and estimated betas
data <- readRDS("SimData/Records_used_in_500_matrices_hypred/SD1.0_data_500_run2.rds")
data <- predictVL(data, b0star=5.8155002, b1star=-0.2442230, b2star=-3.7203965, b3star=0.1260099)


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=0, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_reverse_SD1_ME0_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.025, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_reverse_SD1_ME.025_rand.rds")



set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.05, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_reverse_SD1_ME.05_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.075, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_reverse_SD1_ME.075_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.1, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_reverse_SD1_ME.1_rand.rds")



set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.125, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_reverse_SD1_ME.125_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.15, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_reverse_SD1_ME.15_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.175, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_reverse_SD1_ME.175_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.2, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_reverse_SD1_ME.2_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.225, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_reverse_SD1_ME.225_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.25, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_reverse_SD1_ME.25_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.5, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_reverse_SD1_ME.5_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.75, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_reverse_SD1_ME.75_rand.rds")


##################################################################################################################################

##### Reverse scenarios using SD=0 and estimated betas
data <- readRDS("SimData/Records_used_in_500_matrices_hypred/SD0_data_500_run2.rds")
data <- predictVL(data, b0star=5.8155002, b1star=-0.2442230, b2star=-3.7203965, b3star=0.1260099)



set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=0, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_reverse_SD0_ME0_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.025, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_reverse_SD0_ME.025_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.05, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_reverse_SD0_ME.05_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.075, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_reverse_SD0_ME.075_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.1, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_reverse_SD0_ME.1_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.125, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_reverse_SD0_ME.125_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.15, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_reverse_SD0_ME.15_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.175, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_reverse_SD0_ME.175_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.2, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_reverse_SD0_ME.2_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.225, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_reverse_SD0_ME.225_rand.rds")



set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.25, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_reverse_SD0_ME.25_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.5, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_reverse_SD0_ME.5_rand.rds")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.75, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


saveRDS(result1, file="Results/SimResults/Results_Hypred_reverse_SD0_ME.75_rand.rds")

######################################################################################################################################

