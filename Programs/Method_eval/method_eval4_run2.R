#==============================================================================
# FILENAME: method_eval4_run2.R
# PROJECT: 	Pooled testing in HIV
# PURPOSE: This program calls the main evaluation function to evaluate the pooled testing methods
#          in different scenarios, 4 different scenarios each evaluated in 13 different MEs using 
#          the SD=0 data; does not evaluate the hypred method, that is a separate program
#          
#          
# AUTHOR: Adam Brand


# INPUT datasets: Uganda_SimData_SD0.rds
#                 SD0_data_500_run2_perm.rds

# OUTPUT datasets: 52 results table located here: Pooled_Testing_HIV/Results/SimResults
#                  13 file names beginning with Results_AGAIG_SD0
#                  13 file names beginning with Results_reverse_SD0
#                  13 file names beginning with Results_noassoc_SD0
#                  13 file names beginning with Results_misspec_SD0


# R VERSION: 3.6.1
#==============================================================================
#Notes: 





# =============================================================================


### Evaluating the methods (except for Hypred) in the AGAIG and reverse scenarios using the 
### simulated data with SD=0

library(tidyverse)
library(dplyr)
library(here)

#### sourcing the main source file with all the needed functions
source("Programs/method_eval_source.R")
##############################################################################################


##### AGAIG scenarios using SD=0 and estimated betas
data <- readRDS("SimData/Uganda_SimData_SD0.rds")
#### these betas were estimated using the sim_data_betas program
#### the betas are estimated using the training set with SD=1.0. Using the training set with SD=0
#### provided 'too accurate' results due to the lack of variation
data <- predictVL(data, b0star=0.4414945, b1star=0.1170119, b2star=1.9283680, b3star=0.1369999)



set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=0, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_AGAIG_SD0_ME0_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.025, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_AGAIG_SD0_ME.025_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.05, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_AGAIG_SD0_ME.05_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.075, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_AGAIG_SD0_ME.075_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.1, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_AGAIG_SD0_ME.1_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.125, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_AGAIG_SD0_ME.125_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.15, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_AGAIG_SD0_ME.15_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.175, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_AGAIG_SD0_ME.175_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.2, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_AGAIG_SD0_ME.2_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.225, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_AGAIG_SD0_ME.225_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.25, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_AGAIG_SD0_ME.25_rand.rds"))

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.5, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_AGAIG_SD0_ME.5_rand.rds"))

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.75, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_AGAIG_SD0_ME.75_rand.rds"))



##### reverse scenarios using SD=0 and estimated betas from the reverse training set
data <- readRDS("SimData/Uganda_SimData_SD0.rds")
data <- predictVL(data, b0star=5.8155002, b1star=-0.2442230, b2star=-3.7203965, b3star=0.1260099)


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=0, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_reverse_SD0_ME0_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.025, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_reverse_SD0_ME.025_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.05, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_reverse_SD0_ME.05_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.075, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_reverse_SD0_ME.075_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.1, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_reverse_SD0_ME.1_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.125, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_reverse_SD0_ME.125_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.15, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_reverse_SD0_ME.15_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.175, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_reverse_SD0_ME.175_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.2, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_reverse_SD0_ME.2_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.225, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_reverse_SD0_ME.225_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.25, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_reverse_SD0_ME.25_rand.rds"))

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.5, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_reverse_SD0_ME.5_rand.rds"))

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.75, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_reverse_SD0_ME.75_rand.rds"))


##### no prediction association scenarios using SD=0 and estimated betas
data <- readRDS("SimData/Records_used_in_500_matrices_hypred/SD0_data_500_run2_perm.rds")
#### these betas were estimated using the sim_data_betas program
#### the betas are estimated using the training set with SD=1.0. Using the training set with SD=0
#### provided 'too accurate' results due to the lack of variation
data <- predictVL(data, b0star=0.4414945, b1star=0.1170119, b2star=1.9283680, b3star=0.1369999)



set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=0, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_noassoc_SD0_ME0_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.025, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_noassoc_SD0_ME.025_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.05, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_noassoc_SD0_ME.05_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.075, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_noassoc_SD0_ME.075_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.1, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_noassoc_SD0_ME.1_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.125, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_noassoc_SD0_ME.125_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.15, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_noassoc_SD0_ME.15_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.175, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_noassoc_SD0_ME.175_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.2, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_noassoc_SD0_ME.2_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.225, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_noassoc_SD0_ME.225_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.25, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_noassoc_SD0_ME.25_rand.rds"))

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.5, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_noassoc_SD0_ME.5_rand.rds"))

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.75, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_noassoc_SD0_ME.75_rand.rds"))




##### misspecified model scenarios using SD=0 and estimated betas
data <- readRDS("SimData/Uganda_SimData_SD0.rds")
#### these betas were estimated using the sim_data_betas program
#### the betas are estimated using the training set with SD=1.0 without the binary predictor. Using the training set with SD=0
#### provided 'too accurate' results due to the lack of variation
data <- predictVL(data, b0star=0.623, b1star=0.133, b2star=0, b3star=0)



set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=0, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_misspec_SD0_ME0_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.025, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_misspec_SD0_ME.025_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.05, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_misspec_SD0_ME.05_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.075, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_misspec_SD0_ME.075_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.1, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_misspec_SD0_ME.1_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.125, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_misspec_SD0_ME.125_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.15, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_misspec_SD0_ME.15_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.175, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_misspec_SD0_ME.175_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.2, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_misspec_SD0_ME.2_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.225, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_misspec_SD0_ME.225_rand.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.25, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_misspec_SD0_ME.25_rand.rds"))

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.5, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_misspec_SD0_ME.5_rand.rds"))

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.75, tstperd=5, lowlimit=50, filltyp="rnd")


saveRDS(result1, file=here("Results","SimResults","Results_misspec_SD0_ME.75_rand.rds"))






