#==============================================================================
# FILENAME: mwthod_eval3_run2.R
# PROJECT: 	Pooled testing in HIV
# PURPOSE: evaluating ALL pooled testing methods on the real data from Uganda with 13 different MEs
#          
#          
# AUTHOR: Adam Brand


# INPUT datasets: test_set_final.rds

# OUTPUT dataset: All 26 result tables located here Pooled_Testing_HIV\Results\UgandaResults
#                 13 filenames begin with Uganda_hypred_ME
#                 13 filenames begin with Uganda_ME


# R VERSION: 3.6.1
#==============================================================================
#Notes: 





# =============================================================================




## Evaluating the pooled testing methods using the real Uganda data
## evaluating the Hypred method with the first statements and the other methods below.
## As the real data are the observed values, no SD is added.
## However, we do add ME to the pools as we do not have observed pool averages

library(tidyverse)
library(dplyr)
library(here)

#### sourcing the main source file with all the needed functions
source("Programs/method_eval_source.R")

##############################################################################################


##### Evaluation using the real Uganda data. The first set of statements evaluates the Hypred method
##### which uses individual testing for the top risk tier, and MiniPred for the middle and bottom tiers
##### the below dataset does not exist inthe repository as we did not include the real data
data <- readRDS("UgandaData/test_set_final.rds")
set.seed(1212)
## there are 3607 records in the test set, we choose 3600 to have a number divisible by 100
select <- sample.int(n=length(data$VL), size=3300, replace=FALSE)
data <- data[c(select),]


set.seed(18)
result1 <- hypred_uganda(reps=33, data=data, matsize=10, prec=10, precrd=20,
                         cutoff=1000, SE=0, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


saveRDS(result1, file=here("Results","UgandaResults","Uganda_hypred_ME0.rds"))


set.seed(18)
result1 <- hypred_uganda(reps=33, data=data, matsize=10, prec=10, precrd=20,
                         cutoff=1000, SE=0.025, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


saveRDS(result1, file=here("Results","UgandaResults","Uganda_hypred_ME.025.rds"))


set.seed(18)
result1 <- hypred_uganda(reps=33, data=data, matsize=10, prec=10, precrd=20,
                         cutoff=1000, SE=0.05, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


saveRDS(result1, file=here("Results","UgandaResults","Uganda_hypred_ME.05.rds"))


set.seed(18)
result1 <- hypred_uganda(reps=33, data=data, matsize=10, prec=10, precrd=20,
                         cutoff=1000, SE=0.075, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


saveRDS(result1, file=here("Results","UgandaResults","Uganda_hypred_ME.075.rds"))


set.seed(18)
result1 <- hypred_uganda(reps=33, data=data, matsize=10, prec=10, precrd=20,
                         cutoff=1000, SE=0.1, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


saveRDS(result1, file=here("Results","UgandaResults","Uganda_hypred_ME.1.rds"))



set.seed(18)
result1 <- hypred_uganda(reps=33, data=data, matsize=10, prec=10, precrd=20,
                         cutoff=1000, SE=.125, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


saveRDS(result1, file=here("Results","UgandaResults","Uganda_hypred_ME.125.rds"))


set.seed(18)
result1 <- hypred_uganda(reps=33, data=data, matsize=10, prec=10, precrd=20,
                         cutoff=1000, SE=.15, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


saveRDS(result1, file=here("Results","UgandaResults","Uganda_hypred_ME.15.rds"))


set.seed(18)
result1 <- hypred_uganda(reps=33, data=data, matsize=10, prec=10, precrd=20,
                         cutoff=1000, SE=.175, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


saveRDS(result1, file=here("Results","UgandaResults","Uganda_hypred_ME.175.rds"))


set.seed(18)
result1 <- hypred_uganda(reps=33, data=data, matsize=10, prec=10, precrd=20,
                         cutoff=1000, SE=.2, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


saveRDS(result1, file=here("Results","UgandaResults","Uganda_hypred_ME.2.rds"))


set.seed(18)
result1 <- hypred_uganda(reps=33, data=data, matsize=10, prec=10, precrd=20,
                         cutoff=1000, SE=.225, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


saveRDS(result1, file=here("Results","UgandaResults","Uganda_hypred_ME.225.rds"))



set.seed(18)
result1 <- hypred_uganda(reps=33, data=data, matsize=10, prec=10, precrd=20,
                         cutoff=1000, SE=.25, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


saveRDS(result1, file=here("Results","UgandaResults","Uganda_hypred_ME.25.rds"))


set.seed(18)
result1 <- hypred_uganda(reps=33, data=data, matsize=10, prec=10, precrd=20,
                         cutoff=1000, SE=.5, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


saveRDS(result1, file=here("Results","UgandaResults","Uganda_hypred_ME.5.rds"))


set.seed(18)
result1 <- hypred_uganda(reps=33, data=data, matsize=10, prec=10, precrd=20,
                         cutoff=1000, SE=.75, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


saveRDS(result1, file=here("Results","UgandaResults","Uganda_hypred_ME.75.rds"))


######################################################################################################

#### the below statements evaluate methods other than HyPred using the Uganda data
set.seed(18)
result1 <- pool.alg.cov(reps=33, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=0, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


saveRDS(result1, file=here("Results","UgandaResults","Uganda_ME0.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=33, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=0.025, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


saveRDS(result1, file=here("Results","UgandaResults","Uganda_ME.025.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=33, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=0.05, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


saveRDS(result1, file=here("Results","UgandaResults","Uganda_ME.05.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=33, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=0.075, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


saveRDS(result1, file=here("Results","UgandaResults","Uganda_ME.075.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=33, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=0.1, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


saveRDS(result1, file=here("Results","UgandaResults","Uganda_ME.1.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=33, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=0.125, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


saveRDS(result1, file=here("Results","UgandaResults","Uganda_ME.125.rds"))



set.seed(18)
result1 <- pool.alg.cov(reps=33, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=0.15, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


saveRDS(result1, file=here("Results","UgandaResults","Uganda_ME.15.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=33, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=0.175, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


saveRDS(result1, file=here("Results","UgandaResults","Uganda_ME.175.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=33, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=0.2, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


saveRDS(result1, file=here("Results","UgandaResults","Uganda_ME.2.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=33, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=0.225, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


saveRDS(result1, file=here("Results","UgandaResults","Uganda_ME.225.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=33, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.25, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


saveRDS(result1, file=here("Results","UgandaResults","Uganda_ME.25.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=33, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.5, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


saveRDS(result1, file=here("Results","UgandaResults","Uganda_ME.5.rds"))


set.seed(18)
result1 <- pool.alg.cov(reps=33, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.75, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


saveRDS(result1, file=here("Results","UgandaResults","Uganda_ME.75.rds"))
