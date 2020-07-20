### Evaluating the HyPred method in the AGAIG and reverse scenarios using the 
### simulated data with SD=1.0 and SD=0

library(tidyverse)
library(dplyr)

# set working directory of location of source program
setwd("C:/Users/Barny/Documents/GitHub/Pooled_Testing_HIV/Programs")
source("method_eval_source.R")
# set working directory to location to store the data
setwd("C:/Users/Barny/Documents/GitHub/Pooled_Testing_HIV/SimData/Records_used_in_500_matrices_hypred")

##############################################################################################


##### AGAIG scenarios using SD=1.0 and estimated betas
data <- read.table("SD1.0_data_500_run2.R")
data <- predictVL(data, b0star=0.4414945, b1star=0.1170119, b2star=1.9283680, b3star=0.1369999)


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.05, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


write.table(result1, file="C:/Users/Barny/Documents/GitHub/Pooled_Testing_HIV/Results/SimResults/round2/Results_Hypred_AGAIG_SD1_ME.05_rand.R")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.12, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


write.table(result1, file="C:/Users/Barny/Documents/GitHub/Pooled_Testing_HIV/Results/SimResults/round2/Results_Hypred_AGAIG_SD1_ME.12_rand.R")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.25, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


write.table(result1, file="C:/Users/Barny/Documents/GitHub/Pooled_Testing_HIV/Results/SimResults/round2/Results_Hypred_AGAIG_SD1_ME.25_rand.R")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.5, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


write.table(result1, file="C:/Users/Barny/Documents/GitHub/Pooled_Testing_HIV/Results/SimResults/round2/Results_Hypred_AGAIG_SD1_ME.5_rand.R")


##################################################################################################################################

##### AGAIG scenarios using SD=0 and estimated betas
data <- read.table("SD0_data_500_run2.R")
data <- predictVL(data, b0star=0.4414945, b1star=0.1170119, b2star=1.9283680, b3star=0.1369999)


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=0, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


write.table(result1, file="C:/Users/Barny/Documents/GitHub/Pooled_Testing_HIV/Results/SimResults/round2/Results_Hypred_AGAIG_SD0_ME0_rand.R")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.25, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


write.table(result1, file="C:/Users/Barny/Documents/GitHub/Pooled_Testing_HIV/Results/SimResults/round2/Results_Hypred_AGAIG_SD0_ME.25_rand.R")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.5, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


write.table(result1, file="C:/Users/Barny/Documents/GitHub/Pooled_Testing_HIV/Results/SimResults/round2/Results_Hypred_AGAIG_SD0_ME.5_rand.R")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.75, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


write.table(result1, file="C:/Users/Barny/Documents/GitHub/Pooled_Testing_HIV/Results/SimResults/round2/Results_Hypred_AGAIG_SD0_ME.75_rand.R")

######################################################################################################################################

##### Reverse scenarios using SD=1.0 and estimated betas
data <- read.table("SD1.0_data_500_run2.R")
data <- predictVL(data, b0star=5.8155002, b1star=-0.2442230, b2star=-3.7203965, b3star=0.1260099)


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.05, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


write.table(result1, file="C:/Users/Barny/Documents/GitHub/Pooled_Testing_HIV/Results/SimResults/round2/Results_Hypred_reverse_SD1_ME.05_rand.R")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.12, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


write.table(result1, file="C:/Users/Barny/Documents/GitHub/Pooled_Testing_HIV/Results/SimResults/round2/Results_Hypred_reverse_SD1_ME.12_rand.R")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.25, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


write.table(result1, file="C:/Users/Barny/Documents/GitHub/Pooled_Testing_HIV/Results/SimResults/round2/Results_Hypred_reverse_SD1_ME.25_rand.R")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.5, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


write.table(result1, file="C:/Users/Barny/Documents/GitHub/Pooled_Testing_HIV/Results/SimResults/round2/Results_Hypred_reverse_SD1_ME.5_rand.R")


##################################################################################################################################

##### Reverse scenarios using SD=0 and estimated betas
data <- read.table("SD0_data_500_run2.R")
data <- predictVL(data, b0star=5.8155002, b1star=-0.2442230, b2star=-3.7203965, b3star=0.1260099)



set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=0, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


write.table(result1, file="C:/Users/Barny/Documents/GitHub/Pooled_Testing_HIV/Results/SimResults/round2/Results_Hypred_reverse_SD0_ME0_rand.R")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.25, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


write.table(result1, file="C:/Users/Barny/Documents/GitHub/Pooled_Testing_HIV/Results/SimResults/round2/Results_Hypred_reverse_SD0_ME.25_rand.R")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.5, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


write.table(result1, file="C:/Users/Barny/Documents/GitHub/Pooled_Testing_HIV/Results/SimResults/round2/Results_Hypred_reverse_SD0_ME.5_rand.R")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.75, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


write.table(result1, file="C:/Users/Barny/Documents/GitHub/Pooled_Testing_HIV/Results/SimResults/round2/Results_Hypred_reverse_SD0_ME.75_rand.R")

######################################################################################################################################

