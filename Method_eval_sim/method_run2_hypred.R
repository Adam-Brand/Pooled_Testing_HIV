#### evaluting the methods; sourcing the functions I need through method_eval_source

############ This is a rerun of method_eval1, using the corner matrix fill for the linreg and LRSOE methods

library(tidyverse)
library(dplyr)

setwd("C:/Users/Barny/Dropbox/KI_Project_4/Programs")
source("method_eval_source_new_hypred.R")
setwd("C:/Users/Barny/Dropbox/KI_Project_4/Data/SimData/Round_2_data")

##############################################################################################
# METHOD_EVAL1

##### AGAIG scenarios using SD=1.0 and set betas
data <- read.table("SD1.0_data_500_run2.R")
data <- predictVL(data, b0star=0.4414945, b1star=0.1170119, b2star=1.9283680, b3star=0.1369999)


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.05, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/HyPred_SD1.0/Results_Hypred_AGAIG_ME.05_rand.R")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.12, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/HyPred_SD1.0/Results_Hypred_AGAIG_ME.12_rand.R")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.25, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/HyPred_SD1.0/Results_Hypred_AGAIG_ME.25_rand.R")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.5, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/HyPred_SD1.0/Results_Hypred_AGAIG_ME.5_rand.R")


##################################################################################################################################

data <- read.table("SD0_data_500_run2.R")
data <- predictVL(data, b0star=0.4414945, b1star=0.1170119, b2star=1.9283680, b3star=0.1369999)


debug(hypred)
set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=0, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/HyPred_SD0/Results_Hypred_AGAIG_ME0_rand.R")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.25, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/HyPred_SD0/Results_Hypred_AGAIG_ME.25_rand.R")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.5, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/HyPred_SD0/Results_Hypred_AGAIG_ME.5_rand.R")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.75, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/HyPred_SD0/Results_Hypred_AGAIG_ME.75_rand.R")

######################################################################################################################################

##### AGAIG scenarios using SD=1.0 and set betas
data <- read.table("SD1.0_data_500_run2.R")
data <- predictVL(data, b0star=5.8155002, b1star=-0.2442230, b2star=-3.7203965, b3star=0.1260099)


debug(hypred)
set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.05, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/HyPred_SD1.0/Results_Hypred_reverse_ME.05_rand.R")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.12, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/HyPred_SD1.0/Results_Hypred_reverse_ME.12_rand.R")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.25, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/HyPred_SD1.0/Results_Hypred_reverse_ME.25_rand.R")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.5, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/HyPred_SD1.0/Results_Hypred_reverse_ME.5_rand.R")


##################################################################################################################################

data <- read.table("SD0_data_500_run2.R")
data <- predictVL(data, b0star=5.8155002, b1star=-0.2442230, b2star=-3.7203965, b3star=0.1260099)



set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=0, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/HyPred_SD0/Results_Hypred_reverse_ME0_rand.R")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.25, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/HyPred_SD0/Results_Hypred_reverse_ME.25_rand.R")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.5, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/HyPred_SD0/Results_Hypred_reverse_ME.5_rand.R")


set.seed(18)
result1 <- hypred(reps=500, data=data, matsize=10, prec=10, precrd=20,
                  cutoff=1000, SE=.75, tstperd=5, lowlimit=50, filltyp="rnd", top_percent=0.1, bot_percent=0.6)


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/HyPred_SD0/Results_Hypred_reverse_ME.75_rand.R")

######################################################################################################################################

