#### evaluting the methods; sourcing the functions I need through method_eval_source

############ This is a rerun of method_eval1, using the corner matrix fill for the linreg and LRSOE methods

library(tidyverse)
library(dplyr)

setwd("C:/Users/Barny/Dropbox/KI_Project_4/Programs")
source("method_eval_source_new.R")
setwd("C:/Users/Barny/Dropbox/KI_Project_4/Data/SimData")

##############################################################################################
# METHOD_EVAL1

##### AGAIG scenarios using SD=1.0 and set betas
data <- read.table("Uganda_SimData_SD0.R")
data <- predictVL(data, b0star=0.4414945, b1star=0.1170119, b2star=1.9283680, b3star=0.1369999)



set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=0, tstperd=5, lowlimit=50, filltyp="corner")


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/SD0_est_betas/Results_AGAIG_ME0_corner.R")


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.25, tstperd=5, lowlimit=50, filltyp="corner")


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/SD0_est_betas/Results_AGAIG_ME.25_corner.R")

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.5, tstperd=5, lowlimit=50, filltyp="corner")


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/SD0_est_betas/Results_AGAIG_ME.5_corner.R")

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.75, tstperd=5, lowlimit=50, filltyp="corner")


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/SD0_est_betas/Results_AGAIG_ME.75_corner.R")


##### No cont scenarios using SD=1.0 and set betas
data <- read.table("Uganda_SimData_SD0.R")
data <- predictVL(data, b0star=0.6857145, b1star=0, b2star=2.2288671, b3star= 0)


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=0, tstperd=5, lowlimit=50, filltyp="corner")


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/SD0_est_betas/Results_no_cont_ME0_corner.R")


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.25, tstperd=5, lowlimit=50, filltyp="corner")


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/SD0_est_betas/Results_no_cont_ME.25_corner.R")

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.5, tstperd=5, lowlimit=50, filltyp="corner")


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/SD0_est_betas/Results_no_cont_ME.5_corner.R")

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.75, tstperd=5, lowlimit=50, filltyp="corner")


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/SD0_est_betas/Results_no_cont_ME.75_corner.R")


##### No binary scenarios using SD=1.0 and set betas
data <- read.table("Uganda_SimData_SD0.R")
data <- predictVL(data, b0star=0.6229597, b1star=0.1327599, b2star=0, b3star= 0)


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=0, tstperd=5, lowlimit=50, filltyp="corner")


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/SD0_est_betas/Results_no_bin_ME0_corner.R")


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.25, tstperd=5, lowlimit=50, filltyp="corner")


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/SD0_est_betas/Results_no_bin_ME.25_corner.R")

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.5, tstperd=5, lowlimit=50, filltyp="corner")


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/SD0_est_betas/Results_no_bin_ME.5_corner.R")

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.75, tstperd=5, lowlimit=50, filltyp="corner")


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/SD0_est_betas/Results_no_bin_ME.75_corner.R")


##### reverse scenarios using SD=1.0 and set betas
data <- read.table("Uganda_SimData_SD0.R")
data <- predictVL(data, b0star=5.8155002, b1star=-0.2442230, b2star=-3.7203965, b3star=0.1260099)


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=0, tstperd=5, lowlimit=50, filltyp="corner")


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/SD0_est_betas/Results_reverse_ME0_corner.R")


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.25, tstperd=5, lowlimit=50, filltyp="corner")


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/SD0_est_betas/Results_reverse_ME.25_corner.R")

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.5, tstperd=5, lowlimit=50, filltyp="corner")


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/SD0_est_betas/Results_reverse_ME.5_corner.R")

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.75, tstperd=5, lowlimit=50, filltyp="corner")


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/SD0_est_betas/Results_reverse_ME.75_corner.R")
