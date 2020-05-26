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
data <- read.table("Uganda_SimData_SD1.0.R")
data <- predictVL(data, b0star=0.5, b1star=0.1, b2star=2, b3star=0.1)


sum(data$VL >= 1000)/length(data$VL)


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.05, tstperd=5, lowlimit=50, filltyp="corner")


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/SD1.0_set_betas/Results_AGAIG_ME.05_corner.R")


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.12, tstperd=5, lowlimit=50, filltyp="corner")


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/SD1.0_set_betas/Results_AGAIG_ME.12_corner.R")

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.25, tstperd=5, lowlimit=50, filltyp="corner")


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/SD1.0_set_betas/Results_AGAIG_ME.25_corner.R")

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.5, tstperd=5, lowlimit=50, filltyp="corner")


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/SD1.0_set_betas/Results_AGAIG_ME.5_corner.R")


##### No cont scenarios using SD=1.0 and set betas
data <- read.table("Uganda_SimData_SD1.0.R")
data <- predictVL(data, b0star=0.5, b1star=0, b2star=3.5, b3star= 0)


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.05, tstperd=5, lowlimit=50, filltyp="corner")


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/SD1.0_set_betas/Results_no_cont_ME.05_corner.R")


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.12, tstperd=5, lowlimit=50, filltyp="corner")


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/SD1.0_set_betas/Results_no_cont_ME.12_corner.R")

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.25, tstperd=5, lowlimit=50, filltyp="corner")


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/SD1.0_set_betas/Results_no_cont_ME.25_corner.R")

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.5, tstperd=5, lowlimit=50, filltyp="corner")


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/SD1.0_set_betas/Results_no_cont_ME.5_corner.R")


##### No binary scenarios using SD=1.0 and set betas
data <- read.table("Uganda_SimData_SD1.0.R")
data <- predictVL(data, b0star=0.5, b1star=0.2, b2star=0, b3star= 0)


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.05, tstperd=5, lowlimit=50, filltyp="corner")


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/SD1.0_set_betas/Results_no_bin_ME.05_corner.R")


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.12, tstperd=5, lowlimit=50, filltyp="corner")


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/SD1.0_set_betas/Results_no_bin_ME.12_corner.R")

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.25, tstperd=5, lowlimit=50, filltyp="corner")


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/SD1.0_set_betas/Results_no_bin_ME.25_corner.R")

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.5, tstperd=5, lowlimit=50, filltyp="corner")


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/SD1.0_set_betas/Results_no_bin_ME.5_corner.R")


##### reverse scenarios using SD=1.0 and set betas
data <- read.table("Uganda_SimData_SD1.0.R")
data <- predictVL(data, b0star=6, b1star=-0.1, b2star=-2, b3star=-0.1)


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.05, tstperd=5, lowlimit=50, filltyp="corner")


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/SD1.0_set_betas/Results_reverse_ME.05_corner.R")


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.12, tstperd=5, lowlimit=50, filltyp="corner")


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/SD1.0_set_betas/Results_reverse_ME.12_corner.R")

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.25, tstperd=5, lowlimit=50, filltyp="corner")


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/SD1.0_set_betas/Results_reverse_ME.25_corner.R")



set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.5, tstperd=5, lowlimit=50, filltyp="corner")


write.table(result1, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/2nd round/SD1.0_set_betas/Results_reverse_ME.5_corner.R")
