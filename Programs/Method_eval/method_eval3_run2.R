#### evaluting the methods; sourcing the functions I need through method_eval_source

############ This is a rerun of method_eval1, using the corner matrix fill for the linreg and LRSOE methods

library(tidyverse)
library(dplyr)

# set working directory of location of source program
setwd("")
source("method_eval_source.R")
# set working directory to location to store the data
setwd("")

##############################################################################################
# METHOD_EVAL1

##### AGAIG scenarios using SD=1.0 and estimated betas
data <- read.table("Uganda_SimData_SD1.0.R")
data <- predictVL(data, b0star=0.4414945, b1star=0.1170119, b2star=1.9283680, b3star=0.1369999)



set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.05, tstperd=5, lowlimit=50, filltyp="rnd")


write.table(result1, file="Results_AGAIG_ME.05_rand.R")


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.12, tstperd=5, lowlimit=50, filltyp="rnd")


write.table(result1, file="Results_AGAIG_ME.12_rand.R")

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.25, tstperd=5, lowlimit=50, filltyp="rnd")


write.table(result1, file="Results_AGAIG_ME.25_rand.R")

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.5, tstperd=5, lowlimit=50, filltyp="rnd")


write.table(result1, file="Results_AGAIG_ME.5_rand.R")


##### No cont scenarios using SD=1.0 and set betas
data <- read.table("Uganda_SimData_SD1.0.R")
data <- predictVL(data, b0star=0.6857145, b1star=0, b2star=2.2288671, b3star= 0)


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.05, tstperd=5, lowlimit=50, filltyp="rnd")


write.table(result1, file="Results_no_cont_ME.05_rand.R")


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.12, tstperd=5, lowlimit=50, filltyp="rnd")


write.table(result1, file="Results_no_cont_ME.12_rand.R")

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.25, tstperd=5, lowlimit=50, filltyp="rnd")


write.table(result1, file="Results_no_cont_ME.25_rand.R")

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.5, tstperd=5, lowlimit=50, filltyp="rnd")


write.table(result1, file="Results_no_cont_ME.5_rand.R")


##### No binary scenarios using SD=1.0 and set betas
data <- read.table("Uganda_SimData_SD1.0.R")
data <- predictVL(data, b0star=0.6229597, b1star=0.1327599, b2star=0, b3star= 0)


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.05, tstperd=5, lowlimit=50, filltyp="rnd")


write.table(result1, file="Results_no_bin_ME.05_rand.R")


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.12, tstperd=5, lowlimit=50, filltyp="rnd")


write.table(result1, file="Results_no_bin_ME.12_rand.R")

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.25, tstperd=5, lowlimit=50, filltyp="rnd")


write.table(result1, file="Results_no_bin_ME.25_rand.R")

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.5, tstperd=5, lowlimit=50, filltyp="rnd")


write.table(result1, file="Results_no_bin_ME.5_rand.R")


##### reverse scenarios using SD=1.0 and set betas
data <- read.table("Uganda_SimData_SD1.0.R")
data <- predictVL(data, b0star=5.8155002, b1star=-0.2442230, b2star=-3.7203965, b3star=0.1260099)


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.05, tstperd=5, lowlimit=50, filltyp="rnd")


write.table(result1, file="Results_reverse_ME.05_rand.R")


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.12, tstperd=5, lowlimit=50, filltyp="rnd")


write.table(result1, file="Results_reverse_ME.12_rand.R")

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.25, tstperd=5, lowlimit=50, filltyp="rnd")


write.table(result1, file="Results_reverse_ME.25_rand.R")

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.5, tstperd=5, lowlimit=50, filltyp="rnd")


write.table(result1, file="Results_reverse_ME.5_rand.R")
