### Evaluating the methods (except for Hypred) in the AGAIG and reverse scenarios using the 
### simulated data with SD=0

library(tidyverse)
library(dplyr)

# set working directory of location of source program
setwd("")
source("method_eval_source.R")
# set working directory to location to store the data
setwd("")

##############################################################################################


##### AGAIG scenarios using SD=0 and estimated betas
data <- read.table("Uganda_SimData_SD0.R")
#### these betas were estimated using the sim_data_betas program
#### the betas are estimated using the training set with SD=1.0. Using the training set with SD=0
#### provided 'too accurate' results due to the lack of variation
data <- predictVL(data, b0star=0.4414945, b1star=0.1170119, b2star=1.9283680, b3star=0.1369999)



set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=0, tstperd=5, lowlimit=50, filltyp="rnd")


write.table(result1, file="Results_AGAIG_SD0_ME0_rand.R")


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.25, tstperd=5, lowlimit=50, filltyp="rnd")


write.table(result1, file="Results_AGAIG_SD0_ME.25_rand.R")

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.5, tstperd=5, lowlimit=50, filltyp="rnd")


write.table(result1, file="Results_AGAIG_SD0_ME.5_rand.R")

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.75, tstperd=5, lowlimit=50, filltyp="rnd")


write.table(result1, file="Results_AGAIG_SD0_ME.75_rand.R")



##### reverse scenarios using SD=0 and estimated betas from the reverse training set
data <- read.table("Uganda_SimData_SD0.R")
data <- predictVL(data, b0star=5.8155002, b1star=-0.2442230, b2star=-3.7203965, b3star=0.1260099)


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=0, tstperd=5, lowlimit=50, filltyp="rnd")


write.table(result1, file="Results_reverse_SD0_ME0_rand.R")


set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.25, tstperd=5, lowlimit=50, filltyp="rnd")


write.table(result1, file="Results_reverse_SD0_ME.25_rand.R")

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.5, tstperd=5, lowlimit=50, filltyp="rnd")


write.table(result1, file="Results_reverse_SD0_ME.5_rand.R")

set.seed(18)
result1 <- pool.alg.cov(reps=500, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.75, tstperd=5, lowlimit=50, filltyp="rnd")


write.table(result1, file="Results_reverse_SD0_ME.75_rand.R")
