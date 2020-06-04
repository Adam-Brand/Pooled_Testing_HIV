## Evaluating the pooled testing methods using the real Uganda data
## evaluating the Hypred method with the first statements and the other methods below.
## As the real data are the observed values, no SD is added.
## However, we do add ME to the pools as we do not have observed pool averages

library(tidyverse)
library(dplyr)

# set working directory of location of source program
setwd("")
source("method_eval_source.R")
# set working directory to location to store the data
setwd("")

##############################################################################################


##### Evaluation using the real Uganda data. The first set of statements evaluates the Hypred method
##### which uses individual testing for the top risk tier, and MiniPred for the middle and bottom tiers
##### the below dataset does not exist inthe repository as we did not include the real data
data <- read.table("test_set_final.R")
set.seed(1212)
## there are 3607 records in the test set, we choose 3600 to have a number divisible by 100
select <- sample.int(n=length(data$VL), size=3600, replace=FALSE)
data <- data[c(select),]


set.seed(18)
result1 <- hypred_uganda(reps=36, data=data, matsize=10, prec=10, precrd=20,
                         cutoff=1000, SE=0, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


write.table(result1, file="Uganda_hypred_ME0.R")


set.seed(18)
result1 <- hypred_uganda(reps=36, data=data, matsize=10, prec=10, precrd=20,
                         cutoff=1000, SE=0.05, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


write.table(result1, file="Uganda_hypred_ME.05.R")



set.seed(18)
result1 <- hypred_uganda(reps=36, data=data, matsize=10, prec=10, precrd=20,
                         cutoff=1000, SE=.12, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


write.table(result1, file="Uganda_hypred_ME.12.R")



set.seed(18)
result1 <- hypred_uganda(reps=36, data=data, matsize=10, prec=10, precrd=20,
                         cutoff=1000, SE=.25, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


write.table(result1, file="Uganda_hypred_ME.25.R")


set.seed(18)
result1 <- hypred_uganda(reps=36, data=data, matsize=10, prec=10, precrd=20,
                         cutoff=1000, SE=.5, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


write.table(result1, file="Uganda_hypred_ME.5.R")


set.seed(18)
result1 <- hypred_uganda(reps=36, data=data, matsize=10, prec=10, precrd=20,
                         cutoff=1000, SE=.75, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


write.table(result1, file="Uganda_hypred_ME.75.R")


######################################################################################################

#### the below statements evaluate methods other than HyPred using the Uganda data
set.seed(18)
result1 <- pool.alg.cov(reps=36, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=0, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


write.table(result1, file="Uganda_ME0.R")



set.seed(18)
result1 <- pool.alg.cov(reps=36, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.05, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


write.table(result1, file="Uganda_ME.05.R")



set.seed(18)
result1 <- pool.alg.cov(reps=36, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.12, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


write.table(result1, file="Uganda_ME.12.R")



set.seed(18)
result1 <- pool.alg.cov(reps=36, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.25, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


write.table(result1, file="Uganda_ME.25.R")


set.seed(18)
result1 <- pool.alg.cov(reps=36, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.5, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


write.table(result1, file="Uganda_ME.5.R")


set.seed(18)
result1 <- pool.alg.cov(reps=36, data=data, matsize=10, prec=10, precrd=20,
                        cutoff=1000, SE=.75, tstperd=5, lowlimit=50, filltyp="rnd", Uganda=TRUE)


write.table(result1, file="Uganda_ME.75.R")