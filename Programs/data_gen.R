#==============================================================================
# FILENAME: data_gen.R
# PROJECT: 	Pooled testing in HIV
# PURPOSE: generate simulated dataset for pooled testing method evaluation 
#          the data are simulated to follow the VL distribution in Uganda
# AUTHOR: Adam Brand


# OUTPUT: 6 datasets
#         - dataset for method evaluation with SD=1.0
#         - dataset for method evaluation with SD=0
#         - simulated training set with SD=1.0
#         - simulated training set with SD=0
#         - reverse training dataset with SD=1.0
#         - reverse training dataset with SD=0

# R VERSION: 3.6.1
#==============================================================================
#Notes: 





# =============================================================================

##### Data generation program AB 2020-03-03
##### Testing
library(shiny)
library(dplyr)
library(here)

# the below source program will only work with access to the real Uganda data
# this will produce errors when running without the uganda data, but
# the simulation data will still be generated. It was used to compare the VL distributions
# of the real data to the simulated data to ensure similarity
source(here("Programs","raw_data_clean.R"))

# This function simulates viral load (VL) values under 500 copies/mL
# n is the number of values to simulate
Under500VL <- function(n){
  a <- round(.85*n)
  b <- round(.05*n)
  c <- n - a - b
  # 85% ~ U(0,50); 5% ~ U(50,100); 10% ~ U(100,500)
  y <- c(runif(a,0,50), runif(b,50,100),runif(c,100,500))
  # returns a vector of viral load values under 500 copies/mL
  return(y)}

# Simulates VL values above 500 copies/mL
# n is the number of values to simulate
Over500VL <- function(n){
  a <- round(.93*n)
  b <- n - a
  # 93% ~ Gamma(1.6, 0.5); 7% ~ Gamma(3.18, 0.5); both shifted 2.7 right
  y <- c(replicate(a, 10^(2.7 + rgamma(1,shape = 1.6, scale = 0.5))),
         replicate(b, 10^(2.7 + rgamma(1,shape = 3.18, scale = 0.5))))
  # returns a vector of viral load values over 500 copies/mL
  return(y)}

# Converts a failure prevalence at any cutoff into a failure prevalence above 500
prevover500 <- function(prevfail, cutoff){
  # This uses the cumulative distribution functions of viral load values to 
  # convert prevalence at any cutoff into prevalence over 500.
  # Note:  this depends on where the original failure cutoff is.
  if (50 <= cutoff & cutoff < 100){
    prevover500 = (prevfail - .10 - (.05*(100-cutoff)/50))/(1 - .1 - (.05*(100-cutoff)/50))}
  else if (100 <= cutoff & cutoff < 500){
    prevover500 = (prevfail - (.10*(500-cutoff)/400))/(1 - (.10*(500-cutoff)/400))}
  else if (cutoff == 500){
    prevover500 = prevfail}
  else if (cutoff > 500){
    prevover500 = prevfail/((.93*(1 - pgamma(log10(cutoff)-2.7, shape = 1.6, scale = 0.5)) +
                               .07*(1 - pgamma(log10(cutoff)-2.7, shape = 3.18, scale = 0.5))))}
  else {prevover500 = "error"}
  # returns an euivalent failure prevalence for a cutoff of 500
  return(prevover500)
}


# Creates a vector of realistically skewed viral load values
# n is the number of VL's to simulate
# prevover500 is the prevalence above the cutoff (assumed 500)
# cutoff is the cutoff defining treatment failure
# prevfail is a boolean.  If true, the function will interpret
# prevover500 as the prevalence over the cutoff
createpop <- function(n, prevover500, cutoff, prevfail=FALSE){
  if (prevfail==FALSE){
    prevover500 = prevover500}
  # If prevfail=TRUE, this converts the prevover500
  # into an actual prevalence over 500 copies/mL
  else {prevover500 <- prevover500(prevover500, cutoff)}
  # Calculates the number over 500
  a <- round(n*prevover500)
  # Calculates the number under 500
  b <- n - a
  # Simulates VL's based on mixed distribution by May et al. (2020)
  over500 <- Over500VL(a)
  under500 <- Under500VL(b)
  pop <- c(over500, under500)
  # reutrns vector of VL's
  return(pop)
}

#Next 2 lines checking the create pop function at our Uganda failure prevalence of 0.05664436
set.seed(121212)
check <- createpop(100000, prevover500=0.05664436, cutoff=1000, prevfail=TRUE)
sum(check >= 1000)/length(check)
for (i in 1:length(check)){
  if (check[i] <=50){check[i]=50}
}

#### checking our simulated values against the Uganda values
#### Uganda VLs come from raw_data_clean, VLdata$VL10

par(mfrow=c(1,2))
##Uganda data
hist(VLdata$VL10,
     main="Distribution of Uganda VLs",
     xlab="Log 10 Viral Load",
     xlim=c(0,10),
     ylim=c(0,4),
     breaks=25,
     prob=TRUE)


lines(density(VLdata$VL10))

## Simulated data
hist(log10(check),
     main="Distribution of simulated VLs",
     xlab="Log 10 Viral Load",
     xlim=c(0,10),
     ylim=c(0,4),
     breaks=25,
     prob=TRUE)

lines(density(log10(check)))

## Code to run the shinypopfit program.
## this shiny program, shinypopfit, was written to visualize similarities between
## a skewed, population viral load distribution created by our function createpop
## and a distribution of VLs generated by a model with an intercept (B0), a continuous
## variable (B1), a binary variable (B2) and the interaction (B3) on the log 10 scale
## Using this app along with diagnostics below, we attempted to simulate VLs that mimicked
## the VLs in Uganda using a supposed prediction model with simulated predictors.

## to run the app, set the working directory to the location of the shinypopfit folder
runApp(here("Programs","shinypopfit"))

# Simulates a dataframe with adherence, PFS and viral load for n patients
# based on the covariate model, described in detail later
# Adherence data ~ Beta(shape1, shape2)
# b0 - b4 are the coefficients of our covariate model
# pf.prev is the prevalence of previous failure
# sd is the standard deviation for person-to-person variation
VL.data.adh <- function(n, shape1, shape2, b0, b1, b2, b3, pf.prev, sd){
  # simulates n previous failure statuses using independent Bernoullis
  # where p is the input for pf.prev
  pf <- rbinom(n=n, size=1, prob=pf.prev)
  # Simulates adherence data for n patients based on distr. pars.
  adhere <- (100 - 100*rbeta(n=n, shape1 = shape1, shape2=shape2))
  # simulated person to person variablilty
  var <- rnorm(n=n, mean=0, sd = sd)
  # generates log10 VL from our covariate model
  log.VL <- b0 + b1*adhere + b2*pf + b3*adhere*pf + var
  # getting the model VLs without variability
  model.VL <- b0 + b1*adhere + b2*pf + b3*adhere*pf
  # converts to actual viral load values
  VL <- 10^log.VL
  data <- cbind(VL, log.VL, model.VL, adhere, pf)
  # reutrns a data frame with  VL, adherence and PFS data
  # for each patient
  return(data.frame(data))
}



# Generates 1,000,000 patient data at 5.664436% failure prevalence over 1000 copies/mL
# 
set.seed(1212)
simdata <- VL.data.adh(n=1000000, shape1=70, shape2=1.5, b0=0.5, b1=.1, b2=2, b3=.1, pf.prev=.1, sd=1.0)

# truncating values at around 3,000,000 to better mimic actual values
# also include a little error, so they aren't exactly the same
for (i in 1:length(simdata$VL)){
  if (simdata$VL[i] > 3000000) {simdata$VL[i] = 3000000 + rnorm(n=1, mean=0, sd = 10000)}
}


# creating a categorical classification of VL
for (i in 1:length(simdata$VL)){
  if (simdata$VL[i] <= 50) {simdata$VLcat[i]="low"}
  if (simdata$VL[i] > 50 & simdata$VL[i]<1000) {simdata$VLcat[i]="medium"}
  if (simdata$VL[i] >= 1000 & simdata$VL[i] < 10000) {simdata$VLcat[i]="high"}
  if (simdata$VL[i] >= 10000) {simdata$VLcat[i]="vhigh"}
}


#### getting a smaller dataset as training sim data to fit a model to
#### we will then evaluate the methods using the estimated betas fit to this model
set.seed(2121)
simdata_train <- VL.data.adh(n=5000, shape1=70, shape2=1.5, b0=0.5, b1=.1, b2=2, b3=.1, pf.prev=.1, sd=1.0)


for (i in 1:length(simdata_train$VL)){
  if (simdata_train$VL[i] > 3000000) {simdata_train$VL[i] = 3000000 + rnorm(n=1, mean=0, sd = 10000)}
}

for (i in 1:length(simdata_train$VL)){
  if (simdata_train$VL[i] <= 50) {simdata_train$VLcat[i]="low"}
  if (simdata_train$VL[i] > 50 & simdata_train$VL[i]<1000) {simdata_train$VLcat[i]="medium"}
  if (simdata_train$VL[i] >= 1000 & simdata_train$VL[i] < 10000) {simdata_train$VLcat[i]="high"}
  if (simdata_train$VL[i] >= 10000) {simdata_train$VLcat[i]="vhigh"}
}

## getting the max adhere for the reverse scenario, setting to just over 0.1, so 
## everyone has a positive adhere when we subtract their actual adhere
max_adhere <- max(simdata_train$adhere) + 0.1

#### creating a reverse dataset for the reverse scenario
simdata_rev <- simdata_train

### reversing the direction of the covariates
for (i in 1:length(simdata_rev$VL)){
  simdata_rev$pf[i] <- 1 - simdata_train$pf[i]
  simdata_rev$adhere[i] <- max_adhere - simdata_train$adhere[i]
}

### checking simdata against Uganda data
par(mfrow=c(1,2))
##Uganda data
hist(VLdata$VL10,
     main="Distribution of Uganda VLs",
     xlab="Log 10 Viral Load",
     xlim=c(0,10),
     ylim=c(0,4),
     breaks=25,
     freq=FALSE)


lines(density(VLdata$VL10))

## Simulated data
hist(log10(simdata$VL),
     main="Distribution of simulated VLs",
     xlab="Log 10 Viral Load",
     xlim=c(0,10),
     ylim=c(0,4),
     breaks=25,
     freq=FALSE)

lines(density(simdata$VL))

sum(simdata$VL >=1000)/length(simdata$VL)

summary(simdata$VL)
summary(VLdata$VL)


saveRDS(simdata, file=here("SimData","Uganda_SimData_SD1.0.rds"))
saveRDS(simdata_train, file=here("SimData","Uganda_SimData_train_SD1.0.rds"))
saveRDS(simdata_rev, file=here("SimData","Uganda_SimData_train_SD1.0_reverse.rds"))


### checking that the model isnt too overly predictive. We set the sd at 1.0 on log10 scale to get 0.327 R2
check <- readRDS("SimData/Uganda_SimData_SD1.0.rds")
plot(check$log.VL ~ check$model.VL)

### R squared between model and actual data = 0.32 on log10 scale, 0.47 on the raw scale. SD=1.0
cor(check$log.VL, check$model.VL)^2
check$model10 <- 10^check$model.VL
cor(check$VL, check$model10)^2

##prevalence of treatment failure is 5.7%
sum(check$VL>=1000)/length(check$VL)

### checking the difference in treatment failures between true and observed using measurement error
add_me <- function(data, ME){
  n <- length(data[,1])
  me <- rnorm(n, mean=0, sd=ME)
  obs <- 10^(me+data$log.VL)
  return(cbind(data, me, obs))
}

check.05 <- add_me(check, ME=.05)
check.1 <- add_me(check, ME=.1)
check.15 <- add_me(check, ME=.15)
check.2 <- add_me(check, ME=.2)
check.25 <- add_me(check, ME=.25)
check.5 <- add_me(check, ME=.5)
check.75 <- add_me(check, ME=.75)

assay_sens <- function(data){
  denom <- sum(data$VL>=1000)
  num <- sum(data$VL>=1000 & data$obs>=1000)
  sens <- num/denom
  return(sens)
}

assay_spec <- function(data){
  denom <- sum(data$VL<1000)
  num <- sum(data$VL<1000 & data$obs<1000)
  spec <- num/denom
  return(spec)
}
## 0.9778
assay_sens(check.05)
# 0.9986
assay_spec(check.05)

## 0.9580
assay_sens(check.1)
# 0.9970
assay_spec(check.1)

## 0.9370
assay_sens(check.15)
# 0.9954
assay_spec(check.15)

## 0.9205
assay_sens(check.2)
# 0.9937
assay_spec(check.2)

## 0.9024
assay_sens(check.25)
# 0.9918
assay_spec(check.25)

## 0.8274
assay_sens(check.5)
# 0.9801
assay_spec(check.5)

## 0.7717
assay_sens(check.75)
# 0.9634
assay_spec(check.75)


##### getting and outputting data with no vairation; variation will all be due to ME
set.seed(1212)
simdata2 <- VL.data.adh(n=1000000, shape1=56, shape2=2, b0=0.5, b1=.1, b2=2, b3=.1, pf.prev=.1, sd=0)

for (i in 1:length(simdata2$VL)){
  if (simdata2$VL[i] > 3000000) {simdata2$VL[i] = 3000000 + rnorm(n=1, mean=0, sd = 10000)}
}

for (i in 1:length(simdata2$VL)){
  if (simdata2$VL[i] <= 50) {simdata2$VLcat[i]="low"}
  if (simdata2$VL[i] > 50 & simdata2$VL[i]<1000) {simdata2$VLcat[i]="medium"}
  if (simdata2$VL[i] >= 1000 & simdata2$VL[i] < 10000) {simdata2$VLcat[i]="high"}
  if (simdata2$VL[i] >= 10000) {simdata2$VLcat[i]="vhigh"}
}


set.seed(2121)
simdata2_train <- VL.data.adh(n=5000, shape1=56, shape2=2, b0=0.5, b1=.1, b2=2, b3=.1, pf.prev=.1, sd=0)


for (i in 1:length(simdata2_train$VL)){
  if (simdata2_train$VL[i] > 3000000) {simdata2_train$VL[i] = 3000000 + rnorm(n=1, mean=0, sd = 10000)}
}

#### this creates categories for VL for the sbrl model
for (i in 1:length(simdata2_train$VL)){
  if (simdata2_train$VL[i] <= 50) {simdata2_train$VLcat[i]="low"}
  if (simdata2_train$VL[i] > 50 & simdata2_train$VL[i]<1000) {simdata2_train$VLcat[i]="medium"}
  if (simdata2_train$VL[i] >= 1000 & simdata2_train$VL[i] < 10000) {simdata2_train$VLcat[i]="high"}
  if (simdata2_train$VL[i] >= 10000) {simdata2_train$VLcat[i]="vhigh"}
}


## getting the max adhere for the reverse scenario, setting to just over 0.1, so 
## everyone has a positive adhere when we subtract their actual adhere
max_adhere <- max(simdata2_train$adhere) + 0.1

#### creating a reverse dataset for the reverse scenario
simdata2_rev <- simdata2_train

for (i in 1:length(simdata2_rev$VL)){
  simdata2_rev$pf[i] <- 1 - simdata2_train$pf[i]
  simdata2_rev$adhere[i] <- max_adhere - simdata2_train$adhere[i]
}

sum(simdata2$VL >=1000)/length(simdata2$VL)
sum(simdata2_train$VL >=1000)/length(simdata2_train$VL)


saveRDS(simdata2, file=here("SimData","Uganda_SimData_SD0.rds"))
saveRDS(simdata2_train, file=here("SimData","Uganda_SimData_train_SD0.rds"))
saveRDS(simdata2_rev,  file=here("SimData","Uganda_SimData_train_SD0_reverse.rds"))


check <- readRDS("SimData/Uganda_SimData_SD0.rds")

### R squared between model and actual data
cor(check$log.VL, check$model.VL)^2


##prevalence of treatment failure is 5.8%
sum(check$VL>=1000)/length(check$VL)


check.05 <- add_me(check, ME=.05)
check.1 <- add_me(check, ME=.1)
check.15 <- add_me(check, ME=.15)
check.2 <- add_me(check, ME=.2)
check.25 <- add_me(check, ME=.25)
check.5 <- add_me(check, ME=.5)
check.75 <- add_me(check, ME=.75)


## 0.9674
assay_sens(check.05)
# 0.9979
assay_spec(check.05)

## 0.9367
assay_sens(check.1)
# 0.9958
assay_spec(check.1)

## 0.9069
assay_sens(check.15)
# 0.9934
assay_spec(check.15)

## 0.8804
assay_sens(check.2)
# 0.9916
assay_spec(check.2)

## 0.8542
assay_sens(check.25)
# 0.9900
assay_spec(check.25)

## 0.7655
assay_sens(check.5)
# 0.9845
assay_spec(check.5)

## 0.7046
assay_sens(check.75)
# 0.9793
assay_spec(check.75)
