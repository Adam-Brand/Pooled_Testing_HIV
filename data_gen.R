##### Data generation program AB 2020-03-03
library(shiny)
library(dplyr)
setwd("C:/Users/Barny/Dropbox/KI_Project_4/Programs")
source("raw_data_clean.R")

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

# par(mfrow=c(1,2))
# ##Uganda data
# hist(VLdata$VL10,
#      main="Distribution of Uganda VLs",
#      xlab="Log 10 Viral Load",
#      xlim=c(0,10),
#      ylim=c(0,4),
#      breaks=25,
#      prob=TRUE)
# 
# 
# lines(density(VLdata$VL10))
# 
# ## Simulated data
# hist(log10(check),
#      main="Distribution of simulated VLs",
#      xlab="Log 10 Viral Load",
#      xlim=c(0,10),
#      ylim=c(0,4),
#      breaks=25,
#      prob=TRUE)
# 
# lines(density(log10(check)))

## code to run the shinypopfit program

# setwd("C:/Users/Barny/Dropbox/KI_Project_4/Programs")
# runApp("shinypopfit")

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

#Here is the code that generated the simulation data for the thesis

# Generates 1 million patient data at 1% failure prevalence over 1000 copies/mL
#set.seed(8)
#data125 <- VL.data.adh(n=1000000, shape1=8, shape2 = 0.05, b0=1, b1=.05, b2=1, b3=.05, pf.prev=.25, sd=0.5)
#write.table(data125, file="data125.shapes.8.0.0.05betas1..05.1..05sd.5")

# Generates 1 million patient data at 5% failure prevalence over 1000 copies/mL
#set.seed(11)
#data525 <- VL.data.adh(n=1000000, shape1=7, shape2 = 0.3, b0=1, b1=.05, b2=1, b3=.05, pf.prev=.25, sd=0.5)
#write.table(data525, file="data525.shapes.7.0.0.3betas1..05.1..05sd.5")

# Generates 1 million patient data at 10% failure prevalence over 1000 copies/mL
#set.seed(12)
#data1025 <- VL.data.adh(n=1000000, shape1=4, shape2 = 0.33, b0=1, b1=.05, b2=1, b3=.05, pf.prev=.25, sd=0.5)
#write.table(data1025, file="data1025.shapes.4.0.0.33betas1..05.1..05sd.5")

# Generates 1 million patient data at 20% failure prevalence over 1000 copies/mL
#set.seed(15)
#data2025 <- VL.data.adh(n=1000000, shape1=6, shape2 = 1.05, b0=1, b1=.05, b2=1, b3=.05, pf.prev=.25, sd=0.5)
#write.table(data2025, file="data2025.shapes.6.0.1.05betas1..05.1..05sd.5")


# Generates 50,000 patient data at 5.664436% failure prevalence over 1000 copies/mL
# 
set.seed(1212)
simdata <- VL.data.adh(n=1000000, shape1=70, shape2=1.5, b0=0.5, b1=.1, b2=2, b3=.1, pf.prev=.1, sd=1.0)
## setting lower limit of detection (50)
# for (i in 1:length(simdata$VL)){
#   if (simdata$VL[i] <= 50) {simdata$VL[i]=50}
# }

for (i in 1:length(simdata$VL)){
  if (simdata$VL[i] > 3000000) {simdata$VL[i] = 3000000 + rnorm(n=1, mean=0, sd = 10000)}
}


check <- simdata %>%
  arrange(
    desc(model.VL)
  )

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
# setting lower limit of detection (50)
# for (i in 1:length(simdata_train$VL)){
#   if (simdata_train$VL[i] <= 50) {simdata_train$VL[i]=50}
# }

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

#setwd("C:/Users/adabra/Dropbox/KI_Project_4/Data/SimData")
setwd("C:/Users/Barny/Dropbox/KI_Project_4/Data/SimData")
 write.table(simdata, file="Uganda_SimData_SD1.0.R")
 write.table(simdata_train, file="Uganda_SimData_train_SD1.0.R")
 write.table(simdata_rev, file="Uganda_SimData_train_SD1.0_reverse.R")


### checking that the model isnt too overly predictive. We set the sd at 1.0 on log10 scale to get 0.327 R2
check <- read.table("Uganda_SimData_SD1.0.R")
plot(check$log.VL ~ check$model.VL)

### R squared between model and actual data = 0.327. SD=1.0
cor(check$log.VL, check$model.VL)^2



##### getting and outputting data with no vairation; variation will all be due to ME
set.seed(1212)
simdata2 <- VL.data.adh(n=1000000, shape1=56, shape2=2, b0=0.5, b1=.1, b2=2, b3=.1, pf.prev=.1, sd=0)
# setting lower limit of detection (50)
# for (i in 1:length(simdata2$VL)){
#   if (simdata2$VL[i] <= 50) {simdata2$VL[i]=50}
# }

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
## setting lower limit of detection (50)
# for (i in 1:length(simdata2_train$VL)){
#   if (simdata2_train$VL[i] <= 50) {simdata2_train$VL[i]=50}
# }

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
 write.table(simdata2, file="Uganda_SimData_SD0.R")
 write.table(simdata2_train, file="Uganda_SimData_train_SD0.R")
 write.table(simdata2_rev, file="Uganda_SimData_train_SD0_reverse.R")

