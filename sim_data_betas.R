######### ###############################################

#### this program is to get estimated betas for the simulated data using ridge regression

#### we're going to estimate betas with the full model, ommitting the continuous variable
#### and ommitting the binary variable

#### we will estimate betas for each of our 2 datasets (SD=1.0 and SD=0)

library(tidyverse)
library(broom)
library(glmnet)
library(ggplot2)
library(useful)
library(sbrl)

setwd("C:/Users/Barny/Dropbox/KI_Project_4/Data/SimData")
simdata1.0 <- read.table("Uganda_SimData_train_SD1.0.R")
simdata0 <- read.table("Uganda_SimData_train_SD0.R")

simdata1.0_rev <- read.table("Uganda_SimData_train_SD1.0_reverse.R")
simdata0_rev <- read.table("Uganda_SimData_train_SD0_reverse.R")

simdata1.0$int <- simdata1.0$adhere*simdata1.0$pf
simdata0$int <- simdata0$adhere*simdata0$pf

simdata1.0_rev$int <- simdata1.0_rev$adhere*simdata1.0_rev$pf
simdata0_rev$int <- simdata0_rev$adhere*simdata0_rev$pf

y1.0 <- simdata1.0$log.VL
y0 <- simdata0$log.VL

x1.0 <- simdata1.0 %>% select(adhere, pf, int) %>% data.matrix()
x0 <- simdata0 %>% select(adhere, pf, int) %>% data.matrix()

x1.0_rev <- simdata1.0_rev %>% select(adhere, pf, int) %>% data.matrix()
x0_rev <- simdata0_rev %>% select(adhere, pf, int) %>% data.matrix()

#### setting a variety of lambdas
lambdas <- 10^seq(3, -4, by = -.1)

#### fitting ridge regression models with various lambdas
cv.fit1.0 <- cv.glmnet(x1.0, y1.0, alpha = 0, lambda = lambdas)
plot(cv.fit1.0)

cv.fit0 <- cv.glmnet(x0, y0, alpha = 0, lambda = lambdas)
plot(cv.fit0)

cv.fit1.0_rev <- cv.glmnet(x1.0_rev, y1.0, alpha = 0, lambda = lambdas)
plot(cv.fit1.0_rev)

cv.fit0_rev <- cv.glmnet(x0_rev, y0, alpha = 0, lambda = lambdas)
plot(cv.fit0_rev)

### setting my optimal lambda to the lambda with lowest MSE
opt.lambda1.0 <- cv.fit1.0$lambda.min
opt.lambda1.0

opt.lambda0 <- cv.fit0$lambda.min
opt.lambda0

opt.lambda1.0_rev <- cv.fit1.0_rev$lambda.min
opt.lambda1.0_rev

opt.lambda0_rev <- cv.fit0_rev$lambda.min
opt.lambda0_rev


### getting the betas from the optimal fit
opt.fit1.0 <- glmnet(x1.0, y1.0, alpha = 0, lambda = opt.lambda1.0)
summary(opt.fit1.0)
betas_ridge1.0 <- as.matrix(coef(opt.fit1.0))


### getting the betas from the optimal fit
### as may be expected, the ridge regression estimates are almost exactly the 
### actual betas when the SD=0, so we will not use these
### we will use the estimated betas from the SD=1.0 data, so the betas are variable
opt.fit0 <- glmnet(x0, y0, alpha = 0, lambda = opt.lambda0)
summary(opt.fit0)
betas_ridge0 <- as.matrix(coef(opt.fit0))


### getting the betas from the optimal fit from the reverse scenario
opt.fit1.0_rev <- glmnet(x1.0_rev, y1.0, alpha = 0, lambda = opt.lambda1.0_rev)
summary(opt.fit1.0_rev)
betas_ridge1.0_rev <- as.matrix(coef(opt.fit1.0_rev))

#### now, using only the SD=1.0 data, we will estimate betas with misspecified models

x1.0_no_bin <- simdata1.0 %>% select(adhere) %>% data.matrix()
x1.0_no_cont <- simdata1.0 %>% select(pf) %>% data.matrix()

x1.0_no_bin_rev <- simdata1.0_rev %>% select(adhere) %>% data.matrix()
x1.0_no_cont_rev <- simdata1.0_rev %>% select(pf) %>% data.matrix()

#### because the misspecified models have only 1 covariate, ridge doesn't apply
#### we will use OLS for the misspecified models

fit_ols_no_bin <- glm(y1.0 ~ x1.0_no_bin)
fit_ols_no_cont <- glm(y1.0 ~ x1.0_no_cont)

fit_ols_no_bin_rev <- glm(y1.0 ~ x1.0_no_bin_rev)
fit_ols_no_cont_rev <- glm(y1.0 ~ x1.0_no_cont_rev)

betas_ols_no_bin <- as.matrix(coef(fit_ols_no_bin))
betas_ols_no_cont <- as.matrix(coef(fit_ols_no_cont))

betas_ols_no_bin_rev <- as.matrix(coef(fit_ols_no_bin_rev))
betas_ols_no_cont_rev <- as.matrix(coef(fit_ols_no_cont_rev))


###### getting the classification rules using the sbrl paclage (Bayesian rule lists)
###### the idea is to get the rules using the training set, separate subjects into 3 risk categories
###### and then use the predicted risk categories as a covariate in the prediction model
###### we will use exclusively the SD=1.0 data

simdata_sbrl <- simdata1.0

for (i in 1:length(simdata_sbrl$VL)){
  if (simdata_sbrl$adhere[i] < 3){simdata_sbrl$adhcat[i]="low"}
  if (simdata_sbrl$adhere[i] >= 3 & simdata_sbrl$adhere[i] < 8){simdata_sbrl$adhcat[i]="medium"}
  if (simdata_sbrl$adhere[i] >= 8){simdata_sbrl$adhcat[i]="high"}
}

simdata_sbrl$pf <- as.factor(simdata_sbrl$pf)
simdata_sbrl$adhcat <- as.factor(simdata_sbrl$adhcat)
simdata_sbrl$label <- simdata_sbrl$VLcat

sim_train_sbrl <-simdata_sbrl %>% select(pf, adhcat, label) %>% data.frame()


sbrl_model <- sbrl(sim_train_sbrl, #### this is the training dataset to fit the model
                   iters=20000, #### number of iterations for each MCMC chain
                   pos_sign=c("vhigh","high"), #### the sign for the positive labels in the label column
                   neg_sign=c("low","medium"), #### the sign for the negative labels in the label column
                   rule_minlen=1, ### min number of cardinality for each rule (number of conditions for each rule)
                   rule_maxlen=3,   ### max number of cardinality for each rule
                   minsupport_pos=0.03, ### a number between 0 and 1, for the minimum percentage support for the positive observations
                   minsupport_neg=0.80, #### a number between 0 and 1, for the minimum percentage support for the negative obs
                   lambda=10, ### expected length of rule list
                   eta=2, #### expected cardinality of the rules in optimal rule list
                   nchain=25 ### number of MCMC chains
)

print(sbrl_model)
# The rules list is : 
#   If      {pf=1} (rule[9]) then positive probability = 0.74226804
# else  (default rule)  then positive probability = 0.50210224
















