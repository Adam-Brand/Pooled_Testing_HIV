#==============================================================================
# FILENAME: sim_data_betas.R
# PROJECT: 	Pooled testing in HIV
# PURPOSE: get the estimated model betas from the simulated training sets
# AUTHOR: Adam Brand


# INPUT datasets: 4 datasets
#         - Uganda_SimData_train_SD0.R
#         - Uganda_SimData_train_SD0_reverse.R
#         - Uganda_SimData_train_SD1.0.R
#         - Uganda_SimData_train_SD1.0_reverse.R

# R VERSION: 3.6.1
#==============================================================================
#Notes: 

#### this program is to estimate the betas for the simulated data using ridge regression

#### we're going to estimate betas with the full model, ommitting the continuous variable
#### and ommitting the binary variable

#### we will estimate betas for each of our 2 datasets (SD=1.0 and SD=0)

#### IMPORTANT: for the paper we only used the estimated betas using the full model from the 
#               data with SD=1.0. When SD=0 the estimates were too perfect reflecting the lack 
#               of individual variability. We used 2 sets of betas, one from the full model, SD=1.0
#               from the standard training set, Uganda_SimData_train_SD1.0.R. The other used the 
#               reverse training set, Uganda_SimData_train_SD1.0_reverse.R.


# =============================================================================


library(tidyverse)
library(broom)
library(glmnet)
library(ggplot2)
library(useful)
library(sbrl)

# set working directory to location of data
setwd("")
simdata1.0 <- read.table("Uganda_SimData_train_SD1.0.R")
simdata0 <- read.table("Uganda_SimData_train_SD0.R")

simdata1.0_rev <- read.table("Uganda_SimData_train_SD1.0_reverse.R")
simdata0_rev <- read.table("Uganda_SimData_train_SD0_reverse.R")

# creating the interaction variabe for the full model
simdata1.0$int <- simdata1.0$adhere*simdata1.0$pf
simdata0$int <- simdata0$adhere*simdata0$pf

simdata1.0_rev$int <- simdata1.0_rev$adhere*simdata1.0_rev$pf
simdata0_rev$int <- simdata0_rev$adhere*simdata0_rev$pf

# setting the response variable
y1.0 <- simdata1.0$log.VL
y0 <- simdata0$log.VL

# setting the predictor matrix
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
## USED IN PAPER; no seed was set, so results may vary. The betas we used are recorded in the method_eval programs
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
## USED IN PAPER; no seed was set, so results may vary. The betas we used are recorded in the method_eval programs
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

















