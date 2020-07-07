#==============================================================================
# FILENAME:model_pred.R
# PROJECT: 	Pooled Testing HIV
# PURPOSE: use a variety of regression techniques to try to predict HIV VL using the training set(s) 
# AUTHOR: Adam Brand

# CREATED:	20200302

# INPUT DATA: Uganda training set (not in respository)				
	

# R VERSION: R version 3.6.1 (2019-07-05) -- "Action of the Toes"
#==============================================================================
# Notes:
#
#
#
#
#
#
#
# =============================================================================
#setwd("C:/Users/adabra/Dropbox/KI_Project_4/Programs")

library(tidyverse)
library(broom)
library(glmnet)
library(ggplot2)
library(useful)
library(sbrl)
library(pROC)

# set working directory to location of the Uganda training set
setwd("")

train.data2 <- read.table("train_data.R")

### code to find how many NAs are included for each variable
sum(is.na(train.data2$vllog))
sum(is.na(train.data2$sex))
sum(is.na(train.data2$ageyrs))
sum(is.na(train.data2$hub_name))
sum(is.na(train.data2$basewho))
sum(is.na(train.data2$basereg))
sum(is.na(train.data2$bascd4log)) #15 NAs
sum(is.na(train.data2$basvllog)) #282 NAs
sum(is.na(train.data2$prev_fail))
sum(is.na(train.data2$lastVLlog))
sum(is.na(train.data2$lastVL_t))
sum(is.na(train.data2$enroll_t))

### with less than 300 records with an NA, getting rid of any records with NA

#### subsetting to get rid of the dub_fail outcome for modelling the single fail outcome
### we don't want to exclude records with NA for dub_fail
train.data3 <- train.data2[,-which(names(train.data2) == "dub_fail")]
train.data3 <- train.data3[complete.cases(train.data3),]



### defining the response variable
y <- train.data3$vllog
y_bin <- train.data3$fail_draw

#### defining the predictor matrix
#### reference population is female, age 35 or younger, from a region not listed below, basewho=0, basereg was not one of the 4 common, with prev_fail=0,
#### whose last VL was over 6 months, and who was enrolled over 12 months ago
#### I am including only the 4 base regimens represented commonly in the training set (at least 50 records)
#### also only including regions with at least 50 records
x <- train.data3 %>% select(male, age50, age65,
                            buyamba, kabira, kakuuto, kalisizo, kasaali, kasasa,
                            kibaale, kifamba, kyebe, lwamaggwa, 
                            lwanda, lyantonde, nabigasa,
                            basewho1, basewho2, basewho3, basewho4, 
                            reg_CBVEFV, reg_CBVNVP, reg_D4T3TCEFV,
                            reg_D4T3TCNVP,
                            prev_fail1, prev_fail2, prev_fail3,
                            bascd4log, basvllog, lastVLlog, lastVL_t, enroll_t) %>% data.matrix()


### defining a sequence of tuning parameter, lambda, to fit the model over to find the best lambda
lambdas <- 10^seq(3, -4, by = -.1)


### MODELLING continuous VL fitting ridge regression models - one for each lambda; it looks like ridge regression doesn't like NAs

### fitting ridge regression with cross validation which evaluates performance via cross validation at each lambda 
### used to find our tuning parameter lambda
set.seed(1212)
cv.fit <- cv.glmnet(x, y, alpha = 0, lambda = lambdas)
plot(cv.fit)


### setting my optimal lambda to the lambda with lowest MSE
opt.lambda <- cv.fit$lambda.min
opt.lambda


### getting the betas from the optimal fit
opt.fit <- glmnet(x, y, alpha = 0, lambda = opt.lambda)
summary(opt.fit)
betas_ridge <- as.matrix(coef(opt.fit))


### predicting outcomes on the training set based on model with optimal lambda
y_pred <- predict(opt.fit, s=opt.lambda, newx=x)

### getting R squared for y and predicted y values - not good at 0.112; correlation at 0.34
cor(y, y_pred)^2


### scatterplot of predicted and actual VL values
plot(y_pred, y)


#### putting actual and predicted values into a matrix to look at correlation of true VL failures
act_pred <- data.frame(cbind(y,y_pred))

act_pred_rnk <- act_pred %>%
              mutate(rank_order_y = min_rank(desc(y))) %>%
              mutate(rank_order_y_pred = min_rank(desc(y_pred))) %>%
              arrange(rank_order_y)

#### plotting the actual and predicted ranks
plot(act_pred_rnk$rank_order_y, act_pred_rnk$rank_order_y_pred)

# R squared for the ranks, even worse at 0.053
cor(act_pred_rnk$rank_order_y, act_pred_rnk$rank_order_y_pred)^2

#### getting the auc for the continuous predictions - good at 0.84
auc(y, y_pred)


#### MODELLING binary outcome with ridge regression

#### use cross validation to fit a logistic model
cv.fit_bin <- cv.glmnet(x, y_bin, family="binomial", alpha = 0, lambda = lambdas)
plot(cv.fit_bin)

### grabing the optimal lambda
opt.lambda_bin <- cv.fit_bin$lambda.min
opt.lambda_bin

### grabbing the fit with best lambda
opt.fit_bin <- glmnet(x, y_bin, family="binomial", alpha = 0, lambda = opt.lambda_bin)
summary(opt.fit_bin)

### getting the predicted probabilities for each obs
probabilities <- opt.fit_bin %>% predict(newx=x, type="response")

probs <- as.vector(probabilities)

#### getting the AUC
#### AUC is 0.78
auc(y_bin, probs)

### rule to classify all probabilities over 50% to a failure, and no failure else
predicted_classes <- ifelse(probabilities > 0.5, 1, 0)

## we predict 16 failures
sum(predicted_classes)

### combining the actual outcomes with the predicted outcomes to compare
act_pred_bin <- data.frame(cbind(y_bin,predicted_classes, probabilities, y))
act_pred_bin <- act_pred_bin %>% 
              rename(
                    pred_resp = s0,
                    pred_prob = s0.1
                    )

### we predict 20 failures
sum(act_pred_bin$pred_resp)

#### we predict 15 out of 732 failures and predict 5 that weren't failures
check <- act_pred_bin[act_pred_bin$y_bin==1,]
sum(check$pred_resp==1)

check2 <- act_pred_bin[act_pred_bin$pred_resp==1,]


### comparing the ranks of the outcome with the ranks of probability
act_pred_rnk_bin <- act_pred_bin %>%
  mutate(rank_order_y = min_rank(desc(y))) %>%
  mutate(rank_order_y_pred = min_rank(desc(probabilities))) %>%
  arrange(rank_order_y)

### plotting the ranks of the continuous outcome versus the ranks of the predicted probabilties
plot(act_pred_rnk_bin$rank_order_y, act_pred_rnk_bin$rank_order_y_pred)
###correlation is .24 and R squared is .06
cor(act_pred_rnk_bin$rank_order_y, act_pred_rnk_bin$rank_order_y_pred)^2




##### MODELLING STANDARD OLS

fit_ols <- glm(y ~ x)
summary(fit_ols)
betas_ols <- as.matrix(coef(fit_ols))


#### comparing the OLS coefficients to the ridge regression coefficients
coeff <- cbind(betas_ridge, betas_ols)
coeff <- data.frame(coeff)
coeff$diff <- betas_ridge - betas_ols





### MODELLING continuous VL fitting lasso - one for each lambda; 

### fitting ridge regression with cross validation which evaluates performance via cross validation at each lambda 
### used to find our tuning parameter lambda
cv.fit_lasso <- cv.glmnet(x, y, alpha = 1, lambda = lambdas)
plot(cv.fit_lasso)


### setting my optimal lambda to the lambda with lowest MSE
cv.fit_lasso$lambda.min
opt.lambda <- cv.fit_lasso$lambda.1se
opt.lambda


### getting the betas from the optimal fit
opt.fit_lasso <- glmnet(x, y, alpha = 1, lambda = opt.lambda)
summary(opt.fit_lasso)
betas_lasso <- as.matrix(coef(opt.fit_lasso))


### predicting outcomes on the training set based on model with optimal lambda
y_pred_lasso <- predict(opt.fit_lasso, s=opt.lambda, newx=x)

### getting R squared for y and predicted y values - not good at 0.091; correlation at 0.34
cor(y, y_pred_lasso)^2

### scatterplot of predicted and actual VL values
plot(y_pred_lasso, y)


#####putting all of the beta coefficients for each model into 1 dataframe
betas <- data.frame(cbind(betas_ols, betas_ridge, betas_lasso))


#### putting actual and predicted values into a matrix to look at correlation of true VL failures
act_pred_lasso <- data.frame(cbind(y,y_pred_lasso))

act_pred_rnk_lasso <- act_pred_lasso %>%
  mutate(rank_order_y = min_rank(desc(y))) %>%
  mutate(rank_order_y_pred_lasso = min_rank(desc(y_pred_lasso))) %>%
  arrange(rank_order_y)

#### plotting the actual and predicted ranks
plot(act_pred_rnk_lasso$rank_order_y, act_pred_rnk_lasso$rank_order_y_pred_lasso)

# R squared for the ranks, even worse at 0.035
cor(act_pred_rnk_lasso$rank_order_y, act_pred_rnk_lasso$rank_order_y_pred_lasso)^2


#### MODELLING double failure as the outcome; putting dub_fail back into the matrix with Lasso
#### to serve as our outcome variable
#### we end up with fewer records, because some subjects have NA for the double failure outcome and
#### we exclude those
train.data4 <- train.data2[complete.cases(train.data2),]

### setting lambdas for the double failure models. Previous lmabdas did not converge
lambdas_dub <- 10^seq(5, -5, by = -.2)

##outcome variable
y_dub <- train.data4$dub_fail

## design matrix; we are omitting previous VL measure and time since previous VL measureto model the double failure
x_dub <- train.data4 %>% select(male, age35, age50, age65,
                            buyamba, kabira, kakuuto, kalisizo, kasaali, kasasa,
                            kibaale, kifamba, kyebe, lwamaggwa, 
                            lwanda, lyantonde, nabigasa,
                            basewho1, basewho2, basewho3, basewho4, 
                            reg_CBVEFV, reg_CBVNVP, reg_D4T3TCEFV,
                            reg_D4T3TCNVP,
                            prev_fail1, prev_fail2, prev_fail3,
                            bascd4log, basvllog, 
                            #lastVLlog, lastVL_t, 
                            enroll_t) %>% data.matrix()


### searching through to find the optimal lambda
cv.fit_lasso_dub <- cv.glmnet(x_dub, y_dub, family="binomial", alpha = 1, lambda = lambdas_dub)
plot(cv.fit_lasso_dub)


### setting my optimal lambda to within 1 SE of the lambda with lowest MSE
cv.fit_lasso_dub$lambda.min
opt.lambda_dub <- cv.fit_lasso_dub$lambda.1se
opt.lambda_dub


### getting the betas from the optimal fit
opt.fit_dub <- glmnet(x_dub, y_dub, alpha = 1, lambda = opt.lambda_dub)
summary(opt.fit_dub)
betas_lasso_dub <- as.matrix(coef(opt.fit_dub))


### getting the predicted probabilities for each obs
probabilities_dub <- opt.fit_dub %>% predict(newx=x_dub, type="response")

### rule to classify all probabilities over 50% to a failure, and no failure else
predicted.classes_dub <- ifelse(probabilities_dub > 0.5, 1, 0)

## we predict 0 failures
sum(predicted.classes_dub)

#### putting actual and predicted values into a matrix to look at correlation of true VL failures
act_pred_dub <- data.frame(cbind(y_dub,predicted.classes_dub, probabilities_dub))
act_pred_dub <- act_pred_dub %>% 
  rename(
    pred_resp = s0,
    pred_prob = s0.1
  )


#### Above: lasso using double failure as a binary outcome did not work. Separates patients
#### into 2 groups, neither of which are predicted to fail



##### MODELLING double failures using ridge regression
### searching through to find the optimal lambda
cv.fit_lasso_dub_rig <- cv.glmnet(x_dub, y_dub, family="binomial", alpha = 0, lambda = lambdas_dub)
plot(cv.fit_lasso_dub_rig)

### setting my optimal lambda to within 1 SE of the lambda with lowest MSE
opt.lambda_dub_rig <- cv.fit_lasso_dub$lambda.min
opt.lambda_dub_rig

### getting the betas from the optimal fit
opt.fit_dub_rig <- glmnet(x_dub, y_dub, alpha = 0, lambda = opt.lambda_dub_rig)
summary(opt.fit_dub_rig)
betas_ridge_dub <- as.matrix(coef(opt.fit_dub_rig))

### getting the predicted probabilities for each obs
probabilities_dub_rig <- opt.fit_dub_rig %>% predict(newx=x_dub, type="response")

### rule to classify all probabilities over 50% to a failure, and no failure else
predicted.classes_dub_rig <- ifelse(probabilities_dub_rig > 0.5, 1, 0)

## we predict 2 failures
sum(predicted.classes_dub_rig)

#### putting actual and predicted values into a matrix to look at correlation of true VL failures
act_pred_dub_rig <- data.frame(cbind(y_dub,predicted.classes_dub_rig, probabilities_dub_rig))
act_pred_dub_rig <- act_pred_dub_rig %>% 
  rename(
    pred_resp = s0,
    pred_prob = s0.1
  )

### comparing the ranks of the outcome with the ranks of probability
act_pred_rnk_bin_dub <- act_pred_dub_rig %>%
  mutate(rank_order_y = min_rank(desc(train.data4$vllog))) %>%
  mutate(rank_order_y_pred = min_rank(desc(probabilities_dub_rig))) %>%
  arrange(rank_order_y)


### plotting the ranks of the continuous outcome versus the ranks of the predicted probabilties
plot(act_pred_rnk_bin_dub$rank_order_y, act_pred_rnk_bin_dub$rank_order_y_pred)
###correlation is .24 and R squared is .06
cor(act_pred_rnk_bin_dub$rank_order_y, act_pred_rnk_bin_dub$rank_order_y_pred)^2

#### Above: ridge using double failure as a binary outcome did not work. Ranks of probabilities and
#### ranks of actual VL have very low R squared, 0.03


#############################################################################################################


