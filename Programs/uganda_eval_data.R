#==============================================================================
# FILENAME: Uganda_eval_data.R
# PROJECT: 	Pooled Testing HIV
# PURPOSE: use the betas from the training set to predict VLs for the test set, and get dataset into
#          useable format for the method evaluation programs
# AUTHOR: Adam Brand

# CREATED:	20200428
# UPDATED: 	

# INPUT DATA: test_data2020-04-28.R				
# OUTPUT: 	

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

library(tidyverse)
library(broom)
library(glmnet)
library(ggplot2)
library(useful)
library(sbrl)
library(pROC)

setwd("")

test <- read.table("test_data.R")
train <- read.table("train_data.R")

##### getting the model betas and the predicted cutoff values from the train set

#### subsetting to get rid of the dub_fail outcome for modelling the single fail outcome
### we don't want to exclude records with NA for dub_fail
train2 <- train[,-which(names(train) == "dub_fail")]
train2 <- train2[complete.cases(train2),]

test2 <- test[,-which(names(test) == "dub_fail")]
test2 <- test2[complete.cases(test2),]


### defining the response variable
y <- train2$vllog
y_bin <- train2$fail_draw

y_bin_test <- test2$fail_draw

#### defining the predictor matrix
#### reference population is female, age 16 or younger, from a region not listed below, basewho=0, basereg was not one of the 4 common, with prev_fail=0,
#### whose last VL was over 6 months, and who was enrolled over 12 months ago
#### I am including only the 4 base regimens represented commonly in the training set (at least 50 records)
#### also only including regions with at least 50 records
x <- train2 %>% select(male, age35, age50, age65,
                            buyamba, kabira, kakuuto, kalisizo, kasaali, kasasa,
                            kibaale, kifamba, kyebe, lwamaggwa, 
                            lwanda, lyantonde, nabigasa,
                            basewho1, basewho2, basewho3, basewho4, 
                            reg_CBVEFV, reg_CBVNVP, reg_D4T3TCEFV,
                            reg_D4T3TCNVP,
                            prev_fail1, prev_fail2, prev_fail3,
                            bascd4log, basvllog, lastVLlog, lastVL_t, enroll_t) %>% data.matrix()

x_test <- test2 %>% select(male, age35, age50, age65,
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

### getting R squared for y and predicted y values - not good at 0.114; correlation at 0.34
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

# R squared for the ranks, even worse at 0.056
cor(act_pred_rnk$rank_order_y, act_pred_rnk$rank_order_y_pred)^2


##### getting the quantile cutoffs per the train data predicted VLs

cutoffs <- quantile(y_pred, probs=c(0.60, 0.90))


#### creating variable for the test set, which is the predicted VL value

test_pred <- predict(opt.fit, s=opt.lambda, newx=x_test)


#### putting together final testing dataset

test_final <- test2 %>% select(study_id, male, age35, age50, age65,
                           buyamba, kabira, kakuuto, kalisizo, kasaali, kasasa,
                           kibaale, kifamba, kyebe, lwamaggwa, 
                           lwanda, lyantonde, nabigasa,
                           basewho1, basewho2, basewho3, basewho4, 
                           reg_CBVEFV, reg_CBVNVP, reg_D4T3TCEFV,
                           reg_D4T3TCNVP,
                           prev_fail1, prev_fail2, prev_fail3,
                           bascd4log, basvllog, lastVLlog, lastVL_t, enroll_t, vl, vllog) %>% data.matrix()

test_final <- data.frame(cbind(test_final, test_pred))

test_final <- test_final %>% 
                        rename(
                        VL = vl,
                        log.VL = vllog,
                        predVL=X1
                        )


#### adding variable for predicted classification based on the prediction and cutoffs from the training set

for (i in 1:length(test_final$study_id)){
  if (test_final$predVL[i] > cutoffs[2]){test_final$pred_cat[i]="high"}
  if (test_final$predVL[i] <= cutoffs[2] & test_final$predVL[i] > cutoffs[1]){test_final$pred_cat[i]="mid"}
  if (test_final$predVL[i] <= cutoffs[1]){test_final$pred_cat[i]="low"}
}

for (i in 1:length(test_final$study_id)){
  if (test_final$VL[i] <= 50){test_final$VL[i]=50}
}

##### putting predVZL back on the actual scale

for (i in 1:length(test_final$study_id)){
  test_final$predVL[i] <- 10^test_final$predVL[i]
}


#### above we used the cutoffs at 90% and 60% from the predicted VLs from the training set to get cutoff log VL
#### values to classify subjects for the hypred method. We were shooting for 10%, 30%, 60%
#### using the cutoffs from the training set, we wind up with 5%, 45%, 50%

##### writing the final test dataset
### make sure to include the filepath for the location you want the data
write.table(test_final, "test_set_final.R")


########## getting the AUC in the test set

#### MODELLING binary outcome with ridge regression of the training set

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

#### getting the AUC in the training set
#### AUC is 0.78
auc(y_bin, probs)

probs_test <- opt.fit_bin %>% predict(newx=x_test, type="response")

probs_test <- as.vector(probs_test)

#### getting the AUC in the test set
#### AUC is 0.92
auc(y_bin_test, probs_test)




