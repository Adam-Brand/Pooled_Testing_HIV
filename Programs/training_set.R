#==============================================================================
# FILENAME: training_set.R
# PROJECT: 	Pooled Testing HIV
# PURPOSE: derive the training set from clean_data 
# AUTHOR: Adam Brand

# CREATED:	20200302
# UPDATED: 	

# INPUT DATA: clean_data.R				
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
### set working directory to the location of the training_set.R program; location
### should also have the raw_data_clean.R program
setwd("")
source("raw_data_clean.R")
library(tidyverse)
library(broom)
library(glmnet)
library(ggplot2)
library(useful)
library(sbrl)
### running the data cleaning program


##### Program which splits the data into a training set and a train set

### deriving the train set; minimum 6 month tes date was 2004-12-03
### maximum 6 month train date was 2011-12-19. 

## declaring the 6 month window for the training set; grabbing all dates prior to 2010-04-01
date1 <- as.Date("2000-01-01")
date2 <- as.Date("2010-06-30")

##### selecting all VL measures within a 6 month window - there are no repeated subjects

# starting with the 6 month VLs; grabbing only observations with a 6 month VL in the date window defined by date1 and date2
# and have a basevl value
# below I am going to repeat these steps for values at 12 - 72 months, but will not repeat each comment
train.data6 <- mydata[(!is.na(mydata$visitdate6) & !is.na(mydata$vl_mth6) & mydata$visitdate6 <= date2 & mydata$visitdate6 >= date1),]
# creating variable for previous vl (prev_vl) andtime since prev_vl (prev_vl_t)
train.data6$prev_vl <- train.data6$basevl
train.data6$prev_vl_t <- 6
# keeping only records with a previous vl value
train.data6 <- train.data6[!is.na(train.data6$prev_vl),]

# creating variable prev_fail which is the number of previous treatment failures; because the last record is baseline, nobody will
# have an observed previous treatment failure
train.data6$prev_fail=0

# creating an indicator variabe for whether the VL is from visit date 6 or not
train.data6$visit <- 6

# creating a variable for VL which doesn't depend on visit (both raw and log10)
train.data6$vl <- train.data6$vl_mth6
train.data6$vllog <- train.data6$vl6

#creating a variable for experiencing a double failure at 6 months - defined as 2 consecutive, post-baseline failures
# by def'n there cannt be at 6 months, so all records set to 0
train.data6$dub_fail <- NA




# doing the same for 12 month values
train.data12 <- mydata[(!is.na(mydata$visitdate12) & !is.na(mydata$vl_mth12) & mydata$visitdate12 <= date2 & mydata$visitdate12 >= date1),]

# maiing sure subjects have at least 1 previous VL measure
for (i in 1:length(train.data12$study_id)){
  if (!is.na(train.data12$vl_mth6[i]))
  {train.data12$prev_vl[i]=train.data12$vl_mth6[i]
  train.data12$prev_vl_t[i]=6}
  else if (!is.na(train.data12$basevl[i]))
  {train.data12$prev_vl[i]=train.data12$basevl[i]
  train.data12$prev_vl_t[i]=12}
  else {train.data12$prev_vl[i]=NA}
}

train.data12 <- train.data12[!is.na(train.data12$prev_vl),]

for (i in 1:length(train.data12$study_id)) {
  if (!is.na(train.data12$basevl[i]) & train.data12$basevl[i] < 1000 & !is.na(train.data12$vl_mth6[i]) & train.data12$vl_mth6[i] >= 1000)
  {train.data12$prev_fail[i]=1}
  else {train.data12$prev_fail[i]=0}
}

train.data12$visit=12

train.data12$vl <- train.data12$vl_mth12
train.data12$vllog <- train.data12$vl12

### getting the double failure status at month 12
for (i in 1:length(train.data12$study_id)) {
  if (!is.na(train.data12$prev_vl[i]) & train.data12$prev_vl_t[i]<12 & train.data12$prev_vl[i]>=1000 & train.data12$vl_mth12[i]>=1000)
  {train.data12$dub_fail[i] = 1}
  else if (!is.na(train.data12$prev_vl[i]) & train.data12$prev_vl_t[i]<12 & (train.data12$prev_vl[i]<1000 | train.data12$vl_mth12[i]<1000))
  {train.data12$dub_fail[i] = 0}
  else {train.data12$dub_fail[i] = NA}
}

check <- train.data12[is.na(train.data12$dub_fail),]
check2 <- train.data12[!is.na(train.data12$dub_fail),]
check3 <- train.data12[train.data12$dub_fail==1 & !is.na(train.data12$dub_fail),]

# repeating for 18 month values
train.data18 <- mydata[(!is.na(mydata$visitdate18) & !is.na(mydata$vl_mth18) & mydata$visitdate18 <= date2 & mydata$visitdate18 >= date1),]
sum(duplicated(train.data18$study_id))

# maiking sure subjects have at least 1 previous VL measure
for (i in 1:length(train.data18$study_id)){
  if (!is.na(train.data18$vl_mth12[i]))
  {train.data18$prev_vl[i]=train.data18$vl_mth12[i]
  train.data18$prev_vl_t[i]=6}
  else if (!is.na(train.data18$vl_mth6[i]))
  {train.data18$prev_vl[i]=train.data18$vl_mth6[i]
  train.data18$prev_vl_t[i]=12}
  else if (!is.na(train.data18$basevl[i]))
  {train.data18$prev_vl[i]=train.data18$basevl[i]
  train.data18$prev_vl_t[i]=18}
  else {train.data18$prev_vl[i]=NA}
}

train.data18 <- train.data18[!is.na(train.data18$prev_vl),]

counter <- NULL
for (i in 1:length(train.data18$study_id)) {
  counter[i] <- 0
  # checking for previous failures at 6 months
  if (!is.na(train.data18$basevl[i]) & train.data18$basevl[i] < 1000 & !is.na(train.data18$vl_mth6[i]) & train.data18$vl_mth6[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 12 months
  if (((!is.na(train.data18$vl_mth6[i]) & train.data18$vl_mth6[i] < 1000) | 
       (is.na(train.data18$vl_mth6[i]) & !is.na(train.data18$basevl[i]) & train.data18$basevl[i] < 1000)) 
      & !is.na(train.data18$vl_mth12[i]) & train.data18$vl_mth12[i] >= 1000)
  {counter[i] = counter[i]+1}
  train.data18$prev_fail[i]=counter[i]
}

train.data18$visit=18

train.data18$vl <- train.data18$vl_mth18
train.data18$vllog <- train.data18$vl18

### getting the double failure status at month 18
for (i in 1:length(train.data18$study_id)) {
  if (!is.na(train.data18$prev_vl[i]) & train.data18$prev_vl_t[i]<18 & train.data18$prev_vl[i]>=1000 & train.data18$vl_mth18[i]>=1000)
  {train.data18$dub_fail[i] = 1}
  else if (!is.na(train.data18$prev_vl[i]) & train.data18$prev_vl_t[i]<18 & (train.data18$prev_vl[i]<1000 | train.data18$vl_mth18[i]<1000))
  {train.data18$dub_fail[i] = 0}
  else {train.data18$dub_fail[i] = NA}
}

check <- train.data18[is.na(train.data18$dub_fail),]
check2 <- train.data18[!is.na(train.data18$dub_fail),]
check3 <- train.data18[train.data18$dub_fail==1 & !is.na(train.data18$dub_fail),]


# repeating for 24 month values
train.data24 <- mydata[(!is.na(mydata$visitdate24) & !is.na(mydata$vl_mth24) & mydata$visitdate24 <= date2 & mydata$visitdate24 >= date1),]

# maiing sure subjects have at least 1 previous VL measure
for (i in 1:length(train.data24$study_id)){
  if (!is.na(train.data24$vl_mth18[i])) {
    train.data24$prev_vl[i]=train.data24$vl_mth18[i]
    train.data24$prev_vl_t[i]=6
  }
  else if (!is.na(train.data24$vl_mth12[i]))
  {train.data24$prev_vl[i]=train.data24$vl_mth12[i]
  train.data24$prev_vl_t[i]=12}
  else if (!is.na(train.data24$vl_mth6[i]))
  {train.data24$prev_vl[i]=train.data24$vl_mth6[i]
  train.data24$prev_vl_t[i]=18}
  else if (!is.na(train.data24$basevl[i]))
  {train.data24$prev_vl[i]=train.data24$basevl[i]
  train.data24$prev_vl_t[i]=24}
  else {train.data24$prev_vl[i]=NA}
}

train.data24 <- train.data24[!is.na(train.data24$prev_vl),]

counter <- NULL
for (i in 1:length(train.data24$study_id)) {
  counter[i] <- 0
  ### code to check for failures at 6 months
  if (!is.na(train.data24$basevl[i]) & train.data24$basevl[i] < 1000 & !is.na(train.data24$vl_mth6[i]) & train.data24$vl_mth6[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  ## checking for failures at 12 months
  if (((!is.na(train.data24$vl_mth6[i]) & train.data24$vl_mth6[i] < 1000) | 
       (is.na(train.data24$vl_mth6[i]) & !is.na(train.data24$basevl[i]) & train.data24$basevl[i] < 1000)) 
      & !is.na(train.data24$vl_mth12[i]) & train.data24$vl_mth12[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 18 months
  if (((!is.na(train.data24$vl_mth12[i]) & train.data24$vl_mth12[i] < 1000) |
       (is.na(train.data24$vl_mth12[i]) & !is.na(train.data24$vl_mth6[i]) & train.data24$vl_mth6[i] < 1000) |
       (is.na(train.data24$vl_mth12[i]) & is.na(train.data24$vl_mth6[i]) & !is.na(train.data24$basevl[i]) & train.data24$basevl[i] < 1000))
      & !is.na(train.data24$vl_mth18[i]) & train.data24$vl_mth18[i] >= 1000)
  {counter[i] = counter[i]+1}
  train.data24$prev_fail[i]=counter[i]
}

train.data24$visit=24

train.data24$vl <- train.data24$vl_mth24
train.data24$vllog <- train.data24$vl24

### getting the double failure status at month 24
for (i in 1:length(train.data24$study_id)) {
  if (!is.na(train.data24$prev_vl[i]) & train.data24$prev_vl_t[i]<24 & train.data24$prev_vl[i]>=1000 & train.data24$vl_mth24[i]>=1000)
  {train.data24$dub_fail[i] = 1}
  else if (!is.na(train.data24$prev_vl[i]) & train.data24$prev_vl_t[i]<24 & (train.data24$prev_vl[i]<1000 | train.data24$vl_mth24[i]<1000))
  {train.data24$dub_fail[i] = 0}
  else {train.data24$dub_fail[i] = NA}
}

check <- train.data24[is.na(train.data24$dub_fail),]
check2 <- train.data24[!is.na(train.data24$dub_fail),]
check3 <- train.data24[train.data24$dub_fail==1 & !is.na(train.data24$dub_fail),]

# repeating for 30 month values
train.data30 <- mydata[(!is.na(mydata$visitdate30) & !is.na(mydata$vl_mth30) & mydata$visitdate30 <= date2 & mydata$visitdate30 >= date1),]

# maiking sure subjects have at least 1 previous VL measure
for (i in 1:length(train.data30$study_id)){
  if (!is.na(train.data30$vl_mth24[i])){
    train.data30$prev_vl[i]=train.data30$vl_mth24[i]
    train.data30$prev_vl_t[i]=6
  }
  
  else if (!is.na(train.data30$vl_mth18[i])) {
    train.data30$prev_vl[i]=train.data30$vl_mth18[i]
    train.data30$prev_vl_t[i]=12
  }
  else if (!is.na(train.data30$vl_mth12[i]))
  {train.data30$prev_vl[i]=train.data30$vl_mth12[i]
  train.data30$prev_vl_t[i]=18}
  else if (!is.na(train.data30$vl_mth6[i]))
  {train.data30$prev_vl[i]=train.data30$vl_mth6[i]
  train.data30$prev_vl_t[i]=24}
  else if (!is.na(train.data30$basevl[i]))
  {train.data30$prev_vl[i]=train.data30$basevl[i]
  train.data30$prev_vl_t[i]=30}
  else {train.data30$prev_vl[i]=NA}
}

train.data30 <- train.data30[!is.na(train.data30$prev_vl),]

counter <- NULL
for (i in 1:length(train.data30$study_id)) {
  counter[i] <- 0
  ### code to check for failures at 6 months
  if (!is.na(train.data30$basevl[i]) & train.data30$basevl[i] < 1000 & !is.na(train.data30$vl_mth6[i]) & train.data30$vl_mth6[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  ## checking for failures at 12 months
  if (((!is.na(train.data30$vl_mth6[i]) & train.data30$vl_mth6[i] < 1000) | 
       (is.na(train.data30$vl_mth6[i]) & !is.na(train.data30$basevl[i]) & train.data30$basevl[i] < 1000)) 
      & !is.na(train.data30$vl_mth12[i]) & train.data30$vl_mth12[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 18 months
  if (((!is.na(train.data30$vl_mth12[i]) & train.data30$vl_mth12[i] < 1000) |
       (is.na(train.data30$vl_mth12[i]) & !is.na(train.data30$vl_mth6[i]) & train.data30$vl_mth6[i] < 1000) |
       (is.na(train.data30$vl_mth12[i]) & is.na(train.data30$vl_mth6[i]) & !is.na(train.data30$basevl[i]) & train.data30$basevl[i] < 1000))
      & !is.na(train.data30$vl_mth18[i]) & train.data30$vl_mth18[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 24 months
  if (((!is.na(train.data30$vl_mth18[i]) & train.data30$vl_mth18[i] < 1000) |
       (is.na(train.data30$vl_mth18[i]) & !is.na(train.data30$vl_mth12[i]) & train.data30$vl_mth12[i] < 1000) |
       (is.na(train.data30$vl_mth18[i]) & is.na(train.data30$vl_mth12[i]) & !is.na(train.data30$vl_mth6[i]) & train.data30$vl_mth6[i] < 1000) |
       (is.na(train.data30$vl_mth18[i]) & is.na(train.data30$vl_mth12[i]) & is.na(train.data30$vl_mth6[i]) & !is.na(train.data30$basevl[i]) & train.data30$basevl[i] < 1000))
      & !is.na(train.data30$vl_mth24[i]) & train.data30$vl_mth24[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  train.data30$prev_fail[i]=counter[i]
}

train.data30$visit=30

train.data30$vl <- train.data30$vl_mth30
train.data30$vllog <- train.data30$vl30

### getting the double failure status at month 30
for (i in 1:length(train.data30$study_id)) {
  if (!is.na(train.data30$prev_vl[i]) & train.data30$prev_vl_t[i]<30 & train.data30$prev_vl[i]>=1000 & train.data30$vl_mth30[i]>=1000)
  {train.data30$dub_fail[i] = 1}
  else if (!is.na(train.data30$prev_vl[i]) & train.data30$prev_vl_t[i]<30 & (train.data30$prev_vl[i]<1000 | train.data30$vl_mth30[i]<1000))
  {train.data30$dub_fail[i] = 0}
  else {train.data30$dub_fail[i] = NA}
}

check <- train.data30[is.na(train.data30$dub_fail),]
check2 <- train.data30[!is.na(train.data30$dub_fail),]
check3 <- train.data30[train.data30$dub_fail==1 & !is.na(train.data30$dub_fail),]

# repeating for 36 month values
train.data36 <- mydata[(!is.na(mydata$visitdate36) & !is.na(mydata$vl_mth36) & mydata$visitdate36 <= date2 & mydata$visitdate36 >= date1),]

# maiking sure subjects have at least 1 previous VL measure
for (i in 1:length(train.data36$study_id)){
  if (!is.na(train.data36$vl_mth30[i])){
    train.data36$prev_vl[i]=train.data36$vl_mth30[i]
    train.data36$prev_vl_t[i]=6
  }
  
  else if (!is.na(train.data36$vl_mth24[i])){
    train.data36$prev_vl[i]=train.data36$vl_mth24[i]
    train.data36$prev_vl_t[i]=12
  }
  
  else if (!is.na(train.data36$vl_mth18[i])) {
    train.data36$prev_vl[i]=train.data36$vl_mth18[i]
    train.data36$prev_vl_t[i]=18
  }
  else if (!is.na(train.data36$vl_mth12[i]))
  {train.data36$prev_vl[i]=train.data36$vl_mth12[i]
  train.data36$prev_vl_t[i]=24}
  else if (!is.na(train.data36$vl_mth6[i]))
  {train.data36$prev_vl[i]=train.data36$vl_mth6[i]
  train.data36$prev_vl_t[i]=30}
  else if (!is.na(train.data36$basevl[i]))
  {train.data36$prev_vl[i]=train.data36$basevl[i]
  train.data36$prev_vl_t[i]=36}
  else {train.data36$prev_vl[i]=NA}
}

train.data36 <- train.data36[!is.na(train.data36$prev_vl),]

counter <- NULL
for (i in 1:length(train.data36$study_id)) {
  counter[i] <- 0
  ### code to check for failures at 6 months
  if (!is.na(train.data36$basevl[i]) & train.data36$basevl[i] < 1000 & !is.na(train.data36$vl_mth6[i]) & train.data36$vl_mth6[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  ## checking for failures at 12 months
  if (((!is.na(train.data36$vl_mth6[i]) & train.data36$vl_mth6[i] < 1000) | 
       (is.na(train.data36$vl_mth6[i]) & !is.na(train.data36$basevl[i]) & train.data36$basevl[i] < 1000)) 
      & !is.na(train.data36$vl_mth12[i]) & train.data36$vl_mth12[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 18 months
  if (((!is.na(train.data36$vl_mth12[i]) & train.data36$vl_mth12[i] < 1000) |
       (is.na(train.data36$vl_mth12[i]) & !is.na(train.data36$vl_mth6[i]) & train.data36$vl_mth6[i] < 1000) |
       (is.na(train.data36$vl_mth12[i]) & is.na(train.data36$vl_mth6[i]) & !is.na(train.data36$basevl[i]) & train.data36$basevl[i] < 1000))
      & !is.na(train.data36$vl_mth18[i]) & train.data36$vl_mth18[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 24 months
  if (((!is.na(train.data36$vl_mth18[i]) & train.data36$vl_mth18[i] < 1000) |
       (is.na(train.data36$vl_mth18[i]) & !is.na(train.data36$vl_mth12[i]) & train.data36$vl_mth12[i] < 1000) |
       (is.na(train.data36$vl_mth18[i]) & is.na(train.data36$vl_mth12[i]) & !is.na(train.data36$vl_mth6[i]) & train.data36$vl_mth6[i] < 1000) |
       (is.na(train.data36$vl_mth18[i]) & is.na(train.data36$vl_mth12[i]) & is.na(train.data36$vl_mth6[i]) & !is.na(train.data36$basevl[i]) & train.data36$basevl[i] < 1000))
      & !is.na(train.data36$vl_mth24[i]) & train.data36$vl_mth24[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 30 months
  if (((!is.na(train.data36$vl_mth24[i]) & train.data36$vl_mth24[i] < 1000) |
       (is.na(train.data36$vl_mth24[i]) & !is.na(train.data36$vl_mth18[i]) & train.data36$vl_mth18[i] < 1000) |
       (is.na(train.data36$vl_mth24[i]) & is.na(train.data36$vl_mth18[i]) & !is.na(train.data36$vl_mth12[i]) & train.data36$vl_mth12[i] < 1000) |
       (is.na(train.data36$vl_mth24[i]) & is.na(train.data36$vl_mth18[i]) & is.na(train.data36$vl_mth12[i]) & !is.na(train.data36$vl_mth6[i]) & train.data36$vl_mth6[i] < 1000) |
       (is.na(train.data36$vl_mth24[i]) & is.na(train.data36$vl_mth18[i]) & is.na(train.data36$vl_mth12[i]) & is.na(train.data36$vl_mth6[i]) & !is.na(train.data36$basevl[i]) & train.data36$basevl[i] < 1000))
      & !is.na(train.data36$vl_mth30[i]) & train.data36$vl_mth30[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  train.data36$prev_fail[i]=counter[i]
}

train.data36$visit=36

train.data36$vl <- train.data36$vl_mth36
train.data36$vllog <- train.data36$vl36

### getting the double failure status at month 36
for (i in 1:length(train.data36$study_id)) {
  if (!is.na(train.data36$prev_vl[i]) & train.data36$prev_vl_t[i]<36 & train.data36$prev_vl[i]>=1000 & train.data36$vl_mth36[i]>=1000)
  {train.data36$dub_fail[i] = 1}
  else if (!is.na(train.data36$prev_vl[i]) & train.data36$prev_vl_t[i]<36 & (train.data36$prev_vl[i]<1000 | train.data36$vl_mth36[i]<1000))
  {train.data36$dub_fail[i] = 0}
  else {train.data36$dub_fail[i] = NA}
}

check <- train.data36[is.na(train.data36$dub_fail),]
check2 <- train.data36[!is.na(train.data36$dub_fail),]
check3 <- train.data36[train.data36$dub_fail==1 & !is.na(train.data36$dub_fail),]



# repeating for 42 month values
train.data42 <- mydata[(!is.na(mydata$visitdate42) & !is.na(mydata$vl_mth42) & mydata$visitdate42 <= date2 & mydata$visitdate42 >= date1),]

# maiking sure subjects have at least 1 previous VL measure
for (i in 1:length(train.data42$study_id)){
  if (!is.na(train.data42$vl_mth36[i])) {
    train.data42$prev_vl[i]=train.data42$vl_mth36[i]
    train.data42$prev_vl_t[i]=6
  }
  
  else if (!is.na(train.data42$vl_mth30[i])){
    train.data42$prev_vl[i]=train.data42$vl_mth30[i]
    train.data42$prev_vl_t[i]=12
  }
  
  else if (!is.na(train.data42$vl_mth24[i])){
    train.data42$prev_vl[i]=train.data42$vl_mth24[i]
    train.data42$prev_vl_t[i]=18
  }
  
  else if (!is.na(train.data42$vl_mth18[i])) {
    train.data42$prev_vl[i]=train.data42$vl_mth18[i]
    train.data42$prev_vl_t[i]=24
  }
  else if (!is.na(train.data42$vl_mth12[i]))
  {train.data42$prev_vl[i]=train.data42$vl_mth12[i]
  train.data42$prev_vl_t[i]=30}
  else if (!is.na(train.data42$vl_mth6[i]))
  {train.data42$prev_vl[i]=train.data42$vl_mth6[i]
  train.data42$prev_vl_t[i]=36}
  else if (!is.na(train.data42$basevl[i]))
  {train.data42$prev_vl[i]=train.data42$basevl[i]
  train.data42$prev_vl_t[i]=42}
  else {train.data42$prev_vl[i]=NA}
}

train.data42 <- train.data42[!is.na(train.data42$prev_vl),]

counter <- NULL
for (i in 1:length(train.data42$study_id)) {
  counter[i] <- 0
  ### code to check for failures at 6 months
  if (!is.na(train.data42$basevl[i]) & train.data42$basevl[i] < 1000 & !is.na(train.data42$vl_mth6[i]) & train.data42$vl_mth6[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  ## checking for failures at 12 months
  if (((!is.na(train.data42$vl_mth6[i]) & train.data42$vl_mth6[i] < 1000) | 
       (is.na(train.data42$vl_mth6[i]) & !is.na(train.data42$basevl[i]) & train.data42$basevl[i] < 1000)) 
      & !is.na(train.data42$vl_mth12[i]) & train.data42$vl_mth12[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 18 months
  if (((!is.na(train.data42$vl_mth12[i]) & train.data42$vl_mth12[i] < 1000) |
       (is.na(train.data42$vl_mth12[i]) & !is.na(train.data42$vl_mth6[i]) & train.data42$vl_mth6[i] < 1000) |
       (is.na(train.data42$vl_mth12[i]) & is.na(train.data42$vl_mth6[i]) & !is.na(train.data42$basevl[i]) & train.data42$basevl[i] < 1000))
      & !is.na(train.data42$vl_mth18[i]) & train.data42$vl_mth18[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 24 months
  if (((!is.na(train.data42$vl_mth18[i]) & train.data42$vl_mth18[i] < 1000) |
       (is.na(train.data42$vl_mth18[i]) & !is.na(train.data42$vl_mth12[i]) & train.data42$vl_mth12[i] < 1000) |
       (is.na(train.data42$vl_mth18[i]) & is.na(train.data42$vl_mth12[i]) & !is.na(train.data42$vl_mth6[i]) & train.data42$vl_mth6[i] < 1000) |
       (is.na(train.data42$vl_mth18[i]) & is.na(train.data42$vl_mth12[i]) & is.na(train.data42$vl_mth6[i]) & !is.na(train.data42$basevl[i]) & train.data42$basevl[i] < 1000))
      & !is.na(train.data42$vl_mth24[i]) & train.data42$vl_mth24[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 30 months
  if (((!is.na(train.data42$vl_mth24[i]) & train.data42$vl_mth24[i] < 1000) |
       (is.na(train.data42$vl_mth24[i]) & !is.na(train.data42$vl_mth18[i]) & train.data42$vl_mth18[i] < 1000) |
       (is.na(train.data42$vl_mth24[i]) & is.na(train.data42$vl_mth18[i]) & !is.na(train.data42$vl_mth12[i]) & train.data42$vl_mth12[i] < 1000) |
       (is.na(train.data42$vl_mth24[i]) & is.na(train.data42$vl_mth18[i]) & is.na(train.data42$vl_mth12[i]) & !is.na(train.data42$vl_mth6[i]) & train.data42$vl_mth6[i] < 1000) |
       (is.na(train.data42$vl_mth24[i]) & is.na(train.data42$vl_mth18[i]) & is.na(train.data42$vl_mth12[i]) & is.na(train.data42$vl_mth6[i]) & !is.na(train.data42$basevl[i]) & train.data42$basevl[i] < 1000))
      & !is.na(train.data42$vl_mth30[i]) & train.data42$vl_mth30[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 36 months
  if (((!is.na(train.data42$vl_mth30[i]) & train.data42$vl_mth30[i] < 1000) |
       (is.na(train.data42$vl_mth30[i]) & !is.na(train.data42$vl_mth24[i]) & train.data42$vl_mth24[i] < 1000) |
       (is.na(train.data42$vl_mth30[i]) & is.na(train.data42$vl_mth24[i]) & !is.na(train.data42$vl_mth24[i]) & train.data42$vl_mth24[i] < 1000) |
       (is.na(train.data42$vl_mth30[i]) & is.na(train.data42$vl_mth24[i]) & is.na(train.data42$vl_mth24[i]) & !is.na(train.data42$vl_mth12[i]) & train.data42$vl_mth12[i] < 1000) |
       (is.na(train.data42$vl_mth30[i]) & is.na(train.data42$vl_mth24[i]) & is.na(train.data42$vl_mth18[i]) & is.na(train.data42$vl_mth12[i]) & !is.na(train.data42$vl_mth6[i]) & train.data42$vl_mth6[i] < 1000) |
       (is.na(train.data42$vl_mth30[i]) & is.na(train.data42$vl_mth24[i]) & is.na(train.data42$vl_mth18[i]) & is.na(train.data42$vl_mth12[i]) & is.na(train.data42$vl_mth6[i]) & !is.na(train.data42$basevl[i]) & train.data42$basevl[i] < 1000))
      & !is.na(train.data42$vl_mth36[i]) & train.data42$vl_mth36[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  train.data42$prev_fail[i]=counter[i]
}

train.data42$visit=42

train.data42$vl <- train.data42$vl_mth42
train.data42$vllog <- train.data42$vl42

### getting the double failure status at month 42
for (i in 1:length(train.data42$study_id)) {
  if (!is.na(train.data42$prev_vl[i]) & train.data42$prev_vl_t[i]<42 & train.data42$prev_vl[i]>=1000 & train.data42$vl_mth42[i]>=1000)
  {train.data42$dub_fail[i] = 1}
  else if (!is.na(train.data42$prev_vl[i]) & train.data42$prev_vl_t[i]<42 & (train.data42$prev_vl[i]<1000 | train.data42$vl_mth42[i]<1000))
  {train.data42$dub_fail[i] = 0}
  else {train.data42$dub_fail[i] = NA}
}

check <- train.data42[is.na(train.data42$dub_fail),]
check2 <- train.data42[!is.na(train.data42$dub_fail),]
check3 <- train.data42[train.data42$dub_fail==1 & !is.na(train.data42$dub_fail),]


# repeating for 48 month values
train.data48 <- mydata[(!is.na(mydata$visitdate48) & !is.na(mydata$vl_mth48) & mydata$visitdate48 <= date2 & mydata$visitdate48 >= date1),]

# maiking sure subjects have at least 1 previous VL measure
for (i in 1:length(train.data48$study_id)){
  if (!is.na(train.data48$vl_mth42[i])){
    train.data48$prev_vl[i]=train.data48$vl_mth42[i]
    train.data48$prev_vl_t[i]=6
  }
  
  else if (!is.na(train.data48$vl_mth36[i])) {
    train.data48$prev_vl[i]=train.data48$vl_mth36[i]
    train.data48$prev_vl_t[i]=12
  }
  
  else if (!is.na(train.data48$vl_mth30[i])){
    train.data48$prev_vl[i]=train.data48$vl_mth30[i]
    train.data48$prev_vl_t[i]=18
  }
  
  else if (!is.na(train.data48$vl_mth24[i])){
    train.data48$prev_vl[i]=train.data48$vl_mth24[i]
    train.data48$prev_vl_t[i]=24
  }
  
  else if (!is.na(train.data48$vl_mth18[i])) {
    train.data48$prev_vl[i]=train.data48$vl_mth18[i]
    train.data48$prev_vl_t[i]=30
  }
  else if (!is.na(train.data48$vl_mth12[i]))
  {train.data48$prev_vl[i]=train.data48$vl_mth12[i]
  train.data48$prev_vl_t[i]=36}
  else if (!is.na(train.data48$vl_mth6[i]))
  {train.data48$prev_vl[i]=train.data48$vl_mth6[i]
  train.data48$prev_vl_t[i]=42}
  else if (!is.na(train.data48$basevl[i]))
  {train.data48$prev_vl[i]=train.data48$basevl[i]
  train.data48$prev_vl_t[i]=48}
  else {train.data48$prev_vl[i]=NA}
}

train.data48 <- train.data48[!is.na(train.data48$prev_vl),]

counter <- NULL
for (i in 1:length(train.data48$study_id)) {
  counter[i] <- 0
  ### code to check for failures at 6 months
  if (!is.na(train.data48$basevl[i]) & train.data48$basevl[i] < 1000 & !is.na(train.data48$vl_mth6[i]) & train.data48$vl_mth6[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  ## checking for failures at 12 months
  if (((!is.na(train.data48$vl_mth6[i]) & train.data48$vl_mth6[i] < 1000) | 
       (is.na(train.data48$vl_mth6[i]) & !is.na(train.data48$basevl[i]) & train.data48$basevl[i] < 1000)) 
      & !is.na(train.data48$vl_mth12[i]) & train.data48$vl_mth12[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 18 months
  if (((!is.na(train.data48$vl_mth12[i]) & train.data48$vl_mth12[i] < 1000) |
       (is.na(train.data48$vl_mth12[i]) & !is.na(train.data48$vl_mth6[i]) & train.data48$vl_mth6[i] < 1000) |
       (is.na(train.data48$vl_mth12[i]) & is.na(train.data48$vl_mth6[i]) & !is.na(train.data48$basevl[i]) & train.data48$basevl[i] < 1000))
      & !is.na(train.data48$vl_mth18[i]) & train.data48$vl_mth18[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 24 months
  if (((!is.na(train.data48$vl_mth18[i]) & train.data48$vl_mth18[i] < 1000) |
       (is.na(train.data48$vl_mth18[i]) & !is.na(train.data48$vl_mth12[i]) & train.data48$vl_mth12[i] < 1000) |
       (is.na(train.data48$vl_mth18[i]) & is.na(train.data48$vl_mth12[i]) & !is.na(train.data48$vl_mth6[i]) & train.data48$vl_mth6[i] < 1000) |
       (is.na(train.data48$vl_mth18[i]) & is.na(train.data48$vl_mth12[i]) & is.na(train.data48$vl_mth6[i]) & !is.na(train.data48$basevl[i]) & train.data48$basevl[i] < 1000))
      & !is.na(train.data48$vl_mth24[i]) & train.data48$vl_mth24[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 30 months
  if (((!is.na(train.data48$vl_mth24[i]) & train.data48$vl_mth24[i] < 1000) |
       (is.na(train.data48$vl_mth24[i]) & !is.na(train.data48$vl_mth18[i]) & train.data48$vl_mth18[i] < 1000) |
       (is.na(train.data48$vl_mth24[i]) & is.na(train.data48$vl_mth18[i]) & !is.na(train.data48$vl_mth12[i]) & train.data48$vl_mth12[i] < 1000) |
       (is.na(train.data48$vl_mth24[i]) & is.na(train.data48$vl_mth18[i]) & is.na(train.data48$vl_mth12[i]) & !is.na(train.data48$vl_mth6[i]) & train.data48$vl_mth6[i] < 1000) |
       (is.na(train.data48$vl_mth24[i]) & is.na(train.data48$vl_mth18[i]) & is.na(train.data48$vl_mth12[i]) & is.na(train.data48$vl_mth6[i]) & !is.na(train.data48$basevl[i]) & train.data48$basevl[i] < 1000))
      & !is.na(train.data48$vl_mth30[i]) & train.data48$vl_mth30[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 36 months
  if (((!is.na(train.data48$vl_mth30[i]) & train.data48$vl_mth30[i] < 1000) |
       (is.na(train.data48$vl_mth30[i]) & !is.na(train.data48$vl_mth24[i]) & train.data48$vl_mth24[i] < 1000) |
       (is.na(train.data48$vl_mth30[i]) & is.na(train.data48$vl_mth24[i]) & !is.na(train.data48$vl_mth24[i]) & train.data48$vl_mth24[i] < 1000) |
       (is.na(train.data48$vl_mth30[i]) & is.na(train.data48$vl_mth24[i]) & is.na(train.data48$vl_mth24[i]) & !is.na(train.data48$vl_mth12[i]) & train.data48$vl_mth12[i] < 1000) |
       (is.na(train.data48$vl_mth30[i]) & is.na(train.data48$vl_mth24[i]) & is.na(train.data48$vl_mth18[i]) & is.na(train.data48$vl_mth12[i]) & !is.na(train.data48$vl_mth6[i]) & train.data48$vl_mth6[i] < 1000) |
       (is.na(train.data48$vl_mth30[i]) & is.na(train.data48$vl_mth24[i]) & is.na(train.data48$vl_mth18[i]) & is.na(train.data48$vl_mth12[i]) & is.na(train.data48$vl_mth6[i]) & !is.na(train.data48$basevl[i]) & train.data48$basevl[i] < 1000))
      & !is.na(train.data48$vl_mth36[i]) & train.data48$vl_mth36[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 42 months
  if (((!is.na(train.data48$vl_mth36[i]) & train.data48$vl_mth36[i] < 1000) |
       (is.na(train.data48$vl_mth36[i]) & !is.na(train.data48$vl_mth30[i]) & train.data48$vl_mth30[i] < 1000) |
       (is.na(train.data48$vl_mth36[i]) & is.na(train.data48$vl_mth30[i]) & !is.na(train.data48$vl_mth24[i]) & train.data48$vl_mth24[i] < 1000) |
       (is.na(train.data48$vl_mth36[i]) & is.na(train.data48$vl_mth30[i]) & is.na(train.data48$vl_mth24[i]) & !is.na(train.data48$vl_mth18[i]) & train.data48$vl_mth18[i] < 1000) |
       (is.na(train.data48$vl_mth36[i]) & is.na(train.data48$vl_mth30[i]) & is.na(train.data48$vl_mth24[i]) & is.na(train.data48$vl_mth18[i]) & !is.na(train.data48$vl_mth12[i]) & train.data48$vl_mth12[i] < 1000) |
       (is.na(train.data48$vl_mth36[i]) & is.na(train.data48$vl_mth30[i]) & is.na(train.data48$vl_mth24[i]) & is.na(train.data48$vl_mth18[i]) & is.na(train.data48$vl_mth12[i]) & !is.na(train.data48$vl_mth6[i]) & train.data48$vl_mth6[i] < 1000) |
       (is.na(train.data48$vl_mth36[i]) & is.na(train.data48$vl_mth30[i]) & is.na(train.data48$vl_mth24[i]) & is.na(train.data48$vl_mth18[i]) & is.na(train.data48$vl_mth12[i]) & is.na(train.data48$vl_mth6[i]) & !is.na(train.data48$basevl[i]) & train.data48$basevl[i] < 1000))
      & !is.na(train.data48$vl_mth42[i]) & train.data48$vl_mth42[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  train.data48$prev_fail[i]=counter[i]
}

train.data48$visit=48

train.data48$vl <- train.data48$vl_mth48
train.data48$vllog <- train.data48$vl48

### getting the double failure status at month 48
for (i in 1:length(train.data48$study_id)) {
  if (!is.na(train.data48$prev_vl[i]) & train.data48$prev_vl_t[i]<48 & train.data48$prev_vl[i]>=1000 & train.data48$vl_mth48[i]>=1000)
  {train.data48$dub_fail[i] = 1}
  else if (!is.na(train.data48$prev_vl[i]) & train.data48$prev_vl_t[i]<48 & (train.data48$prev_vl[i]<1000 | train.data48$vl_mth48[i]<1000))
  {train.data48$dub_fail[i] = 0}
  else {train.data48$dub_fail[i] = NA}
}

check <- train.data48[is.na(train.data48$dub_fail),]
check2 <- train.data48[!is.na(train.data48$dub_fail),]
check3 <- train.data48[train.data48$dub_fail==1 & !is.na(train.data48$dub_fail),]




# repeating for 54 month values
train.data54 <- mydata[(!is.na(mydata$visitdate54) & !is.na(mydata$vl_mth54) & mydata$visitdate54 <= date2 & mydata$visitdate54 >= date1),]

# making sure subjects have at least 1 previous VL measure
for (i in 1:length(train.data54$study_id)){
  
  if (!is.na(train.data54$vl_mth48[i])){
    train.data54$prev_vl[i]=train.data54$vl_mth48[i]
    train.data54$prev_vl_t[i]=6
  }
  else if (!is.na(train.data54$vl_mth42[i])){
    train.data54$prev_vl[i]=train.data54$vl_mth42[i]
    train.data54$prev_vl_t[i]=12
  }
  
  else if (!is.na(train.data54$vl_mth36[i])) {
    train.data54$prev_vl[i]=train.data54$vl_mth36[i]
    train.data54$prev_vl_t[i]=18
  }
  
  else if (!is.na(train.data54$vl_mth30[i])){
    train.data54$prev_vl[i]=train.data54$vl_mth30[i]
    train.data54$prev_vl_t[i]=24
  }
  
  else if (!is.na(train.data54$vl_mth24[i])){
    train.data54$prev_vl[i]=train.data54$vl_mth24[i]
    train.data54$prev_vl_t[i]=30
  }
  
  else if (!is.na(train.data54$vl_mth18[i])) {
    train.data54$prev_vl[i]=train.data54$vl_mth18[i]
    train.data54$prev_vl_t[i]=36
  }
  else if (!is.na(train.data54$vl_mth12[i]))
  {train.data54$prev_vl[i]=train.data54$vl_mth12[i]
  train.data54$prev_vl_t[i]=42}
  else if (!is.na(train.data54$vl_mth6[i]))
  {train.data54$prev_vl[i]=train.data54$vl_mth6[i]
  train.data54$prev_vl_t[i]=48}
  else if (!is.na(train.data54$basevl[i]))
  {train.data54$prev_vl[i]=train.data54$basevl[i]
  train.data54$prev_vl_t[i]=54}
  else {train.data54$prev_vl[i]=NA}
}

train.data54 <- train.data54[!is.na(train.data54$prev_vl),]

counter <- NULL
for (i in 1:length(train.data54$study_id)) {
  counter[i] <- 0
  ### code to check for failures at 6 months
  if (!is.na(train.data54$basevl[i]) & train.data54$basevl[i] < 1000 & !is.na(train.data54$vl_mth6[i]) & train.data54$vl_mth6[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  ## checking for failures at 12 months
  if (((!is.na(train.data54$vl_mth6[i]) & train.data54$vl_mth6[i] < 1000) | 
       (is.na(train.data54$vl_mth6[i]) & !is.na(train.data54$basevl[i]) & train.data54$basevl[i] < 1000)) 
      & !is.na(train.data54$vl_mth12[i]) & train.data54$vl_mth12[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 18 months
  if (((!is.na(train.data54$vl_mth12[i]) & train.data54$vl_mth12[i] < 1000) |
       (is.na(train.data54$vl_mth12[i]) & !is.na(train.data54$vl_mth6[i]) & train.data54$vl_mth6[i] < 1000) |
       (is.na(train.data54$vl_mth12[i]) & is.na(train.data54$vl_mth6[i]) & !is.na(train.data54$basevl[i]) & train.data54$basevl[i] < 1000))
      & !is.na(train.data54$vl_mth18[i]) & train.data54$vl_mth18[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 24 months
  if (((!is.na(train.data54$vl_mth18[i]) & train.data54$vl_mth18[i] < 1000) |
       (is.na(train.data54$vl_mth18[i]) & !is.na(train.data54$vl_mth12[i]) & train.data54$vl_mth12[i] < 1000) |
       (is.na(train.data54$vl_mth18[i]) & is.na(train.data54$vl_mth12[i]) & !is.na(train.data54$vl_mth6[i]) & train.data54$vl_mth6[i] < 1000) |
       (is.na(train.data54$vl_mth18[i]) & is.na(train.data54$vl_mth12[i]) & is.na(train.data54$vl_mth6[i]) & !is.na(train.data54$basevl[i]) & train.data54$basevl[i] < 1000))
      & !is.na(train.data54$vl_mth24[i]) & train.data54$vl_mth24[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 30 months
  if (((!is.na(train.data54$vl_mth24[i]) & train.data54$vl_mth24[i] < 1000) |
       (is.na(train.data54$vl_mth24[i]) & !is.na(train.data54$vl_mth18[i]) & train.data54$vl_mth18[i] < 1000) |
       (is.na(train.data54$vl_mth24[i]) & is.na(train.data54$vl_mth18[i]) & !is.na(train.data54$vl_mth12[i]) & train.data54$vl_mth12[i] < 1000) |
       (is.na(train.data54$vl_mth24[i]) & is.na(train.data54$vl_mth18[i]) & is.na(train.data54$vl_mth12[i]) & !is.na(train.data54$vl_mth6[i]) & train.data54$vl_mth6[i] < 1000) |
       (is.na(train.data54$vl_mth24[i]) & is.na(train.data54$vl_mth18[i]) & is.na(train.data54$vl_mth12[i]) & is.na(train.data54$vl_mth6[i]) & !is.na(train.data54$basevl[i]) & train.data54$basevl[i] < 1000))
      & !is.na(train.data54$vl_mth30[i]) & train.data54$vl_mth30[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 36 months
  if (((!is.na(train.data54$vl_mth30[i]) & train.data54$vl_mth30[i] < 1000) |
       (is.na(train.data54$vl_mth30[i]) & !is.na(train.data54$vl_mth24[i]) & train.data54$vl_mth24[i] < 1000) |
       (is.na(train.data54$vl_mth30[i]) & is.na(train.data54$vl_mth24[i]) & !is.na(train.data54$vl_mth24[i]) & train.data54$vl_mth24[i] < 1000) |
       (is.na(train.data54$vl_mth30[i]) & is.na(train.data54$vl_mth24[i]) & is.na(train.data54$vl_mth24[i]) & !is.na(train.data54$vl_mth12[i]) & train.data54$vl_mth12[i] < 1000) |
       (is.na(train.data54$vl_mth30[i]) & is.na(train.data54$vl_mth24[i]) & is.na(train.data54$vl_mth18[i]) & is.na(train.data54$vl_mth12[i]) & !is.na(train.data54$vl_mth6[i]) & train.data54$vl_mth6[i] < 1000) |
       (is.na(train.data54$vl_mth30[i]) & is.na(train.data54$vl_mth24[i]) & is.na(train.data54$vl_mth18[i]) & is.na(train.data54$vl_mth12[i]) & is.na(train.data54$vl_mth6[i]) & !is.na(train.data54$basevl[i]) & train.data54$basevl[i] < 1000))
      & !is.na(train.data54$vl_mth36[i]) & train.data54$vl_mth36[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 42 months
  if (((!is.na(train.data54$vl_mth36[i]) & train.data54$vl_mth36[i] < 1000) |
       (is.na(train.data54$vl_mth36[i]) & !is.na(train.data54$vl_mth30[i]) & train.data54$vl_mth30[i] < 1000) |
       (is.na(train.data54$vl_mth36[i]) & is.na(train.data54$vl_mth30[i]) & !is.na(train.data54$vl_mth24[i]) & train.data54$vl_mth24[i] < 1000) |
       (is.na(train.data54$vl_mth36[i]) & is.na(train.data54$vl_mth30[i]) & is.na(train.data54$vl_mth24[i]) & !is.na(train.data54$vl_mth18[i]) & train.data54$vl_mth18[i] < 1000) |
       (is.na(train.data54$vl_mth36[i]) & is.na(train.data54$vl_mth30[i]) & is.na(train.data54$vl_mth24[i]) & is.na(train.data54$vl_mth18[i]) & !is.na(train.data54$vl_mth12[i]) & train.data54$vl_mth12[i] < 1000) |
       (is.na(train.data54$vl_mth36[i]) & is.na(train.data54$vl_mth30[i]) & is.na(train.data54$vl_mth24[i]) & is.na(train.data54$vl_mth18[i]) & is.na(train.data54$vl_mth12[i]) & !is.na(train.data54$vl_mth6[i]) & train.data54$vl_mth6[i] < 1000) |
       (is.na(train.data54$vl_mth36[i]) & is.na(train.data54$vl_mth30[i]) & is.na(train.data54$vl_mth24[i]) & is.na(train.data54$vl_mth18[i]) & is.na(train.data54$vl_mth12[i]) & is.na(train.data54$vl_mth6[i]) & !is.na(train.data54$basevl[i]) & train.data54$basevl[i] < 1000))
      & !is.na(train.data54$vl_mth42[i]) & train.data54$vl_mth42[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 48 months
  if (((!is.na(train.data54$vl_mth42[i]) & train.data54$vl_mth42[i] < 1000) |
       (is.na(train.data54$vl_mth42[i]) & !is.na(train.data54$vl_mth36[i]) & train.data54$vl_mth36[i] < 1000) |
       (is.na(train.data54$vl_mth42[i]) & is.na(train.data54$vl_mth36[i]) & !is.na(train.data54$vl_mth30[i]) & train.data54$vl_mth30[i] < 1000) |
       (is.na(train.data54$vl_mth42[i]) & is.na(train.data54$vl_mth36[i]) & is.na(train.data54$vl_mth30[i]) & !is.na(train.data54$vl_mth24[i]) & train.data54$vl_mth24[i] < 1000) |
       (is.na(train.data54$vl_mth42[i]) & is.na(train.data54$vl_mth36[i]) & is.na(train.data54$vl_mth30[i]) & is.na(train.data54$vl_mth24[i]) & !is.na(train.data54$vl_mth18[i]) & train.data54$vl_mth18[i] < 1000) |
       (is.na(train.data54$vl_mth42[i]) & is.na(train.data54$vl_mth36[i]) & is.na(train.data54$vl_mth30[i]) & is.na(train.data54$vl_mth24[i]) & is.na(train.data54$vl_mth18[i]) & !is.na(train.data54$vl_mth12[i]) & train.data54$vl_mth12[i] < 1000) |
       (is.na(train.data54$vl_mth42[i]) & is.na(train.data54$vl_mth36[i]) & is.na(train.data54$vl_mth30[i]) & is.na(train.data54$vl_mth24[i]) & is.na(train.data54$vl_mth18[i]) & is.na(train.data54$vl_mth12[i]) & !is.na(train.data54$vl_mth6[i]) & train.data54$vl_mth6[i] < 1000) | 
       (is.na(train.data54$vl_mth42[i]) & is.na(train.data54$vl_mth36[i]) & is.na(train.data54$vl_mth30[i]) & is.na(train.data54$vl_mth24[i]) & is.na(train.data54$vl_mth18[i]) & is.na(train.data54$vl_mth12[i]) & is.na(train.data54$vl_mth6[i]) & !is.na(train.data54$basevl[i]) & train.data54$basevl[i] < 1000))
      & !is.na(train.data54$vl_mth48[i]) & train.data54$vl_mth48[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  train.data54$prev_fail[i]=counter[i]
}

train.data54$visit=54

train.data54$vl <- train.data54$vl_mth54
train.data54$vllog <- train.data54$vl54

### getting the double failure status at month 54
for (i in 1:length(train.data54$study_id)) {
  if (!is.na(train.data54$prev_vl[i]) & train.data54$prev_vl_t[i]<54 & train.data54$prev_vl[i]>=1000 & train.data54$vl_mth54[i]>=1000)
  {train.data54$dub_fail[i] = 1}
  else if (!is.na(train.data54$prev_vl[i]) & train.data54$prev_vl_t[i]<54 & (train.data54$prev_vl[i]<1000 | train.data54$vl_mth54[i]<1000))
  {train.data54$dub_fail[i] = 0}
  else {train.data54$dub_fail[i] = NA}
}

check <- train.data54[is.na(train.data54$dub_fail),]
check2 <- train.data54[!is.na(train.data54$dub_fail),]
check3 <- train.data54[train.data54$dub_fail==1 & !is.na(train.data54$dub_fail),]


# repeating for 60 month values
train.data60 <- mydata[(!is.na(mydata$visitdate60) & !is.na(mydata$vl_mth60) & mydata$visitdate60 <= date2 & mydata$visitdate60 >= date1),]

# making sure subjects have at least 1 previous VL measure
for (i in 1:length(train.data60$study_id)){
  
  if (!is.na(train.data60$vl_mth54[i])){
    train.data60$prev_vl[i]=train.data60$vl_mth54[i]
    train.data60$prev_vl_t[i]=6
  }
  
  else if (!is.na(train.data60$vl_mth48[i])){
    train.data60$prev_vl[i]=train.data60$vl_mth48[i]
    train.data60$prev_vl_t[i]=12
  }
  else if (!is.na(train.data60$vl_mth42[i])){
    train.data60$prev_vl[i]=train.data60$vl_mth42[i]
    train.data60$prev_vl_t[i]=18
  }
  
  else if (!is.na(train.data60$vl_mth36[i])) {
    train.data60$prev_vl[i]=train.data60$vl_mth36[i]
    train.data60$prev_vl_t[i]=24
  }
  
  else if (!is.na(train.data60$vl_mth30[i])){
    train.data60$prev_vl[i]=train.data60$vl_mth30[i]
    train.data60$prev_vl_t[i]=30
  }
  
  else if (!is.na(train.data60$vl_mth24[i])){
    train.data60$prev_vl[i]=train.data60$vl_mth24[i]
    train.data60$prev_vl_t[i]=36
  }
  
  else if (!is.na(train.data60$vl_mth18[i])) {
    train.data60$prev_vl[i]=train.data60$vl_mth18[i]
    train.data60$prev_vl_t[i]=42
  }
  else if (!is.na(train.data60$vl_mth12[i]))
  {train.data60$prev_vl[i]=train.data60$vl_mth12[i]
  train.data60$prev_vl_t[i]=48}
  else if (!is.na(train.data60$vl_mth6[i]))
  {train.data60$prev_vl[i]=train.data60$vl_mth6[i]
  train.data60$prev_vl_t[i]=54}
  else if (!is.na(train.data60$basevl[i]))
  {train.data60$prev_vl[i]=train.data60$basevl[i]
  train.data60$prev_vl_t[i]=60}
  else {train.data60$prev_vl[i]=NA}
}

train.data60 <- train.data60[!is.na(train.data60$prev_vl),]

counter <- NULL
for (i in 1:length(train.data60$study_id)) {
  counter[i] <- 0
  ### code to check for failures at 6 months
  if (!is.na(train.data60$basevl[i]) & train.data60$basevl[i] < 1000 & !is.na(train.data60$vl_mth6[i]) & train.data60$vl_mth6[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  ## checking for failures at 12 months
  if (((!is.na(train.data60$vl_mth6[i]) & train.data60$vl_mth6[i] < 1000) | 
       (is.na(train.data60$vl_mth6[i]) & !is.na(train.data60$basevl[i]) & train.data60$basevl[i] < 1000)) 
      & !is.na(train.data60$vl_mth12[i]) & train.data60$vl_mth12[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 18 months
  if (((!is.na(train.data60$vl_mth12[i]) & train.data60$vl_mth12[i] < 1000) |
       (is.na(train.data60$vl_mth12[i]) & !is.na(train.data60$vl_mth6[i]) & train.data60$vl_mth6[i] < 1000) |
       (is.na(train.data60$vl_mth12[i]) & is.na(train.data60$vl_mth6[i]) & !is.na(train.data60$basevl[i]) & train.data60$basevl[i] < 1000))
      & !is.na(train.data60$vl_mth18[i]) & train.data60$vl_mth18[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 24 months
  if (((!is.na(train.data60$vl_mth18[i]) & train.data60$vl_mth18[i] < 1000) |
       (is.na(train.data60$vl_mth18[i]) & !is.na(train.data60$vl_mth12[i]) & train.data60$vl_mth12[i] < 1000) |
       (is.na(train.data60$vl_mth18[i]) & is.na(train.data60$vl_mth12[i]) & !is.na(train.data60$vl_mth6[i]) & train.data60$vl_mth6[i] < 1000) |
       (is.na(train.data60$vl_mth18[i]) & is.na(train.data60$vl_mth12[i]) & is.na(train.data60$vl_mth6[i]) & !is.na(train.data60$basevl[i]) & train.data60$basevl[i] < 1000))
      & !is.na(train.data60$vl_mth24[i]) & train.data60$vl_mth24[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 30 months
  if (((!is.na(train.data60$vl_mth24[i]) & train.data60$vl_mth24[i] < 1000) |
       (is.na(train.data60$vl_mth24[i]) & !is.na(train.data60$vl_mth18[i]) & train.data60$vl_mth18[i] < 1000) |
       (is.na(train.data60$vl_mth24[i]) & is.na(train.data60$vl_mth18[i]) & !is.na(train.data60$vl_mth12[i]) & train.data60$vl_mth12[i] < 1000) |
       (is.na(train.data60$vl_mth24[i]) & is.na(train.data60$vl_mth18[i]) & is.na(train.data60$vl_mth12[i]) & !is.na(train.data60$vl_mth6[i]) & train.data60$vl_mth6[i] < 1000) |
       (is.na(train.data60$vl_mth24[i]) & is.na(train.data60$vl_mth18[i]) & is.na(train.data60$vl_mth12[i]) & is.na(train.data60$vl_mth6[i]) & !is.na(train.data60$basevl[i]) & train.data60$basevl[i] < 1000))
      & !is.na(train.data60$vl_mth30[i]) & train.data60$vl_mth30[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 36 months
  if (((!is.na(train.data60$vl_mth30[i]) & train.data60$vl_mth30[i] < 1000) |
       (is.na(train.data60$vl_mth30[i]) & !is.na(train.data60$vl_mth24[i]) & train.data60$vl_mth24[i] < 1000) |
       (is.na(train.data60$vl_mth30[i]) & is.na(train.data60$vl_mth24[i]) & !is.na(train.data60$vl_mth24[i]) & train.data60$vl_mth24[i] < 1000) |
       (is.na(train.data60$vl_mth30[i]) & is.na(train.data60$vl_mth24[i]) & is.na(train.data60$vl_mth24[i]) & !is.na(train.data60$vl_mth12[i]) & train.data60$vl_mth12[i] < 1000) |
       (is.na(train.data60$vl_mth30[i]) & is.na(train.data60$vl_mth24[i]) & is.na(train.data60$vl_mth18[i]) & is.na(train.data60$vl_mth12[i]) & !is.na(train.data60$vl_mth6[i]) & train.data60$vl_mth6[i] < 1000) |
       (is.na(train.data60$vl_mth30[i]) & is.na(train.data60$vl_mth24[i]) & is.na(train.data60$vl_mth18[i]) & is.na(train.data60$vl_mth12[i]) & is.na(train.data60$vl_mth6[i]) & !is.na(train.data60$basevl[i]) & train.data60$basevl[i] < 1000))
      & !is.na(train.data60$vl_mth36[i]) & train.data60$vl_mth36[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 42 months
  if (((!is.na(train.data60$vl_mth36[i]) & train.data60$vl_mth36[i] < 1000) |
       (is.na(train.data60$vl_mth36[i]) & !is.na(train.data60$vl_mth30[i]) & train.data60$vl_mth30[i] < 1000) |
       (is.na(train.data60$vl_mth36[i]) & is.na(train.data60$vl_mth30[i]) & !is.na(train.data60$vl_mth24[i]) & train.data60$vl_mth24[i] < 1000) |
       (is.na(train.data60$vl_mth36[i]) & is.na(train.data60$vl_mth30[i]) & is.na(train.data60$vl_mth24[i]) & !is.na(train.data60$vl_mth18[i]) & train.data60$vl_mth18[i] < 1000) |
       (is.na(train.data60$vl_mth36[i]) & is.na(train.data60$vl_mth30[i]) & is.na(train.data60$vl_mth24[i]) & is.na(train.data60$vl_mth18[i]) & !is.na(train.data60$vl_mth12[i]) & train.data60$vl_mth12[i] < 1000) |
       (is.na(train.data60$vl_mth36[i]) & is.na(train.data60$vl_mth30[i]) & is.na(train.data60$vl_mth24[i]) & is.na(train.data60$vl_mth18[i]) & is.na(train.data60$vl_mth12[i]) & !is.na(train.data60$vl_mth6[i]) & train.data60$vl_mth6[i] < 1000) |
       (is.na(train.data60$vl_mth36[i]) & is.na(train.data60$vl_mth30[i]) & is.na(train.data60$vl_mth24[i]) & is.na(train.data60$vl_mth18[i]) & is.na(train.data60$vl_mth12[i]) & is.na(train.data60$vl_mth6[i]) & !is.na(train.data60$basevl[i]) & train.data60$basevl[i] < 1000))
      & !is.na(train.data60$vl_mth42[i]) & train.data60$vl_mth42[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 48 months
  if (((!is.na(train.data60$vl_mth42[i]) & train.data60$vl_mth42[i] < 1000) |
       (is.na(train.data60$vl_mth42[i]) & !is.na(train.data60$vl_mth36[i]) & train.data60$vl_mth36[i] < 1000) |
       (is.na(train.data60$vl_mth42[i]) & is.na(train.data60$vl_mth36[i]) & !is.na(train.data60$vl_mth30[i]) & train.data60$vl_mth30[i] < 1000) |
       (is.na(train.data60$vl_mth42[i]) & is.na(train.data60$vl_mth36[i]) & is.na(train.data60$vl_mth30[i]) & !is.na(train.data60$vl_mth24[i]) & train.data60$vl_mth24[i] < 1000) |
       (is.na(train.data60$vl_mth42[i]) & is.na(train.data60$vl_mth36[i]) & is.na(train.data60$vl_mth30[i]) & is.na(train.data60$vl_mth24[i]) & !is.na(train.data60$vl_mth18[i]) & train.data60$vl_mth18[i] < 1000) |
       (is.na(train.data60$vl_mth42[i]) & is.na(train.data60$vl_mth36[i]) & is.na(train.data60$vl_mth30[i]) & is.na(train.data60$vl_mth24[i]) & is.na(train.data60$vl_mth18[i]) & !is.na(train.data60$vl_mth12[i]) & train.data60$vl_mth12[i] < 1000) |
       (is.na(train.data60$vl_mth42[i]) & is.na(train.data60$vl_mth36[i]) & is.na(train.data60$vl_mth30[i]) & is.na(train.data60$vl_mth24[i]) & is.na(train.data60$vl_mth18[i]) & is.na(train.data60$vl_mth12[i]) & !is.na(train.data60$vl_mth6[i]) & train.data60$vl_mth6[i] < 1000) | 
       (is.na(train.data60$vl_mth42[i]) & is.na(train.data60$vl_mth36[i]) & is.na(train.data60$vl_mth30[i]) & is.na(train.data60$vl_mth24[i]) & is.na(train.data60$vl_mth18[i]) & is.na(train.data60$vl_mth12[i]) & is.na(train.data60$vl_mth6[i]) & !is.na(train.data60$basevl[i]) & train.data60$basevl[i] < 1000))
      & !is.na(train.data60$vl_mth48[i]) & train.data60$vl_mth48[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 54 months
  if (((!is.na(train.data60$vl_mth48[i]) & train.data60$vl_mth48[i] < 1000) |
       (is.na(train.data60$vl_mth48[i]) & !is.na(train.data60$vl_mth42[i]) & train.data60$vl_mth42[i] < 1000) |
       (is.na(train.data60$vl_mth48[i]) & is.na(train.data60$vl_mth42[i]) & !is.na(train.data60$vl_mth36[i]) & train.data60$vl_mth36[i] < 1000) |
       (is.na(train.data60$vl_mth48[i]) & is.na(train.data60$vl_mth42[i]) & is.na(train.data60$vl_mth36[i]) & !is.na(train.data60$vl_mth30[i]) & train.data60$vl_mth30[i] < 1000) |
       (is.na(train.data60$vl_mth48[i]) & is.na(train.data60$vl_mth42[i]) & is.na(train.data60$vl_mth36[i]) & is.na(train.data60$vl_mth30[i]) & !is.na(train.data60$vl_mth24[i]) & train.data60$vl_mth24[i] < 1000) |
       (is.na(train.data60$vl_mth48[i]) & is.na(train.data60$vl_mth42[i]) & is.na(train.data60$vl_mth36[i]) & is.na(train.data60$vl_mth30[i]) & is.na(train.data60$vl_mth24[i]) & !is.na(train.data60$vl_mth18[i]) & train.data60$vl_mth18[i] < 1000) |
       (is.na(train.data60$vl_mth48[i]) & is.na(train.data60$vl_mth42[i]) & is.na(train.data60$vl_mth36[i]) & is.na(train.data60$vl_mth30[i]) & is.na(train.data60$vl_mth24[i]) & is.na(train.data60$vl_mth18[i]) & !is.na(train.data60$vl_mth12[i]) & train.data60$vl_mth12[i] < 1000) | 
       (is.na(train.data60$vl_mth48[i]) & is.na(train.data60$vl_mth42[i]) & is.na(train.data60$vl_mth36[i]) & is.na(train.data60$vl_mth30[i]) & is.na(train.data60$vl_mth24[i]) & is.na(train.data60$vl_mth18[i]) & is.na(train.data60$vl_mth12[i]) & !is.na(train.data60$vl_mth6[i]) & train.data60$vl_mth6[i] < 1000) |
       (is.na(train.data60$vl_mth48[i]) & is.na(train.data60$vl_mth42[i]) & is.na(train.data60$vl_mth36[i]) & is.na(train.data60$vl_mth30[i]) & is.na(train.data60$vl_mth24[i]) & is.na(train.data60$vl_mth18[i]) & is.na(train.data60$vl_mth12[i]) & is.na(train.data60$vl_mth6[i]) & !is.na(train.data60$basevl[i]) & train.data60$basevl[i] < 1000))
      & !is.na(train.data60$vl_mth54[i]) & train.data60$vl_mth54[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  train.data60$prev_fail[i]=counter[i]
}

train.data60$visit=60

train.data60$vl <- train.data60$vl_mth60
train.data60$vllog <- train.data60$vl60

### getting the double failure status at month 60
for (i in 1:length(train.data60$study_id)) {
  if (!is.na(train.data60$prev_vl[i]) & train.data60$prev_vl_t[i]<60 & train.data60$prev_vl[i]>=1000 & train.data60$vl_mth60[i]>=1000)
  {train.data60$dub_fail[i] = 1}
  else if (!is.na(train.data60$prev_vl[i]) & train.data60$prev_vl_t[i]<60 & (train.data60$prev_vl[i]<1000 | train.data60$vl_mth60[i]<1000))
  {train.data60$dub_fail[i] = 0}
  else {train.data60$dub_fail[i] = NA}
}

check <- train.data60[is.na(train.data60$dub_fail),]
check2 <- train.data60[!is.na(train.data60$dub_fail),]
check3 <- train.data60[train.data60$dub_fail==1 & !is.na(train.data60$dub_fail),]


# repeating for 66 month values
train.data66 <- mydata[(!is.na(mydata$visitdate66) & !is.na(mydata$vl_mth66) & mydata$visitdate66 <= date2 & mydata$visitdate66 >= date1),]

# making sure subjects have at least 1 previous VL measure
for (i in 1:length(train.data66$study_id)){
  
  if (!is.na(train.data66$vl_mth60[i])){
    train.data66$prev_vl[i]=train.data66$vl_mth60[i]
    train.data66$prev_vl_t[i]=6
  }
  
  else if (!is.na(train.data66$vl_mth54[i])){
    train.data66$prev_vl[i]=train.data66$vl_mth54[i]
    train.data66$prev_vl_t[i]=12
  }
  
  else if (!is.na(train.data66$vl_mth48[i])){
    train.data66$prev_vl[i]=train.data66$vl_mth48[i]
    train.data66$prev_vl_t[i]=18
  }
  else if (!is.na(train.data66$vl_mth42[i])){
    train.data66$prev_vl[i]=train.data66$vl_mth42[i]
    train.data66$prev_vl_t[i]=24
  }
  
  else if (!is.na(train.data66$vl_mth36[i])) {
    train.data66$prev_vl[i]=train.data66$vl_mth36[i]
    train.data66$prev_vl_t[i]=30
  }
  
  else if (!is.na(train.data66$vl_mth30[i])){
    train.data66$prev_vl[i]=train.data66$vl_mth30[i]
    train.data66$prev_vl_t[i]=36
  }
  
  else if (!is.na(train.data66$vl_mth24[i])){
    train.data66$prev_vl[i]=train.data66$vl_mth24[i]
    train.data66$prev_vl_t[i]=42
  }
  
  else if (!is.na(train.data66$vl_mth18[i])) {
    train.data66$prev_vl[i]=train.data66$vl_mth18[i]
    train.data66$prev_vl_t[i]=48
  }
  else if (!is.na(train.data66$vl_mth12[i]))
  {train.data66$prev_vl[i]=train.data66$vl_mth12[i]
  train.data66$prev_vl_t[i]=54}
  else if (!is.na(train.data66$vl_mth6[i]))
  {train.data66$prev_vl[i]=train.data66$vl_mth6[i]
  train.data66$prev_vl_t[i]=60}
  else if (!is.na(train.data66$basevl[i]))
  {train.data66$prev_vl[i]=train.data66$basevl[i]
  train.data66$prev_vl_t[i]=66}
  else {train.data66$prev_vl[i]=NA}
}

train.data66 <- train.data66[!is.na(train.data66$prev_vl),]

counter <- NULL
for (i in 1:length(train.data66$study_id)) {
  counter[i] <- 0
  ### code to check for failures at 6 months
  if (!is.na(train.data66$basevl[i]) & train.data66$basevl[i] < 1000 & !is.na(train.data66$vl_mth6[i]) & train.data66$vl_mth6[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  ## checking for failures at 12 months
  if (((!is.na(train.data66$vl_mth6[i]) & train.data66$vl_mth6[i] < 1000) | 
       (is.na(train.data66$vl_mth6[i]) & !is.na(train.data66$basevl[i]) & train.data66$basevl[i] < 1000)) 
      & !is.na(train.data66$vl_mth12[i]) & train.data66$vl_mth12[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 18 months
  if (((!is.na(train.data66$vl_mth12[i]) & train.data66$vl_mth12[i] < 1000) |
       (is.na(train.data66$vl_mth12[i]) & !is.na(train.data66$vl_mth6[i]) & train.data66$vl_mth6[i] < 1000) |
       (is.na(train.data66$vl_mth12[i]) & is.na(train.data66$vl_mth6[i]) & !is.na(train.data66$basevl[i]) & train.data66$basevl[i] < 1000))
      & !is.na(train.data66$vl_mth18[i]) & train.data66$vl_mth18[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 24 months
  if (((!is.na(train.data66$vl_mth18[i]) & train.data66$vl_mth18[i] < 1000) |
       (is.na(train.data66$vl_mth18[i]) & !is.na(train.data66$vl_mth12[i]) & train.data66$vl_mth12[i] < 1000) |
       (is.na(train.data66$vl_mth18[i]) & is.na(train.data66$vl_mth12[i]) & !is.na(train.data66$vl_mth6[i]) & train.data66$vl_mth6[i] < 1000) |
       (is.na(train.data66$vl_mth18[i]) & is.na(train.data66$vl_mth12[i]) & is.na(train.data66$vl_mth6[i]) & !is.na(train.data66$basevl[i]) & train.data66$basevl[i] < 1000))
      & !is.na(train.data66$vl_mth24[i]) & train.data66$vl_mth24[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 30 months
  if (((!is.na(train.data66$vl_mth24[i]) & train.data66$vl_mth24[i] < 1000) |
       (is.na(train.data66$vl_mth24[i]) & !is.na(train.data66$vl_mth18[i]) & train.data66$vl_mth18[i] < 1000) |
       (is.na(train.data66$vl_mth24[i]) & is.na(train.data66$vl_mth18[i]) & !is.na(train.data66$vl_mth12[i]) & train.data66$vl_mth12[i] < 1000) |
       (is.na(train.data66$vl_mth24[i]) & is.na(train.data66$vl_mth18[i]) & is.na(train.data66$vl_mth12[i]) & !is.na(train.data66$vl_mth6[i]) & train.data66$vl_mth6[i] < 1000) |
       (is.na(train.data66$vl_mth24[i]) & is.na(train.data66$vl_mth18[i]) & is.na(train.data66$vl_mth12[i]) & is.na(train.data66$vl_mth6[i]) & !is.na(train.data66$basevl[i]) & train.data66$basevl[i] < 1000))
      & !is.na(train.data66$vl_mth30[i]) & train.data66$vl_mth30[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 36 months
  if (((!is.na(train.data66$vl_mth30[i]) & train.data66$vl_mth30[i] < 1000) |
       (is.na(train.data66$vl_mth30[i]) & !is.na(train.data66$vl_mth24[i]) & train.data66$vl_mth24[i] < 1000) |
       (is.na(train.data66$vl_mth30[i]) & is.na(train.data66$vl_mth24[i]) & !is.na(train.data66$vl_mth24[i]) & train.data66$vl_mth24[i] < 1000) |
       (is.na(train.data66$vl_mth30[i]) & is.na(train.data66$vl_mth24[i]) & is.na(train.data66$vl_mth24[i]) & !is.na(train.data66$vl_mth12[i]) & train.data66$vl_mth12[i] < 1000) |
       (is.na(train.data66$vl_mth30[i]) & is.na(train.data66$vl_mth24[i]) & is.na(train.data66$vl_mth18[i]) & is.na(train.data66$vl_mth12[i]) & !is.na(train.data66$vl_mth6[i]) & train.data66$vl_mth6[i] < 1000) |
       (is.na(train.data66$vl_mth30[i]) & is.na(train.data66$vl_mth24[i]) & is.na(train.data66$vl_mth18[i]) & is.na(train.data66$vl_mth12[i]) & is.na(train.data66$vl_mth6[i]) & !is.na(train.data66$basevl[i]) & train.data66$basevl[i] < 1000))
      & !is.na(train.data66$vl_mth36[i]) & train.data66$vl_mth36[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 42 months
  if (((!is.na(train.data66$vl_mth36[i]) & train.data66$vl_mth36[i] < 1000) |
       (is.na(train.data66$vl_mth36[i]) & !is.na(train.data66$vl_mth30[i]) & train.data66$vl_mth30[i] < 1000) |
       (is.na(train.data66$vl_mth36[i]) & is.na(train.data66$vl_mth30[i]) & !is.na(train.data66$vl_mth24[i]) & train.data66$vl_mth24[i] < 1000) |
       (is.na(train.data66$vl_mth36[i]) & is.na(train.data66$vl_mth30[i]) & is.na(train.data66$vl_mth24[i]) & !is.na(train.data66$vl_mth18[i]) & train.data66$vl_mth18[i] < 1000) |
       (is.na(train.data66$vl_mth36[i]) & is.na(train.data66$vl_mth30[i]) & is.na(train.data66$vl_mth24[i]) & is.na(train.data66$vl_mth18[i]) & !is.na(train.data66$vl_mth12[i]) & train.data66$vl_mth12[i] < 1000) |
       (is.na(train.data66$vl_mth36[i]) & is.na(train.data66$vl_mth30[i]) & is.na(train.data66$vl_mth24[i]) & is.na(train.data66$vl_mth18[i]) & is.na(train.data66$vl_mth12[i]) & !is.na(train.data66$vl_mth6[i]) & train.data66$vl_mth6[i] < 1000) |
       (is.na(train.data66$vl_mth36[i]) & is.na(train.data66$vl_mth30[i]) & is.na(train.data66$vl_mth24[i]) & is.na(train.data66$vl_mth18[i]) & is.na(train.data66$vl_mth12[i]) & is.na(train.data66$vl_mth6[i]) & !is.na(train.data66$basevl[i]) & train.data66$basevl[i] < 1000))
      & !is.na(train.data66$vl_mth42[i]) & train.data66$vl_mth42[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 48 months
  if (((!is.na(train.data66$vl_mth42[i]) & train.data66$vl_mth42[i] < 1000) |
       (is.na(train.data66$vl_mth42[i]) & !is.na(train.data66$vl_mth36[i]) & train.data66$vl_mth36[i] < 1000) |
       (is.na(train.data66$vl_mth42[i]) & is.na(train.data66$vl_mth36[i]) & !is.na(train.data66$vl_mth30[i]) & train.data66$vl_mth30[i] < 1000) |
       (is.na(train.data66$vl_mth42[i]) & is.na(train.data66$vl_mth36[i]) & is.na(train.data66$vl_mth30[i]) & !is.na(train.data66$vl_mth24[i]) & train.data66$vl_mth24[i] < 1000) |
       (is.na(train.data66$vl_mth42[i]) & is.na(train.data66$vl_mth36[i]) & is.na(train.data66$vl_mth30[i]) & is.na(train.data66$vl_mth24[i]) & !is.na(train.data66$vl_mth18[i]) & train.data66$vl_mth18[i] < 1000) |
       (is.na(train.data66$vl_mth42[i]) & is.na(train.data66$vl_mth36[i]) & is.na(train.data66$vl_mth30[i]) & is.na(train.data66$vl_mth24[i]) & is.na(train.data66$vl_mth18[i]) & !is.na(train.data66$vl_mth12[i]) & train.data66$vl_mth12[i] < 1000) |
       (is.na(train.data66$vl_mth42[i]) & is.na(train.data66$vl_mth36[i]) & is.na(train.data66$vl_mth30[i]) & is.na(train.data66$vl_mth24[i]) & is.na(train.data66$vl_mth18[i]) & is.na(train.data66$vl_mth12[i]) & !is.na(train.data66$vl_mth6[i]) & train.data66$vl_mth6[i] < 1000) | 
       (is.na(train.data66$vl_mth42[i]) & is.na(train.data66$vl_mth36[i]) & is.na(train.data66$vl_mth30[i]) & is.na(train.data66$vl_mth24[i]) & is.na(train.data66$vl_mth18[i]) & is.na(train.data66$vl_mth12[i]) & is.na(train.data66$vl_mth6[i]) & !is.na(train.data66$basevl[i]) & train.data66$basevl[i] < 1000))
      & !is.na(train.data66$vl_mth48[i]) & train.data66$vl_mth48[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 54 months
  if (((!is.na(train.data66$vl_mth48[i]) & train.data66$vl_mth48[i] < 1000) |
       (is.na(train.data66$vl_mth48[i]) & !is.na(train.data66$vl_mth42[i]) & train.data66$vl_mth42[i] < 1000) |
       (is.na(train.data66$vl_mth48[i]) & is.na(train.data66$vl_mth42[i]) & !is.na(train.data66$vl_mth36[i]) & train.data66$vl_mth36[i] < 1000) |
       (is.na(train.data66$vl_mth48[i]) & is.na(train.data66$vl_mth42[i]) & is.na(train.data66$vl_mth36[i]) & !is.na(train.data66$vl_mth30[i]) & train.data66$vl_mth30[i] < 1000) |
       (is.na(train.data66$vl_mth48[i]) & is.na(train.data66$vl_mth42[i]) & is.na(train.data66$vl_mth36[i]) & is.na(train.data66$vl_mth30[i]) & !is.na(train.data66$vl_mth24[i]) & train.data66$vl_mth24[i] < 1000) |
       (is.na(train.data66$vl_mth48[i]) & is.na(train.data66$vl_mth42[i]) & is.na(train.data66$vl_mth36[i]) & is.na(train.data66$vl_mth30[i]) & is.na(train.data66$vl_mth24[i]) & !is.na(train.data66$vl_mth18[i]) & train.data66$vl_mth18[i] < 1000) |
       (is.na(train.data66$vl_mth48[i]) & is.na(train.data66$vl_mth42[i]) & is.na(train.data66$vl_mth36[i]) & is.na(train.data66$vl_mth30[i]) & is.na(train.data66$vl_mth24[i]) & is.na(train.data66$vl_mth18[i]) & !is.na(train.data66$vl_mth12[i]) & train.data66$vl_mth12[i] < 1000) | 
       (is.na(train.data66$vl_mth48[i]) & is.na(train.data66$vl_mth42[i]) & is.na(train.data66$vl_mth36[i]) & is.na(train.data66$vl_mth30[i]) & is.na(train.data66$vl_mth24[i]) & is.na(train.data66$vl_mth18[i]) & is.na(train.data66$vl_mth12[i]) & !is.na(train.data66$vl_mth6[i]) & train.data66$vl_mth6[i] < 1000) |
       (is.na(train.data66$vl_mth48[i]) & is.na(train.data66$vl_mth42[i]) & is.na(train.data66$vl_mth36[i]) & is.na(train.data66$vl_mth30[i]) & is.na(train.data66$vl_mth24[i]) & is.na(train.data66$vl_mth18[i]) & is.na(train.data66$vl_mth12[i]) & is.na(train.data66$vl_mth6[i]) & !is.na(train.data66$basevl[i]) & train.data66$basevl[i] < 1000))
      & !is.na(train.data66$vl_mth54[i]) & train.data66$vl_mth54[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 60 months
  if (((!is.na(train.data66$vl_mth54[i]) & train.data66$vl_mth54[i] < 1000) |
       (is.na(train.data66$vl_mth54[i]) & !is.na(train.data66$vl_mth48[i]) & train.data66$vl_mth48[i] < 1000) |
       (is.na(train.data66$vl_mth54[i]) & is.na(train.data66$vl_mth48[i]) & !is.na(train.data66$vl_mth42[i]) & train.data66$vl_mth42[i] < 1000) |
       (is.na(train.data66$vl_mth54[i]) & is.na(train.data66$vl_mth48[i]) & is.na(train.data66$vl_mth42[i]) & !is.na(train.data66$vl_mth36[i]) & train.data66$vl_mth36[i] < 1000) |
       (is.na(train.data66$vl_mth54[i]) & is.na(train.data66$vl_mth48[i]) & is.na(train.data66$vl_mth42[i]) & is.na(train.data66$vl_mth36[i]) & !is.na(train.data66$vl_mth30[i]) & train.data66$vl_mth30[i] < 1000) |
       (is.na(train.data66$vl_mth54[i]) & is.na(train.data66$vl_mth48[i]) & is.na(train.data66$vl_mth42[i]) & is.na(train.data66$vl_mth36[i]) & is.na(train.data66$vl_mth30[i]) & !is.na(train.data66$vl_mth24[i]) & train.data66$vl_mth24[i] < 1000) |
       (is.na(train.data66$vl_mth54[i]) & is.na(train.data66$vl_mth48[i]) & is.na(train.data66$vl_mth42[i]) & is.na(train.data66$vl_mth36[i]) & is.na(train.data66$vl_mth30[i]) & is.na(train.data66$vl_mth24[i]) & !is.na(train.data66$vl_mth18[i]) & train.data66$vl_mth18[i] < 1000) | 
       (is.na(train.data66$vl_mth54[i]) & is.na(train.data66$vl_mth48[i]) & is.na(train.data66$vl_mth42[i]) & is.na(train.data66$vl_mth36[i]) & is.na(train.data66$vl_mth30[i]) & is.na(train.data66$vl_mth24[i]) & is.na(train.data66$vl_mth18[i]) & !is.na(train.data66$vl_mth12[i]) & train.data66$vl_mth12[i] < 1000) |
       (is.na(train.data66$vl_mth54[i]) & is.na(train.data66$vl_mth48[i]) & is.na(train.data66$vl_mth42[i]) & is.na(train.data66$vl_mth36[i]) & is.na(train.data66$vl_mth30[i]) & is.na(train.data66$vl_mth24[i]) & is.na(train.data66$vl_mth18[i]) & is.na(train.data66$vl_mth12[i]) & !is.na(train.data66$vl_mth6[i]) & train.data66$vl_mth6[i] < 1000) |
       (is.na(train.data66$vl_mth54[i]) & is.na(train.data66$vl_mth48[i]) & is.na(train.data66$vl_mth42[i]) & is.na(train.data66$vl_mth36[i]) & is.na(train.data66$vl_mth30[i]) & is.na(train.data66$vl_mth24[i]) & is.na(train.data66$vl_mth18[i]) & is.na(train.data66$vl_mth12[i]) & is.na(train.data66$vl_mth6[i]) & !is.na(train.data66$basevl[i]) & train.data66$basevl[i] < 1000))
      & !is.na(train.data66$vl_mth60[i]) & train.data66$vl_mth60[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  train.data66$prev_fail[i]=counter[i]
}

train.data66$visit=66

train.data66$vl <- train.data66$vl_mth66
train.data66$vllog <- train.data66$vl66

### getting the double failure status at month 66
for (i in 1:length(train.data66$study_id)) {
  if (!is.na(train.data66$prev_vl[i]) & train.data66$prev_vl_t[i]<66 & train.data66$prev_vl[i]>=1000 & train.data66$vl_mth66[i]>=1000)
  {train.data66$dub_fail[i] = 1}
  else if (!is.na(train.data66$prev_vl[i]) & train.data66$prev_vl_t[i]<66 & (train.data66$prev_vl[i]<1000 | train.data66$vl_mth66[i]<1000))
  {train.data66$dub_fail[i] = 0}
  else {train.data66$dub_fail[i] = NA}
}

check <- train.data66[is.na(train.data66$dub_fail),]
check2 <- train.data66[!is.na(train.data66$dub_fail),]
check3 <- train.data66[train.data66$dub_fail==1 & !is.na(train.data66$dub_fail),]



# repeating for 72 month values
train.data72 <- mydata[(!is.na(mydata$visitdate72) & !is.na(mydata$vl_mth72) & mydata$visitdate72 <= date2 & mydata$visitdate72 >= date1),]

# making sure subjects have at least 1 previous VL measure
for (i in 1:length(train.data72$study_id)){
  
  if (!is.na(train.data72$vl_mth66[i])){
    train.data72$prev_vl[i]=train.data72$vl_mth66[i]
    train.data72$prev_vl_t[i]=6
  }
  
  else if (!is.na(train.data72$vl_mth60[i])){
    train.data72$prev_vl[i]=train.data72$vl_mth60[i]
    train.data72$prev_vl_t[i]=12
  }
  
  else if (!is.na(train.data72$vl_mth54[i])){
    train.data72$prev_vl[i]=train.data72$vl_mth54[i]
    train.data72$prev_vl_t[i]=18
  }
  
  else if (!is.na(train.data72$vl_mth48[i])){
    train.data72$prev_vl[i]=train.data72$vl_mth48[i]
    train.data72$prev_vl_t[i]=24
  }
  else if (!is.na(train.data72$vl_mth42[i])){
    train.data72$prev_vl[i]=train.data72$vl_mth42[i]
    train.data72$prev_vl_t[i]=30
  }
  
  else if (!is.na(train.data72$vl_mth36[i])) {
    train.data72$prev_vl[i]=train.data72$vl_mth36[i]
    train.data72$prev_vl_t[i]=36
  }
  
  else if (!is.na(train.data72$vl_mth30[i])){
    train.data72$prev_vl[i]=train.data72$vl_mth30[i]
    train.data72$prev_vl_t[i]=42
  }
  
  else if (!is.na(train.data72$vl_mth24[i])){
    train.data72$prev_vl[i]=train.data72$vl_mth24[i]
    train.data72$prev_vl_t[i]=48
  }
  
  else if (!is.na(train.data72$vl_mth18[i])) {
    train.data72$prev_vl[i]=train.data72$vl_mth18[i]
    train.data72$prev_vl_t[i]=54
  }
  else if (!is.na(train.data72$vl_mth12[i]))
  {train.data72$prev_vl[i]=train.data72$vl_mth12[i]
  train.data72$prev_vl_t[i]=60}
  else if (!is.na(train.data72$vl_mth6[i]))
  {train.data72$prev_vl[i]=train.data72$vl_mth6[i]
  train.data72$prev_vl_t[i]=66}
  else if (!is.na(train.data72$basevl[i]))
  {train.data72$prev_vl[i]=train.data72$basevl[i]
  train.data72$prev_vl_t[i]=72}
  else {train.data72$prev_vl[i]=NA}
}

train.data72 <- train.data72[!is.na(train.data72$prev_vl),]

counter <- NULL
for (i in 1:length(train.data72$study_id)) {
  counter[i] <- 0
  ### code to check for failures at 6 months
  if (!is.na(train.data72$basevl[i]) & train.data72$basevl[i] < 1000 & !is.na(train.data72$vl_mth6[i]) & train.data72$vl_mth6[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  ## checking for failures at 12 months
  if (((!is.na(train.data72$vl_mth6[i]) & train.data72$vl_mth6[i] < 1000) | 
       (is.na(train.data72$vl_mth6[i]) & !is.na(train.data72$basevl[i]) & train.data72$basevl[i] < 1000)) 
      & !is.na(train.data72$vl_mth12[i]) & train.data72$vl_mth12[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 18 months
  if (((!is.na(train.data72$vl_mth12[i]) & train.data72$vl_mth12[i] < 1000) |
       (is.na(train.data72$vl_mth12[i]) & !is.na(train.data72$vl_mth6[i]) & train.data72$vl_mth6[i] < 1000) |
       (is.na(train.data72$vl_mth12[i]) & is.na(train.data72$vl_mth6[i]) & !is.na(train.data72$basevl[i]) & train.data72$basevl[i] < 1000))
      & !is.na(train.data72$vl_mth18[i]) & train.data72$vl_mth18[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 24 months
  if (((!is.na(train.data72$vl_mth18[i]) & train.data72$vl_mth18[i] < 1000) |
       (is.na(train.data72$vl_mth18[i]) & !is.na(train.data72$vl_mth12[i]) & train.data72$vl_mth12[i] < 1000) |
       (is.na(train.data72$vl_mth18[i]) & is.na(train.data72$vl_mth12[i]) & !is.na(train.data72$vl_mth6[i]) & train.data72$vl_mth6[i] < 1000) |
       (is.na(train.data72$vl_mth18[i]) & is.na(train.data72$vl_mth12[i]) & is.na(train.data72$vl_mth6[i]) & !is.na(train.data72$basevl[i]) & train.data72$basevl[i] < 1000))
      & !is.na(train.data72$vl_mth24[i]) & train.data72$vl_mth24[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 30 months
  if (((!is.na(train.data72$vl_mth24[i]) & train.data72$vl_mth24[i] < 1000) |
       (is.na(train.data72$vl_mth24[i]) & !is.na(train.data72$vl_mth18[i]) & train.data72$vl_mth18[i] < 1000) |
       (is.na(train.data72$vl_mth24[i]) & is.na(train.data72$vl_mth18[i]) & !is.na(train.data72$vl_mth12[i]) & train.data72$vl_mth12[i] < 1000) |
       (is.na(train.data72$vl_mth24[i]) & is.na(train.data72$vl_mth18[i]) & is.na(train.data72$vl_mth12[i]) & !is.na(train.data72$vl_mth6[i]) & train.data72$vl_mth6[i] < 1000) |
       (is.na(train.data72$vl_mth24[i]) & is.na(train.data72$vl_mth18[i]) & is.na(train.data72$vl_mth12[i]) & is.na(train.data72$vl_mth6[i]) & !is.na(train.data72$basevl[i]) & train.data72$basevl[i] < 1000))
      & !is.na(train.data72$vl_mth30[i]) & train.data72$vl_mth30[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 36 months
  if (((!is.na(train.data72$vl_mth30[i]) & train.data72$vl_mth30[i] < 1000) |
       (is.na(train.data72$vl_mth30[i]) & !is.na(train.data72$vl_mth24[i]) & train.data72$vl_mth24[i] < 1000) |
       (is.na(train.data72$vl_mth30[i]) & is.na(train.data72$vl_mth24[i]) & !is.na(train.data72$vl_mth24[i]) & train.data72$vl_mth24[i] < 1000) |
       (is.na(train.data72$vl_mth30[i]) & is.na(train.data72$vl_mth24[i]) & is.na(train.data72$vl_mth24[i]) & !is.na(train.data72$vl_mth12[i]) & train.data72$vl_mth12[i] < 1000) |
       (is.na(train.data72$vl_mth30[i]) & is.na(train.data72$vl_mth24[i]) & is.na(train.data72$vl_mth18[i]) & is.na(train.data72$vl_mth12[i]) & !is.na(train.data72$vl_mth6[i]) & train.data72$vl_mth6[i] < 1000) |
       (is.na(train.data72$vl_mth30[i]) & is.na(train.data72$vl_mth24[i]) & is.na(train.data72$vl_mth18[i]) & is.na(train.data72$vl_mth12[i]) & is.na(train.data72$vl_mth6[i]) & !is.na(train.data72$basevl[i]) & train.data72$basevl[i] < 1000))
      & !is.na(train.data72$vl_mth36[i]) & train.data72$vl_mth36[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 42 months
  if (((!is.na(train.data72$vl_mth36[i]) & train.data72$vl_mth36[i] < 1000) |
       (is.na(train.data72$vl_mth36[i]) & !is.na(train.data72$vl_mth30[i]) & train.data72$vl_mth30[i] < 1000) |
       (is.na(train.data72$vl_mth36[i]) & is.na(train.data72$vl_mth30[i]) & !is.na(train.data72$vl_mth24[i]) & train.data72$vl_mth24[i] < 1000) |
       (is.na(train.data72$vl_mth36[i]) & is.na(train.data72$vl_mth30[i]) & is.na(train.data72$vl_mth24[i]) & !is.na(train.data72$vl_mth18[i]) & train.data72$vl_mth18[i] < 1000) |
       (is.na(train.data72$vl_mth36[i]) & is.na(train.data72$vl_mth30[i]) & is.na(train.data72$vl_mth24[i]) & is.na(train.data72$vl_mth18[i]) & !is.na(train.data72$vl_mth12[i]) & train.data72$vl_mth12[i] < 1000) |
       (is.na(train.data72$vl_mth36[i]) & is.na(train.data72$vl_mth30[i]) & is.na(train.data72$vl_mth24[i]) & is.na(train.data72$vl_mth18[i]) & is.na(train.data72$vl_mth12[i]) & !is.na(train.data72$vl_mth6[i]) & train.data72$vl_mth6[i] < 1000) |
       (is.na(train.data72$vl_mth36[i]) & is.na(train.data72$vl_mth30[i]) & is.na(train.data72$vl_mth24[i]) & is.na(train.data72$vl_mth18[i]) & is.na(train.data72$vl_mth12[i]) & is.na(train.data72$vl_mth6[i]) & !is.na(train.data72$basevl[i]) & train.data72$basevl[i] < 1000))
      & !is.na(train.data72$vl_mth42[i]) & train.data72$vl_mth42[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 48 months
  if (((!is.na(train.data72$vl_mth42[i]) & train.data72$vl_mth42[i] < 1000) |
       (is.na(train.data72$vl_mth42[i]) & !is.na(train.data72$vl_mth36[i]) & train.data72$vl_mth36[i] < 1000) |
       (is.na(train.data72$vl_mth42[i]) & is.na(train.data72$vl_mth36[i]) & !is.na(train.data72$vl_mth30[i]) & train.data72$vl_mth30[i] < 1000) |
       (is.na(train.data72$vl_mth42[i]) & is.na(train.data72$vl_mth36[i]) & is.na(train.data72$vl_mth30[i]) & !is.na(train.data72$vl_mth24[i]) & train.data72$vl_mth24[i] < 1000) |
       (is.na(train.data72$vl_mth42[i]) & is.na(train.data72$vl_mth36[i]) & is.na(train.data72$vl_mth30[i]) & is.na(train.data72$vl_mth24[i]) & !is.na(train.data72$vl_mth18[i]) & train.data72$vl_mth18[i] < 1000) |
       (is.na(train.data72$vl_mth42[i]) & is.na(train.data72$vl_mth36[i]) & is.na(train.data72$vl_mth30[i]) & is.na(train.data72$vl_mth24[i]) & is.na(train.data72$vl_mth18[i]) & !is.na(train.data72$vl_mth12[i]) & train.data72$vl_mth12[i] < 1000) |
       (is.na(train.data72$vl_mth42[i]) & is.na(train.data72$vl_mth36[i]) & is.na(train.data72$vl_mth30[i]) & is.na(train.data72$vl_mth24[i]) & is.na(train.data72$vl_mth18[i]) & is.na(train.data72$vl_mth12[i]) & !is.na(train.data72$vl_mth6[i]) & train.data72$vl_mth6[i] < 1000) | 
       (is.na(train.data72$vl_mth42[i]) & is.na(train.data72$vl_mth36[i]) & is.na(train.data72$vl_mth30[i]) & is.na(train.data72$vl_mth24[i]) & is.na(train.data72$vl_mth18[i]) & is.na(train.data72$vl_mth12[i]) & is.na(train.data72$vl_mth6[i]) & !is.na(train.data72$basevl[i]) & train.data72$basevl[i] < 1000))
      & !is.na(train.data72$vl_mth48[i]) & train.data72$vl_mth48[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 54 months
  if (((!is.na(train.data72$vl_mth48[i]) & train.data72$vl_mth48[i] < 1000) |
       (is.na(train.data72$vl_mth48[i]) & !is.na(train.data72$vl_mth42[i]) & train.data72$vl_mth42[i] < 1000) |
       (is.na(train.data72$vl_mth48[i]) & is.na(train.data72$vl_mth42[i]) & !is.na(train.data72$vl_mth36[i]) & train.data72$vl_mth36[i] < 1000) |
       (is.na(train.data72$vl_mth48[i]) & is.na(train.data72$vl_mth42[i]) & is.na(train.data72$vl_mth36[i]) & !is.na(train.data72$vl_mth30[i]) & train.data72$vl_mth30[i] < 1000) |
       (is.na(train.data72$vl_mth48[i]) & is.na(train.data72$vl_mth42[i]) & is.na(train.data72$vl_mth36[i]) & is.na(train.data72$vl_mth30[i]) & !is.na(train.data72$vl_mth24[i]) & train.data72$vl_mth24[i] < 1000) |
       (is.na(train.data72$vl_mth48[i]) & is.na(train.data72$vl_mth42[i]) & is.na(train.data72$vl_mth36[i]) & is.na(train.data72$vl_mth30[i]) & is.na(train.data72$vl_mth24[i]) & !is.na(train.data72$vl_mth18[i]) & train.data72$vl_mth18[i] < 1000) |
       (is.na(train.data72$vl_mth48[i]) & is.na(train.data72$vl_mth42[i]) & is.na(train.data72$vl_mth36[i]) & is.na(train.data72$vl_mth30[i]) & is.na(train.data72$vl_mth24[i]) & is.na(train.data72$vl_mth18[i]) & !is.na(train.data72$vl_mth12[i]) & train.data72$vl_mth12[i] < 1000) | 
       (is.na(train.data72$vl_mth48[i]) & is.na(train.data72$vl_mth42[i]) & is.na(train.data72$vl_mth36[i]) & is.na(train.data72$vl_mth30[i]) & is.na(train.data72$vl_mth24[i]) & is.na(train.data72$vl_mth18[i]) & is.na(train.data72$vl_mth12[i]) & !is.na(train.data72$vl_mth6[i]) & train.data72$vl_mth6[i] < 1000) |
       (is.na(train.data72$vl_mth48[i]) & is.na(train.data72$vl_mth42[i]) & is.na(train.data72$vl_mth36[i]) & is.na(train.data72$vl_mth30[i]) & is.na(train.data72$vl_mth24[i]) & is.na(train.data72$vl_mth18[i]) & is.na(train.data72$vl_mth12[i]) & is.na(train.data72$vl_mth6[i]) & !is.na(train.data72$basevl[i]) & train.data72$basevl[i] < 1000))
      & !is.na(train.data72$vl_mth54[i]) & train.data72$vl_mth54[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 60 months
  if (((!is.na(train.data72$vl_mth54[i]) & train.data72$vl_mth54[i] < 1000) |
       (is.na(train.data72$vl_mth54[i]) & !is.na(train.data72$vl_mth48[i]) & train.data72$vl_mth48[i] < 1000) |
       (is.na(train.data72$vl_mth54[i]) & is.na(train.data72$vl_mth48[i]) & !is.na(train.data72$vl_mth42[i]) & train.data72$vl_mth42[i] < 1000) |
       (is.na(train.data72$vl_mth54[i]) & is.na(train.data72$vl_mth48[i]) & is.na(train.data72$vl_mth42[i]) & !is.na(train.data72$vl_mth36[i]) & train.data72$vl_mth36[i] < 1000) |
       (is.na(train.data72$vl_mth54[i]) & is.na(train.data72$vl_mth48[i]) & is.na(train.data72$vl_mth42[i]) & is.na(train.data72$vl_mth36[i]) & !is.na(train.data72$vl_mth30[i]) & train.data72$vl_mth30[i] < 1000) |
       (is.na(train.data72$vl_mth54[i]) & is.na(train.data72$vl_mth48[i]) & is.na(train.data72$vl_mth42[i]) & is.na(train.data72$vl_mth36[i]) & is.na(train.data72$vl_mth30[i]) & !is.na(train.data72$vl_mth24[i]) & train.data72$vl_mth24[i] < 1000) |
       (is.na(train.data72$vl_mth54[i]) & is.na(train.data72$vl_mth48[i]) & is.na(train.data72$vl_mth42[i]) & is.na(train.data72$vl_mth36[i]) & is.na(train.data72$vl_mth30[i]) & is.na(train.data72$vl_mth24[i]) & !is.na(train.data72$vl_mth18[i]) & train.data72$vl_mth18[i] < 1000) | 
       (is.na(train.data72$vl_mth54[i]) & is.na(train.data72$vl_mth48[i]) & is.na(train.data72$vl_mth42[i]) & is.na(train.data72$vl_mth36[i]) & is.na(train.data72$vl_mth30[i]) & is.na(train.data72$vl_mth24[i]) & is.na(train.data72$vl_mth18[i]) & !is.na(train.data72$vl_mth12[i]) & train.data72$vl_mth12[i] < 1000) |
       (is.na(train.data72$vl_mth54[i]) & is.na(train.data72$vl_mth48[i]) & is.na(train.data72$vl_mth42[i]) & is.na(train.data72$vl_mth36[i]) & is.na(train.data72$vl_mth30[i]) & is.na(train.data72$vl_mth24[i]) & is.na(train.data72$vl_mth18[i]) & is.na(train.data72$vl_mth12[i]) & !is.na(train.data72$vl_mth6[i]) & train.data72$vl_mth6[i] < 1000) |
       (is.na(train.data72$vl_mth54[i]) & is.na(train.data72$vl_mth48[i]) & is.na(train.data72$vl_mth42[i]) & is.na(train.data72$vl_mth36[i]) & is.na(train.data72$vl_mth30[i]) & is.na(train.data72$vl_mth24[i]) & is.na(train.data72$vl_mth18[i]) & is.na(train.data72$vl_mth12[i]) & is.na(train.data72$vl_mth6[i]) & !is.na(train.data72$basevl[i]) & train.data72$basevl[i] < 1000))
      & !is.na(train.data72$vl_mth60[i]) & train.data72$vl_mth60[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  #checking for failures at 66 months
  if (((!is.na(train.data72$vl_mth60[i]) & train.data72$vl_mth60[i] < 1000) |
       (is.na(train.data72$vl_mth60[i]) & !is.na(train.data72$vl_mth54[i]) & train.data72$vl_mth54[i] < 1000) |
       (is.na(train.data72$vl_mth60[i]) & is.na(train.data72$vl_mth54[i]) & !is.na(train.data72$vl_mth48[i]) & train.data72$vl_mth48[i] < 1000) |
       (is.na(train.data72$vl_mth60[i]) & is.na(train.data72$vl_mth54[i]) & is.na(train.data72$vl_mth48[i]) & !is.na(train.data72$vl_mth42[i]) & train.data72$vl_mth42[i] < 1000) |
       (is.na(train.data72$vl_mth60[i]) & is.na(train.data72$vl_mth54[i]) & is.na(train.data72$vl_mth48[i]) & is.na(train.data72$vl_mth42[i]) & !is.na(train.data72$vl_mth36[i]) & train.data72$vl_mth36[i] < 1000) |
       (is.na(train.data72$vl_mth60[i]) & is.na(train.data72$vl_mth54[i]) & is.na(train.data72$vl_mth48[i]) & is.na(train.data72$vl_mth42[i]) & is.na(train.data72$vl_mth36[i]) & !is.na(train.data72$vl_mth30[i]) & train.data72$vl_mth30[i] < 1000) |
       (is.na(train.data72$vl_mth60[i]) & is.na(train.data72$vl_mth54[i]) & is.na(train.data72$vl_mth48[i]) & is.na(train.data72$vl_mth42[i]) & is.na(train.data72$vl_mth36[i]) & is.na(train.data72$vl_mth30[i]) & !is.na(train.data72$vl_mth24[i]) & train.data72$vl_mth24[i] < 1000) | 
       (is.na(train.data72$vl_mth60[i]) & is.na(train.data72$vl_mth54[i]) & is.na(train.data72$vl_mth48[i]) & is.na(train.data72$vl_mth42[i]) & is.na(train.data72$vl_mth36[i]) & is.na(train.data72$vl_mth30[i]) & is.na(train.data72$vl_mth24[i]) & !is.na(train.data72$vl_mth18[i]) & train.data72$vl_mth18[i] < 1000) |
       (is.na(train.data72$vl_mth60[i]) & is.na(train.data72$vl_mth54[i]) & is.na(train.data72$vl_mth48[i]) & is.na(train.data72$vl_mth42[i]) & is.na(train.data72$vl_mth36[i]) & is.na(train.data72$vl_mth30[i]) & is.na(train.data72$vl_mth24[i]) & is.na(train.data72$vl_mth18[i]) & !is.na(train.data72$vl_mth12[i]) & train.data72$vl_mth12[i] < 1000) |
       (is.na(train.data72$vl_mth60[i]) & is.na(train.data72$vl_mth54[i]) & is.na(train.data72$vl_mth48[i]) & is.na(train.data72$vl_mth42[i]) & is.na(train.data72$vl_mth36[i]) & is.na(train.data72$vl_mth30[i]) & is.na(train.data72$vl_mth24[i]) & is.na(train.data72$vl_mth18[i]) & is.na(train.data72$vl_mth12[i]) & !is.na(train.data72$vl_mth6[i]) & train.data72$vl_mth6[i] < 1000) |
       (is.na(train.data72$vl_mth60[i]) & is.na(train.data72$vl_mth54[i]) & is.na(train.data72$vl_mth48[i]) & is.na(train.data72$vl_mth42[i]) & is.na(train.data72$vl_mth36[i]) & is.na(train.data72$vl_mth30[i]) & is.na(train.data72$vl_mth24[i]) & is.na(train.data72$vl_mth18[i]) & is.na(train.data72$vl_mth12[i]) & is.na(train.data72$vl_mth6[i]) & !is.na(train.data72$basevl[i]) & train.data72$basevl[i] < 1000))
      & !is.na(train.data72$vl_mth66[i]) & train.data72$vl_mth66[i] >= 1000)
  {counter[i] = counter[i]+1}
  
  train.data72$prev_fail[i]=counter[i]
}

train.data72$visit=72

train.data72$vl <- train.data72$vl_mth72
train.data72$vllog <- train.data72$vl72

### getting the double failure status at month 72
for (i in 1:length(train.data72$study_id)) {
  if (!is.na(train.data72$prev_vl[i]) & train.data72$prev_vl_t[i]<72 & train.data72$prev_vl[i]>=1000 & train.data72$vl_mth72[i]>=1000)
  {train.data72$dub_fail[i] = 1}
  else if (!is.na(train.data72$prev_vl[i]) & train.data72$prev_vl_t[i]<72 & (train.data72$prev_vl[i]<1000 | train.data72$vl_mth72[i]<1000))
  {train.data72$dub_fail[i] = 0}
  else {train.data72$dub_fail[i] = NA}
}

check <- train.data72[is.na(train.data72$dub_fail),]
check2 <- train.data72[!is.na(train.data72$dub_fail),]
check3 <- train.data72[train.data72$dub_fail==1 & !is.na(train.data72$dub_fail),]


##### putting all of the individual timepoint dataset together

train.data <- rbind(train.data6, train.data12, train.data18, train.data24, train.data30, train.data36, train.data42, train.data48,
                    train.data54, train.data60, train.data66, train.data72)

#### creating variable lastVL which is 50 if VL below 50 and VL if VL > 1000

for (i in 1:length(train.data$study_id)){
  if (is.na(train.data$prev_vl[i])) {train.data$lastVL[i]=NA
                                  train.data$lastVLlog[i]=NA}
  else if (train.data$prev_vl[i]<=50){train.data$lastVL[i]=50
                                    train.data$lastVLlog[i]=log10(50)}
  else if (train.data$prev_vl[i]>50){train.data$lastVL[i]=train.data$prev_vl[i]
                                      train.data$lastVLlog[i]=log10(train.data$prev_vl[i])}
}


#### creating variable for time since last VL; 1 if 6 months or less, 0 else
train.data$lastVL_t <- ifelse(train.data$prev_vl_t<=6,1,0)

#### creating variable for time since enrollment; 1 if 12 months or less, 0 else
train.data$enroll_t <- ifelse(train.data$visit<=12,1,0)


#### creating interaction variable lastVLlog*lastVL_t
for (i in 1:length(train.data$study_id)){
  train.data$lastVL_time_int[i] <- train.data$lastVLlog[i]*train.data$lastVL_t[i]
}

### turning all categorical variables into binary variables for use with glmnet
train.data$male <- ifelse(train.data$sex=="M",1,0)
train.data$female <- ifelse(train.data$sex=="F",1,0)

train.data$basewho0 <- ifelse(train.data$basewho=="0",1,0)
train.data$basewho1 <- ifelse(train.data$basewho=="1",1,0)
train.data$basewho2 <- ifelse(train.data$basewho=="2",1,0)
train.data$basewho3 <- ifelse(train.data$basewho=="3",1,0)
train.data$basewho4 <- ifelse(train.data$basewho=="4",1,0)

train.data$buyamba <- ifelse(train.data$hub_name=="BUYAMBA",1,0)
train.data$kabira <- ifelse(train.data$hub_name=="KABIRA",1,0)
train.data$kakuuto <- ifelse(train.data$hub_name=="KAKUUTO",1,0)
train.data$kaleere <- ifelse(train.data$hub_name=="KALEERE",1,0)
train.data$kalisizo <- ifelse(train.data$hub_name=="KALISIZO",1,0)
train.data$kasaali <- ifelse(train.data$hub_name=="KASAALI",1,0)
train.data$kasasa <- ifelse(train.data$hub_name=="KASASA",1,0)
train.data$kayanja <- ifelse(train.data$hub_name=="KAYANJA",1,0)
train.data$kibaale <- ifelse(train.data$hub_name=="KIBAALE",1,0)
train.data$kifamba <- ifelse(train.data$hub_name=="KIFAMBA",1,0)
train.data$kyabigondo <- ifelse(train.data$hub_name=="KAYBIGONDO",1,0)
train.data$kyebe <- ifelse(train.data$hub_name=="KYEBE",1,0)
train.data$lwamaggwa <- ifelse(train.data$hub_name=="LWAMAGGWA",1,0)
train.data$lwanda <- ifelse(train.data$hub_name=="LWANDA",1,0)
train.data$lyantonde <- ifelse(train.data$hub_name=="LYANTONDE",1,0)
train.data$nabigasa <- ifelse(train.data$hub_name=="NABIGASA",1,0)
train.data$nakatoogo <- ifelse(train.data$hub_name=="NAKATOOGO",1,0)

train.data$reg_ABC3TCEFV <- ifelse(train.data$basereg=="ABC/3TC/EFV",1,0)
train.data$reg_ABC3TCLPVr <- ifelse(train.data$basereg=="ABC/3TC/LPVr",1,0)
train.data$reg_ABC3TCNVP <- ifelse(train.data$basereg=="ABC/3TC/NVP",1,0)
train.data$reg_CBVABC <- ifelse(train.data$basereg=="CBV/ABC",1,0)
train.data$reg_CBVEFV <- ifelse(train.data$basereg=="CBV/EFV",1,0)
train.data$reg_CBVNVP <- ifelse(train.data$basereg=="CBV/NVP",1,0)
train.data$reg_D4T3TCABC <- ifelse(train.data$basereg=="D4T/3TC/ABC",1,0)
train.data$reg_D4T3TCEFV <- ifelse(train.data$basereg=="D4T/3TC/EFV",1,0)
train.data$reg_D4T3TCLPVr <- ifelse(train.data$basereg=="D4T/3TC/LPVr",1,0)
train.data$reg_D4T3TCNVP <- ifelse(train.data$basereg=="D4T/3TC/NVP",1,0)
train.data$reg_TDF3TCEFV <- ifelse(train.data$basereg=="TDF/3TC/EFV",1,0)
train.data$reg_TDF3TCNVP <- ifelse(train.data$basereg=="TDF/3TC/NVP",1,0)
train.data$reg_TVDEFV <- ifelse(train.data$basereg=="TVD/EFV",1,0)

train.data$prev_fail0 <- ifelse(train.data$prev_fail==0,1,0)
train.data$prev_fail1 <- ifelse(train.data$prev_fail==1,1,0)
train.data$prev_fail2 <- ifelse(train.data$prev_fail==2,1,0)
train.data$prev_fail3 <- ifelse(train.data$prev_fail==3,1,0)


#### creating a variable for whether or not the current VL was a failure
train.data$fail_draw <- ifelse(train.data$vl>=1000,1,0)


#### adding a variable for high, low and medium risk groups based on VL we are trying to prdict
#### low risk is VL < 0 (lower limit of detection), medium risk is 0 <= VL <= 10000, high risk is VL > 10000

for (i in 1:length(train.data$study_id)){
  if (train.data$vl[i] < 0){train.data$label[i]="low"}
  else if (train.data$vl[i] >= 0 & train.data$vl[i] <= 10000){train.data$label[i]="medium"}
  else if (train.data$vl[i] > 10000){train.data$label[i]="high"}
}

#### creating a category variable (low, medium, high) for the last VL

for (i in 1:length(train.data$study_id)){
  if (train.data$lastVL[i] <= 50){train.data$lastVLcat[i]="low"}
  else if (train.data$lastVL[i] > 50 & train.data$lastVL[i] <= 10000){train.data$lastVLcat[i]="medium"}
  else if (train.data$lastVL[i] > 10000){train.data$lastVLcat[i]="high"}
}

# subsetting for only the variable we need
train.data2 <- select(train.data, study_id, sex, male, female, ageyrs,  age35, age50, age65, 
                      hub_name, buyamba, kabira, kakuuto, kaleere, kalisizo, kasaali, kasasa, kayanja,
                      kibaale, kifamba, kyabigondo, kyebe, lwamaggwa, lwanda, lyantonde, nabigasa, nakatoogo,
                      basewho, basewho0, basewho1, basewho2, basewho3, basewho4, 
                      basereg, reg_TVDEFV, reg_TDF3TCNVP, reg_TDF3TCEFV, reg_D4T3TCNVP, 
                      reg_D4T3TCLPVr, reg_D4T3TCEFV, reg_D4T3TCABC, reg_CBVNVP, reg_CBVEFV, 
                      reg_CBVABC, reg_ABC3TCNVP, reg_ABC3TCLPVr, reg_ABC3TCEFV,
                      prev_fail, prev_fail0, prev_fail1, prev_fail2, prev_fail3,
                      basecd4, bascd4log, basevl, basvllog, lastVL, lastVLlog,
                      lastVL_t, enroll_t, vl, vllog, fail_draw, dub_fail, lastVLcat, label)

# set working directory to location you want to write the dataset and uncomment below line
#setwd("")
#write.table(train.data2, paste0("train_data", format(Sys.time(),"%Y-%m-%d"), ".R"))





