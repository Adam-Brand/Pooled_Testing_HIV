#==============================================================================
# FILENAME: raw_data_clean.R
# PROJECT: 	Pooled testing in HIV
# PURPOSE: clean and add derived variables to the raw Uganda data 
# AUTHOR: Adam Brand

# CREATED:	20200224
# UPDATED: 	

# INPUT DATA: raw data form Uganda HIV clinic 			
# OUTPUT: output dataset Clean_data2020-03-06	

# R VERSION: 3.6.1
#==============================================================================
#Notes: 





# =============================================================================


#### adding relevant libraries
library(lubridate)
library(ggplot2)
library(dplyr)

#### set working directory to location with data
setwd("C:/Users/Barny/Dropbox/KI_Project_4/Data/data_uganda")

### reading in raw data csv file
raw <- read.csv("Raw_data_version1_2396.csv", header=T)

#### checking the names of all variables in the raw data
names(raw)

##### getting data type of each variable
sapply(raw, class)

### checking for duplicate study_id numbers
length(unique(raw$study_id))
length(unique(raw$study_id)) == nrow(raw)


### writing a function to reutrn number of NAs if numeric and number of missing if factor
miss <- function(x){
  if (is.numeric(x)) {
    sum(is.na(x))}
     else if (is.factor(x)){
    sum(x=="")}
    else if (is.character(x)){
      sum(x=="")}
   else {"NA"}
}

### setting mydata to raw
mydata <- raw


#### turning the date variables into numeric date variables using package lubridate
mydata$cd4basedate <- dmy(raw$cd4basedate)
mydata$artstartdate <- dmy(raw$artstartdate)
mydata$datedied <- dmy(raw$datedied)
mydata$seclinedate <- dmy(raw$seclinedate)
mydata$visitdate6 <- dmy(raw$visitdate6)
mydata$visitdate12 <- dmy(raw$visitdate12)
mydata$visitdate18 <- dmy(raw$visitdate18)
mydata$visitdate24 <- dmy(raw$visitdate24)
mydata$visitdate30 <- dmy(raw$visitdate30)
mydata$visitdate36 <- dmy(raw$visitdate36)
mydata$visitdate42 <- dmy(raw$visitdate42)
mydata$visitdate48 <- dmy(raw$visitdate48)
mydata$visitdate54 <- dmy(raw$visitdate54)
mydata$visitdate60 <- dmy(raw$visitdate60)
mydata$visitdate66 <- dmy(raw$visitdate66)
mydata$visitdate72 <- dmy(raw$visitdate72)

#### checking that the date conversion didn't remove any dates
summary(mydata$cd4basedate)
summary(mydata$artstartdate)
summary(mydata$datedied)
summary(mydata$seclinedate)
summary(mydata$visitdate6)
summary(mydata$visitdate12)
summary(mydata$visitdate18)
summary(mydata$visitdate24)
summary(mydata$visitdate30)
summary(mydata$visitdate36)
summary(mydata$visitdate42)
summary(mydata$visitdate48)
summary(mydata$visitdate54)
summary(mydata$visitdate60)
summary(mydata$visitdate66)
summary(mydata$visitdate72)


#### checking number of missing variables for numeric variables
sapply(mydata, function(x) miss(x))

##### determining the types of entries for numeric, categorical variables
levels(as.factor(mydata$secondline))
levels(as.factor(mydata$ltfu))
levels(as.factor(mydata$basewho))


#### looking at a summary of baseline VL values
summary(mydata$basevl)
length(which(mydata$basevl==-1)) # 80 subjects has undetectable VL at baseline
length(which(mydata$basevl< 1000 )) # 121 subjects have VL below 1000 at baseline

### creating variables for age categories
mydata$age16 <- ifelse(mydata$ageyrs<=16,1,0)
mydata$age35 <- ifelse(mydata$ageyrs>16 & mydata$ageyrs<=35,1,0)
mydata$age50 <- ifelse(mydata$ageyrs>35 & mydata$ageyrs<=50,1,0)
mydata$age65 <- ifelse(mydata$ageyrs>50,1,0)

#creating a variable for base CDD4
for (i in 1:length(mydata$basecd4)) {
  if (is.na(mydata$basecd4[i])) {mydata$bascd4log[i]=NA} 
  else if (mydata$basecd4[i]<=0) {mydata$bascd4log[i]=log10(1)}
  else {mydata$bascd4log[i]=log10(mydata$basecd4[i])}
}

## changing basewho into a factor variable
mydata$basewho <- as.factor(mydata$basewho)

### creating variables for log10 viral load

#setting lower limit of detection
lowdet <- 50

## if HIV VL is < 0, it indicates below LLD, so setting the VL to the LLD
for (i in 1:length(mydata$basevl)) {
  if (is.na(mydata$basevl[i])) {mydata$basvllog[i]=NA} 
  else if (mydata$basevl[i]<=0) {mydata$basvllog[i]=log10(lowdet)}
  else {mydata$basvllog[i]=log10(mydata$basevl[i])}
}

for (i in 1:length(mydata$vl_mth6)) {
  if (is.na(mydata$vl_mth6[i])) {mydata$vl6[i]=NA} 
  else if (mydata$vl_mth6[i]<=0) {mydata$vl6[i]=log10(lowdet)}
  else {mydata$vl6[i]=log10(mydata$vl_mth6[i])}
}

for (i in 1:length(mydata$vl_mth12)) {
  if (is.na(mydata$vl_mth12[i])) {mydata$vl12[i]=NA} 
  else if (mydata$vl_mth12[i]<=0) {mydata$vl12[i]=log10(lowdet)}
  else {mydata$vl12[i]=log10(mydata$vl_mth12[i])}
}

for (i in 1:length(mydata$vl_mth18)) {
  if (is.na(mydata$vl_mth18[i])) {mydata$vl18[i]=NA} 
  else if (mydata$vl_mth18[i]<=0) {mydata$vl18[i]=log10(lowdet)}
  else {mydata$vl18[i]=log10(mydata$vl_mth18[i])}
}

for (i in 1:length(mydata$vl_mth24)) {
  if (is.na(mydata$vl_mth24[i])) {mydata$vl24[i]=NA} 
  else if (mydata$vl_mth24[i]<=0) {mydata$vl24[i]=log10(lowdet)}
  else {mydata$vl24[i]=log10(mydata$vl_mth24[i])}
}

for (i in 1:length(mydata$vl_mth30)) {
  if (is.na(mydata$vl_mth30[i])) {mydata$vl30[i]=NA} 
  else if (mydata$vl_mth30[i]<=0) {mydata$vl30[i]=log10(lowdet)}
  else {mydata$vl30[i]=log10(mydata$vl_mth30[i])}
}

for (i in 1:length(mydata$vl_mth36)) {
  if (is.na(mydata$vl_mth36[i])) {mydata$vl36[i]=NA} 
  else if (mydata$vl_mth36[i]<=0) {mydata$vl36[i]=log10(lowdet)}
  else {mydata$vl36[i]=log10(mydata$vl_mth36[i])}
}

for (i in 1:length(mydata$vl_mth42)) {
  if (is.na(mydata$vl_mth42[i])) {mydata$vl42[i]=NA} 
  else if (mydata$vl_mth42[i]<=0) {mydata$vl42[i]=log10(lowdet)}
  else {mydata$vl42[i]=log10(mydata$vl_mth42[i])}
}

for (i in 1:length(mydata$vl_mth48)) {
  if (is.na(mydata$vl_mth48[i])) {mydata$vl48[i]=NA} 
  else if (mydata$vl_mth48[i]<=0) {mydata$vl48[i]=log10(lowdet)}
  else {mydata$vl48[i]=log10(mydata$vl_mth48[i])}
}

for (i in 1:length(mydata$vl_mth54)) {
  if (is.na(mydata$vl_mth54[i])) {mydata$vl54[i]=NA} 
  else if (mydata$vl_mth54[i]<=0) {mydata$vl54[i]=log10(lowdet)}
  else {mydata$vl54[i]=log10(mydata$vl_mth54[i])}
}

for (i in 1:length(mydata$vl_mth60)) {
  if (is.na(mydata$vl_mth60[i])) {mydata$vl60[i]=NA} 
  else if (mydata$vl_mth60[i]<=0) {mydata$vl60[i]=log10(lowdet)}
  else {mydata$vl60[i]=log10(mydata$vl_mth60[i])}
}

for (i in 1:length(mydata$vl_mth66)) {
  if (is.na(mydata$vl_mth66[i])) {mydata$vl66[i]=NA} 
  else if (mydata$vl_mth66[i]<=0) {mydata$vl66[i]=log10(lowdet)}
  else {mydata$vl66[i]=log10(mydata$vl_mth66[i])}
}

for (i in 1:length(mydata$vl_mth72)) {
  if (is.na(mydata$vl_mth72[i])) {mydata$vl72[i]=NA} 
  else if (mydata$vl_mth72[i]<=0) {mydata$vl72[i]=log10(lowdet)}
  else {mydata$vl72[i]=log10(mydata$vl_mth72[i])}
}



## creating a binary variable, base_supp, which is 1 if VL suppressed at baseline and 0 else
for (i in 1:nrow(mydata)){
  if (mydata$basevl[i] < 1000 & !is.na(mydata$basevl[i])) {mydata$base_supp[i]=1} 
  else if (!is.na(mydata$basevl[i])) {mydata$base_supp[i]=0}
  else {mydata$base_supp[i]=NA}
}

# checking the code
sum(mydata$basevl, na.rm=T)
summary(mydata$base_supp) # 6.22% suppressed at baseline
miss(mydata$base_supp)

### creating the rest of the binary suppression variables for each time point
for (i in 1:nrow(mydata)){
  if (mydata$vl_mth6[i] < 1000 & !is.na(mydata$vl_mth6[i])) {mydata$mth6_supp[i]=1} 
  else if (!is.na(mydata$vl_mth6[i])) {mydata$mth6_supp[i]=0}
  else {mydata$mth6_supp[i]=NA}
}
sum(mydata$vl_mth6, na.rm=T)
summary(mydata$mth6_supp)
summary(mydata$vl_mth6)

for (i in 1:nrow(mydata)){
  if (mydata$vl_mth12[i] < 1000 & !is.na(mydata$vl_mth12[i])) {mydata$mth12_supp[i]=1} 
  else if (!is.na(mydata$vl_mth12[i])) {mydata$mth12_supp[i]=0}
  else {mydata$mth12_supp[i]=NA}
}
sum(mydata$vl_mth12, na.rm=T)
summary(mydata$mth12_supp)
summary(mydata$vl_mth12)

for (i in 1:nrow(mydata)){
  if (mydata$vl_mth18[i] < 1000 & !is.na(mydata$vl_mth18[i])) {mydata$mth18_supp[i]=1} 
  else if (!is.na(mydata$vl_mth18[i])) {mydata$mth18_supp[i]=0}
  else {mydata$mth18_supp[i]=NA}
}
sum(mydata$vl_mth18, na.rm=T)
summary(mydata$mth18_supp)
summary(mydata$vl_mth18)


for (i in 1:nrow(mydata)){
  if (mydata$vl_mth24[i] < 1000 & !is.na(mydata$vl_mth24[i])) {mydata$mth24_supp[i]=1} 
  else if (!is.na(mydata$vl_mth24[i])) {mydata$mth24_supp[i]=0}
  else {mydata$mth24_supp[i]=NA}
}
sum(mydata$vl_mth24, na.rm=T)
summary(mydata$mth24_supp)
summary(mydata$vl_mth24)


for (i in 1:nrow(mydata)){
  if (mydata$vl_mth30[i] < 1000 & !is.na(mydata$vl_mth30[i])) {mydata$mth30_supp[i]=1} 
  else if (!is.na(mydata$vl_mth30[i])) {mydata$mth30_supp[i]=0}
  else {mydata$mth30_supp[i]=NA}
}
sum(mydata$vl_mth30, na.rm=T)
summary(mydata$mth30_supp)
summary(mydata$vl_mth30)


for (i in 1:nrow(mydata)){
  if (mydata$vl_mth36[i] < 1000 & !is.na(mydata$vl_mth36[i])) {mydata$mth36_supp[i]=1} 
  else if (!is.na(mydata$vl_mth36[i])) {mydata$mth36_supp[i]=0}
  else {mydata$mth36_supp[i]=NA}
}
sum(mydata$vl_mth36, na.rm=T)
summary(mydata$mth36_supp)
summary(mydata$vl_mth36)


for (i in 1:nrow(mydata)){
  if (mydata$vl_mth42[i] < 1000 & !is.na(mydata$vl_mth42[i])) {mydata$mth42_supp[i]=1} 
  else if (!is.na(mydata$vl_mth42[i])) {mydata$mth42_supp[i]=0}
  else {mydata$mth42_supp[i]=NA}
}
sum(mydata$vl_mth42, na.rm=T)
summary(mydata$mth42_supp)
summary(mydata$vl_mth42)


for (i in 1:nrow(mydata)){
  if (mydata$vl_mth48[i] < 1000 & !is.na(mydata$vl_mth48[i])) {mydata$mth48_supp[i]=1} 
  else if (!is.na(mydata$vl_mth48[i])) {mydata$mth48_supp[i]=0}
  else {mydata$mth48_supp[i]=NA}
}
sum(mydata$vl_mth48, na.rm=T)
summary(mydata$mth48_supp)
summary(mydata$vl_mth48)


for (i in 1:nrow(mydata)){
  if (mydata$vl_mth54[i] < 1000 & !is.na(mydata$vl_mth54[i])) {mydata$mth54_supp[i]=1} 
  else if (!is.na(mydata$vl_mth54[i])) {mydata$mth54_supp[i]=0}
  else {mydata$mth54_supp[i]=NA}
}
sum(mydata$vl_mth54, na.rm=T)
summary(mydata$mth54_supp)
summary(mydata$vl_mth54)

for (i in 1:nrow(mydata)){
  if (mydata$vl_mth60[i] < 1000 & !is.na(mydata$vl_mth60[i])) {mydata$mth60_supp[i]=1} 
  else if (!is.na(mydata$vl_mth60[i])) {mydata$mth60_supp[i]=0}
  else {mydata$mth60_supp[i]=NA}
}
sum(mydata$vl_mth60, na.rm=T)
summary(mydata$mth60_supp)
summary(mydata$vl_mth60)


for (i in 1:nrow(mydata)){
  if (mydata$vl_mth66[i] < 1000 & !is.na(mydata$vl_mth66[i])) {mydata$mth66_supp[i]=1} 
  else if (!is.na(mydata$vl_mth66[i])) {mydata$mth66_supp[i]=0}
  else {mydata$mth66_supp[i]=NA}
}
sum(mydata$vl_mth66, na.rm=T)
summary(mydata$mth66_supp)
summary(mydata$vl_mth66)


for (i in 1:nrow(mydata)){
  if (mydata$vl_mth72[i] < 1000 & !is.na(mydata$vl_mth72[i])) {mydata$mth72_supp[i]=1} 
  else if (!is.na(mydata$vl_mth72[i])) {mydata$mth72_supp[i]=0}
  else {mydata$mth72_supp[i]=NA}
}
sum(mydata$vl_mth72, na.rm=T)
summary(mydata$mth72_supp)
summary(mydata$vl_mth72)


##### deriving variable earlysupp which is the earliest date the subject had their VL suppressed
for (i in 1:nrow(mydata)){
    if (!is.na(mydata$base_supp[i]) & mydata$base_supp[i]==1 & !is.na(mydata$cd4basedate[i])) {mydata$earlysupp[i]=mydata$cd4basedate[i]}
    else if (!is.na(mydata$mth6_supp[i]) & mydata$mth6_supp[i]==1 & !is.na(mydata$visitdate6[i])) {mydata$earlysupp[i]=mydata$visitdate6[i]}
    else if (!is.na(mydata$mth12_supp[i]) & mydata$mth12_supp[i]==1 & !is.na(mydata$visitdate12[i])) {mydata$earlysupp[i]=mydata$visitdate12[i]}
    else if (!is.na(mydata$mth18_supp[i]) & mydata$mth18_supp[i]==1 & !is.na(mydata$visitdate18[i])) {mydata$earlysupp[i]=mydata$visitdate18[i]}
    else if (!is.na(mydata$mth24_supp[i]) & mydata$mth24_supp[i]==1 & !is.na(mydata$visitdate24[i])) {mydata$earlysupp[i]=mydata$visitdate24[i]}
    else if (!is.na(mydata$mth30_supp[i]) & mydata$mth30_supp[i]==1 & !is.na(mydata$visitdate30[i])) {mydata$earlysupp[i]=mydata$visitdate30[i]}
    else if (!is.na(mydata$mth36_supp[i]) & mydata$mth36_supp[i]==1 & !is.na(mydata$visitdate36[i])) {mydata$earlysupp[i]=mydata$visitdate36[i]}
    else if (!is.na(mydata$mth42_supp[i]) & mydata$mth42_supp[i]==1 & !is.na(mydata$visitdate42[i])) {mydata$earlysupp[i]=mydata$visitdate42[i]}
    else if (!is.na(mydata$mth48_supp[i]) & mydata$mth48_supp[i]==1 & !is.na(mydata$visitdate48[i])) {mydata$earlysupp[i]=mydata$visitdate48[i]}
    else if (!is.na(mydata$mth54_supp[i]) & mydata$mth54_supp[i]==1 & !is.na(mydata$visitdate54[i])) {mydata$earlysupp[i]=mydata$visitdate54[i]}
    else if (!is.na(mydata$mth60_supp[i]) & mydata$mth60_supp[i]==1 & !is.na(mydata$visitdate60[i])) {mydata$earlysupp[i]=mydata$visitdate60[i]}
    else if (!is.na(mydata$mth66_supp[i]) & mydata$mth66_supp[i]==1 & !is.na(mydata$visitdate66[i])) {mydata$earlysupp[i]=mydata$visitdate66[i]}
    else if (!is.na(mydata$mth72_supp[i]) & mydata$mth72_supp[i]==1 & !is.na(mydata$visitdate72[i])) {mydata$earlysupp[i]=mydata$visitdate72[i]}
    else {NA}
}

##### converting the format of earlysupp to Date format
mydata$earlysupp <- as.Date(mydata$earlysupp, origin)


#### creating the previous failure variable which is the earliest date after their first suppresion
#### that they failed treatment, i.e., VL >= 1000
mydata$fail <- as.numeric(NA)
mydata$earlyfaildt <- as.Date(NA)

for (i in 1:length(mydata$study_id)){
  if (!is.na(mydata$earlysupp[i])) {
    if (!is.na(mydata$visitdate6[i]) & mydata$visitdate6[i] > mydata$earlysupp[i] & !is.na(mydata$vl_mth6[i]) & mydata$vl_mth6[i] >= 1000) 
    {mydata$fail[i]=1  
    mydata$earlyfaildt[i]=mydata$visitdate6[i]}
    else if (!is.na(mydata$visitdate12[i]) & mydata$visitdate12[i] > mydata$earlysupp[i] & !is.na(mydata$vl_mth12[i]) & mydata$vl_mth12[i] >= 1000) 
    {mydata$fail[i]=1 
    mydata$earlyfaildt[i]=mydata$visitdate12[i]}
    else if (!is.na(mydata$visitdate18[i]) & mydata$visitdate18[i] > mydata$earlysupp[i] & !is.na(mydata$vl_mth18[i]) & mydata$vl_mth18[i] >= 1000) 
    {mydata$fail[i]=1 
    mydata$earlyfaildt[i]=mydata$visitdate18[i]}
    else if (!is.na(mydata$visitdate24[i]) & mydata$visitdate24[i] > mydata$earlysupp[i] & !is.na(mydata$vl_mth24[i]) & mydata$vl_mth24[i] >= 1000) 
    {mydata$fail[i]=1 
    mydata$earlyfaildt[i]=mydata$visitdate24[i]}
    else if (!is.na(mydata$visitdate30[i]) & mydata$visitdate30[i] > mydata$earlysupp[i] & !is.na(mydata$vl_mth30[i]) & mydata$vl_mth30[i] >= 1000) 
    {mydata$fail[i]=1
    mydata$earlyfaildt[i]=mydata$visitdate30[i]}
    else if (!is.na(mydata$visitdate36[i]) & mydata$visitdate36[i] > mydata$earlysupp[i] & !is.na(mydata$vl_mth36[i]) & mydata$vl_mth36[i] >= 1000) 
    {mydata$fail[i]=1
    mydata$earlyfaildt[i]=mydata$visitdate36[i]}
    else if (!is.na(mydata$visitdate42[i]) & mydata$visitdate42[i] > mydata$earlysupp[i] & !is.na(mydata$vl_mth42[i]) & mydata$vl_mth42[i] >= 1000) 
    {mydata$fail[i]=1
    mydata$earlyfaildt[i]=mydata$visitdate42[i]}
    else if (!is.na(mydata$visitdate48[i]) & mydata$visitdate48[i] > mydata$earlysupp[i] & !is.na(mydata$vl_mth48[i]) & mydata$vl_mth48[i] >= 1000) 
    {mydata$fail[i]=1 
    mydata$earlyfaildt[i]=mydata$visitdate48[i]}
    else if (!is.na(mydata$visitdate54[i]) & mydata$visitdate54[i] > mydata$earlysupp[i] & !is.na(mydata$vl_mth54[i]) & mydata$vl_mth54[i] >= 1000) 
    {mydata$fail[i]=1 
    mydata$earlyfaildt[i]=mydata$visitdate54[i]}
    else if (!is.na(mydata$visitdate60[i]) & mydata$visitdate60[i] > mydata$earlysupp[i] & !is.na(mydata$vl_mth60[i]) & mydata$vl_mth60[i] >= 1000) 
    {mydata$fail[i]=1 
    mydata$earlyfaildt[i]=mydata$visitdate60[i]}
    else if (!is.na(mydata$visitdate66[i]) & mydata$visitdate66[i] > mydata$earlysupp[i] & !is.na(mydata$vl_mth66[i]) & mydata$vl_mth66[i] >= 1000) 
    {mydata$fail[i]=1 
    mydata$earlyfaildt[i]=mydata$visitdate66[i]}
    else if (!is.na(mydata$visitdate72[i]) & mydata$visitdate72[i] > mydata$earlysupp[i] & !is.na(mydata$vl_mth72[i]) & mydata$vl_mth72[i] >= 1000) 
    {mydata$fail[i]=1 
    mydata$earlyfaildt[i]=mydata$visitdate72[i]}
  }
  
  else if (mydata$secondline[i]==1) {mydata$fail[i]=1
  mydata$earlyfaildt[i]=mydata$seclinedate[i]}
  else if (is.na(mydata$fail[i])) {mydata$fail[i]=0}
  
}

for (i in 1:length(mydata$study_id)){
  if (is.na(mydata$fail[i])){mydata$fail[i]=0}
}



##### End of cleaning and derivation code; this writes the dataset to a csv file

##set working directory to location to write clean_data, and uncomment below line
setwd("")
#write.csv2(mydata, paste0("Clean_data", format(Sys.time(),"%Y-%m-%d"), ".csv"))



##### creating frequency distributions of VL levels ommitting baseline
##### this code is intended to visualize the distribution of VLs post-baseline
##### none of these are variables in the dataset, the visualization is used to
##### get an idea of how to simulate VLs like those seen in practice

## putting all post-baseline VL measures into one dataset

VL6 <- data.frame(cbind(mydata$study_id, mydata$vl_mth6))
colnames(VL6) <- c("ID", "VL")

VL12 <- data.frame(cbind(mydata$study_id, mydata$vl_mth12))
colnames(VL12) <- c("ID", "VL")

VL18 <- data.frame(cbind(mydata$study_id, mydata$vl_mth18))
colnames(VL18) <- c("ID", "VL")

VL24 <- data.frame(cbind(mydata$study_id, mydata$vl_mth24))
colnames(VL24) <- c("ID", "VL")

VL30 <- data.frame(cbind(mydata$study_id, mydata$vl_mth30))
colnames(VL30) <- c("ID", "VL")

VL36 <- data.frame(cbind(mydata$study_id, mydata$vl_mth36))
colnames(VL36) <- c("ID", "VL")

VL42 <- data.frame(cbind(mydata$study_id, mydata$vl_mth42))
colnames(VL42) <- c("ID", "VL")

VL48 <- data.frame(cbind(mydata$study_id, mydata$vl_mth48))
colnames(VL48) <- c("ID", "VL")

VL54 <- data.frame(cbind(mydata$study_id, mydata$vl_mth54))
colnames(VL54) <- c("ID", "VL")

VL60 <- data.frame(cbind(mydata$study_id, mydata$vl_mth60))
colnames(VL60) <- c("ID", "VL")

VL66 <- data.frame(cbind(mydata$study_id, mydata$vl_mth66))
colnames(VL66) <- c("ID", "VL")

VL72 <- data.frame(cbind(mydata$study_id, mydata$vl_mth72))
colnames(VL72) <- c("ID", "VL")

VLdata <- data.frame(rbind(VL6,VL12,VL18,VL24,VL30,VL36,VL42,VL48,VL54,VL60,VL66,VL72))
# removing NAs
VLdata <- VLdata[!is.na(VLdata$VL),]

# setting all subjects with undetectable VL to the lower limit of detection set at 50
for (i in 1:length(VLdata$VL)){
    if (VLdata$VL[i]==-1) {VLdata$VL[i]=50}
}

#Creating a variable that is the log10 of VL
VLdata$VL10 <- log10(VLdata$VL)

#proportion of treatment failures in population at 1000 - 0.05664436 treatment failure rate
sum(VLdata$VL >= 1000)/length(VLdata$VL)
sum(VLdata$VL >= 500)/length(VLdata$VL)

summary(VLdata$VL10)

#making a histogram of post-baseline VL values on log10 scale
hist(VLdata$VL10,
     main="Distribution of Uganda VLs",
     xlab="Log 10 Viral Load",
     breaks=25,
     prob=TRUE)

lines(density(VLdata$VL10))

hist(log10(check),
     main="Distribution of simulated VLs",
     xlab="Log 10 Viral Load",
     breaks=25,
     prob=TRUE)

lines(density(log10(check)))






