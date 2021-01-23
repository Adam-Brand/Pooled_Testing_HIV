#==============================================================================
# FILENAME: data_permute.R
# PROJECT: 	Pooled testing in HIV
# PURPOSE: permutes the outcome data from the inout dataset to
#          eliminate the association between the predictors and the outcome
# AUTHOR: Adam Brand

# INPUT: 2 datasets
#        SD0_data_500_run2.rds
#        SD1.0_data_500_run2.rds

# OUTPUT: 2 datasets
#         - SD0_data_500_run2_perm.rds
#         - SD1.0_data_500_run2_perm.rds


# R VERSION: 3.6.1
#==============================================================================
#Notes: 





# =============================================================================

##### Data generation program AB 2020-03-03
##### Testing
library(shiny)
library(dplyr)
library(here)

SD0_data <- readRDS("SimData/Records_used_in_500_matrices_hypred/SD0_data_500_run2.rds")
SD1_data <- readRDS("SimData/Records_used_in_500_matrices_hypred/SD1.0_data_500_run2.rds")

set.seed(191919191)
VL0_perm <- sample(SD0_data$VL, replace=FALSE)
SD0_data_perm <- SD0_data
SD0_data_perm$VL <- VL0_perm


set.seed(191919191)
VL1_perm <- sample(SD1_data$VL, replace=FALSE)
SD1_data_perm <- SD1_data
SD1_data_perm$VL <- VL1_perm

saveRDS(SD0_data_perm, file=here("SimData","Records_used_in_500_matrices_hypred","SD0_data_500_run2_perm.rds"))
saveRDS(SD1_data_perm, file=here("SimData","Records_used_in_500_matrices_hypred","SD1_data_500_run2_perm.rds"))


