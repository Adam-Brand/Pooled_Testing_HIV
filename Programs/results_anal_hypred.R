#==============================================================================
# FILENAME: results_anal_hypred.R
# PROJECT: 	Pooled testing in HIV
# PURPOSE: compile the method evaluation results for HyPred from the simulated data into performance summaries
#          and prints code for latex tables
# AUTHOR: Adam Brand


# INPUT datasets: the results datasets in the 'Results' folder


# R VERSION: 3.6.1
#==============================================================================
#Notes: 





# =============================================================================

library(qwraps2)
library(dplyr)
library(rlang)
library(xtable)
library(PropCIs)
library(DescTools)

# set working directory to location of results
setwd("")

hypred_AGAIG_ME.05 <- read.table("Results_Hypred_AGAIG_ME.05_rand.R")
hypred_AGAIG_ME.12 <- read.table("Results_Hypred_AGAIG_ME.12_rand.R")
hypred_AGAIG_ME.25 <- read.table("Results_Hypred_AGAIG_ME.25_rand.R")
hypred_AGAIG_ME.5 <- read.table("Results_Hypred_AGAIG_ME.5_rand.R")

hypred_reverse_ME.05 <- read.table("Results_Hypred_reverse_ME.05_rand.R")
hypred_reverse_ME.12 <- read.table("Results_Hypred_reverse_ME.12_rand.R")
hypred_reverse_ME.25 <- read.table("Results_Hypred_reverse_ME.25_rand.R")
hypred_reverse_ME.5 <- read.table("Results_Hypred_reverse_ME.5_rand.R")

                ##############################

# set working directory to location of results
setwd("")

hypred_AGAIG_SD0_ME0 <- read.table("Results_Hypred_AGAIG_ME0_rand.R")
hypred_AGAIG_SD0_ME.25 <- read.table("Results_Hypred_AGAIG_ME.25_rand.R")
hypred_AGAIG_SD0_ME.5 <- read.table("Results_Hypred_AGAIG_ME.5_rand.R")
hypred_AGAIG_SD0_ME.75 <- read.table("Results_Hypred_AGAIG_ME.75_rand.R")

hypred_reverse_SD0_ME0 <- read.table("Results_Hypred_reverse_ME0_rand.R")
hypred_reverse_SD0_ME.25 <- read.table("Results_Hypred_reverse_ME.25_rand.R")
hypred_reverse_SD0_ME.5 <- read.table("Results_Hypred_reverse_ME.5_rand.R")
hypred_reverse_SD0_ME.75 <- read.table("Results_Hypred_reverse_ME.75_rand.R")

setwd("C:/Users/Barny/Dropbox/KI_Project_4/Results/Uganda")
hypred_uganda_ME0 <- read.table("Uganda_hypred_ME0.R")
hypred_uganda_ME.05 <- read.table("Uganda_hypred_ME.05.R")
hypred_uganda_ME.12 <- read.table("Uganda_hypred_ME.12.R")
hypred_uganda_ME.25 <- read.table("Uganda_hypred_ME.25.R")
hypred_uganda_ME.5 <- read.table("Uganda_hypred_ME.5.R")
hypred_uganda_ME.75 <- read.table("Uganda_hypred_ME.75.R")


### reforming the data frame into long format with the method name as a variable
reform <- function(df_name){
  temp <- df_name
  
  linreg <- temp[,c(1:6)]
  linreg <- linreg %>%
    rename(
      tsens = linreg_tsens,
      esense = linreg_esens,
      t = linreg_t,
      rds = linreg_rds,
      tprevfail = linreg_tprevfail,
      eprevfail = linreg_eprevfail
    )
  linreg$eff <- (100 - linreg$t)/100
  linreg$method <- "linreg"
  
  mss <- temp[,c(7:12)]
  mss <- mss %>%
    rename(
      tsens = mss_tsens,
      esense = mss_esens,
      t = mss_t,
      rds = mss_rds,
      tprevfail = mss_tprevfail,
      eprevfail = mss_eprevfail
    )
  mss$eff <- (100 - mss$t)/100
  mss$method <- "mss"
  
  lrsoe <- temp[,c(13:18)]
  lrsoe <- lrsoe %>%
    rename(
      tsens = lrsoe_tsens,
      esense = lrsoe_esens,
      t = lrsoe_t,
      rds = lrsoe_rds,
      tprevfail = lrsoe_tprevfail,
      eprevfail = lrsoe_eprevfail
    )
  lrsoe$eff <- (100 - lrsoe$t)/100
  lrsoe$method <- "lrsoe"
  
  mincov <- temp[,c(19:24)]
  mincov <- mincov %>%
    rename(
      tsens = mincov_tsens,
      esense = mincov_esens,
      t = mincov_t,
      rds = mincov_rds,
      tprevfail = mincov_tprevfail,
      eprevfail = mincov_eprevfail
    )
  mincov$eff <- (100 - mincov$t)/100
  mincov$method <- "mincov"
  
  mini <- temp[,c(25:30)]
  mini <- mini %>%
    rename(
      tsens = mini_tsens,
      esense = mini_esens,
      t = mini_t,
      rds = mini_rds,
      tprevfail = mini_tprevfail,
      eprevfail = mini_eprevfail
    )
  mini$eff <- (100 - mini$t)/100
  mini$method <- "mini"
  
  temp2 <- rbind(linreg, mss, lrsoe, mincov, mini)
  
  for (i in 1:length(temp2$esense)){
    if (is.na(temp2$esense[i])){
      temp2$esense[i]=1
    }
  }
  return(temp2)
}


binCI <- function(dataset, pool_method, matsize=10, ci_method="clopper-pearson"){
    temp <- reform(dataset)
    temp <- temp[which(temp$method==pool_method),]
  
    sens <- data.frame(BinomCI(sum(temp$esense*temp$eprevfail, na.rm=TRUE),
        sum(temp$eprevfail),
        conf.level=0.95,
        method=ci_method))

    sens.est <- paste0(round(sens$est*100, digits=1), c("%"))
    sens.LCL <- paste0(c("("), round(sens$lwr.ci*100, digits=1), c("%"), c(","))
    sens.UCL <- paste0(round(sens$upr.ci*100, digits=1), c("%"), c(")"))

    sens.summ <- paste(sens.est, sens.LCL, sens.UCL, sep=" ", collapse=NULL)
  
    eff <- data.frame(BinomCI(sum(temp$eff*(matsize^2), na.rm=TRUE),
                            length(temp$t)*(matsize^2),
                            conf.level=0.95,
                            method=ci_method))
    eff.est <- paste0(round(eff$est*100, digits=1), c("%"))
    eff.LCL <- paste0(c("("), round(eff$lwr.ci*100, digits=1), c("%"), c(","))
    eff.UCL <- paste0(round(eff$upr.ci*100, digits=1), c("%"), c(")"))
    
    eff.summ <- paste(eff.est, eff.LCL, eff.UCL, sep=" ", collapse=NULL)
    
    eff.quant <- quantile(temp$eff)
    Q1 <- paste0(c("("), round(eff.quant[2]*100, digits=1), c("%"), c(","))
    Q3 <- paste0(round(eff.quant[4]*100, digits=1), c("%"), c(")"))
    med <- paste0(round(eff.quant[3]*100, digits=1), c("%"))
    eff.quant.summ <- paste(med, Q1, Q3, sep=" ", collapse=NULL)
    
    min <- paste0(min(temp$eff)*100, c("%"), c(","))
    max <- paste0(max(temp$eff)*100, c("%"))
    eff.min <- paste(min, max, sep=" ", collapse=NULL)
    
    rds.mean <- mean(temp$rds)
    rds.stdev <- sd(temp$rds)
    rds.n = length(temp$rds)
    error <- qt(0.975,df=rds.n-1)*rds.stdev/sqrt(rds.n)
    UL <- rds.mean+error
    LL <- rds.mean-error
    
    rds.summ <- paste0(round(rds.mean, digits=1), c(" ("), 
                       round(LL, digits=1), c(", "), round(UL, digits=1), c(")"))
    
    rds.quant <- quantile(temp$rds)
    rds.Q1 <- paste0(c("("),rds.quant[2],c(","))
    rds.med <- rds.quant[3]
    rds.Q3 <- paste0(rds.quant[4], c(")"))
    rds.summ.quant <- paste(rds.med, rds.Q1, rds.Q3, sep=" ", collapse=NULL)
    
    rds.min <- paste(min(temp$rds), max(temp$rds), sep=", ", collapse=NULL)
    
    temp_blank <- matrix(nrow=1,ncol=1)
    
    
    summ <- data.frame(rbind(temp_blank, sens.summ, temp_blank, eff.summ, eff.quant.summ, eff.min, 
                             temp_blank, rds.summ, rds.summ.quant, rds.min))
    
    
    if (pool_method=="linreg"){names(summ)[1] <- "Linreg"}
    if (pool_method=="mss"){names(summ)[1] <- "MSS"}
    if (pool_method=="lrsoe"){names(summ)[1] <- "LRSOE"}
    if (pool_method=="mincov"){names(summ)[1] <- "MiniCov"}
    if (pool_method=="mini"){names(summ)[1] <- "Mini"}
    
return(summ)
}

large <- function(x){
  paste0('{\\Large ', x, '}')
}
bold <- function(x){
  paste0('{\\bfseries ', x, '}')
}
italic <- function(x){
  paste0('{\\emph{ ', x, '}}')
}



final_table <- function(dataset, caption, matsize=10, ci_method="clopper-pearson"){
  temp <- dataset
  
  linreg <- binCI(dataset=temp, pool_method="linreg")
  mss <- binCI(dataset=temp, pool_method="mss")
  lrsoe <- binCI(dataset=temp, pool_method="lrsoe")
  mincov <- binCI(dataset=temp, pool_method="mincov")
  mini <- binCI(dataset=temp, pool_method="mini")
  
  labels <- NULL
  labels[1] <- ""
  labels[2] <- "Mean (95% CI)"
  labels[3] <- ""
  labels[4] <- "Mean (95% CI)"
  labels[5] <- "Median (Q1, Q3)"
  labels[6] <- "Min, Max"
  labels[7] <- ""
  labels[8] <- "Mean (95% CI)"
  labels[9] <- "Median (Q1, Q3)"
  labels[10] <- "Min, Max"
  
  label1 <- NULL
  label1[1] <- "Sensitivity"
  label1[2] <- ""
  label1[3] <- "Efficiency"
  label1[4] <- ""
  label1[5] <- ""
  label1[6] <- ""
  label1[7] <- "Rounds"
  label1[8] <- ""
  label1[9] <- ""
  label1[10] <- ""
  
  result <- data.frame(cbind(label1, labels, linreg, mss, lrsoe, mincov, mini))
  
  names(result)[1] <- ""
  names(result)[2] <- ""
  
  final <- xtable(result, include.rownames=FALSE, caption=caption)
  align(final) <- "rrr|c|c|c|c|c"
  
  return(print(final, include.rownames=FALSE, sanitize.colnames.function = bold,
               caption.placement="top"))
}



final_table(dataset=hypred_AGAIG_SD0_ME0[hypred_AGAIG_SD0_ME0$section=="mid",], 
              caption = "Hypred Mid Tier, AGAIG, SD=0 ME=0, est betas")
final_table(dataset=hypred_AGAIG_SD0_ME.25[hypred_AGAIG_SD0_ME.25$section=="mid",], 
            caption = "Hypred Mid Tier, AGAIG, SD=0 ME=0.25, est betas")
final_table(dataset=hypred_AGAIG_SD0_ME.5[hypred_AGAIG_SD0_ME.5$section=="mid",], 
            caption = "Hypred Mid Tier, AGAIG, SD=0 ME=0.5, est betas")
final_table(dataset=hypred_AGAIG_SD0_ME.75[hypred_AGAIG_SD0_ME.75$section=="mid",], 
            caption = "Hypred Mid Tier, AGAIG, SD=0 ME=0.75, est betas")


final_table(dataset=hypred_AGAIG_ME.05[hypred_AGAIG_ME.05$section=="mid",], 
            caption = "Hypred Mid Tier, AGAIG, SD=1 ME=0.05, est betas")
final_table(dataset=hypred_AGAIG_ME.12[hypred_AGAIG_ME.12$section=="mid",], 
            caption = "Hypred Mid Tier, AGAIG, SD=1 ME=0.12, est betas")
final_table(dataset=hypred_AGAIG_ME.25[hypred_AGAIG_ME.25$section=="mid",], 
            caption = "Hypred Mid Tier, AGAIG, SD=1 ME=0.25, est betas")
final_table(dataset=hypred_AGAIG_ME.5[hypred_AGAIG_ME.5$section=="mid",], 
            caption = "Hypred Mid Tier, AGAIG, SD=1 ME=0.5, est betas")


final_table(dataset=hypred_AGAIG_SD0_ME0[hypred_AGAIG_SD0_ME0$section=="low",], 
            caption = "Hypred low Tier, AGAIG, SD=0 ME=0, est betas")
final_table(dataset=hypred_AGAIG_SD0_ME.25[hypred_AGAIG_SD0_ME.25$section=="low",], 
            caption = "Hypred low Tier, AGAIG, SD=0 ME=0.25, est betas")
final_table(dataset=hypred_AGAIG_SD0_ME.5[hypred_AGAIG_SD0_ME.5$section=="low",], 
            caption = "Hypred low Tier, AGAIG, SD=0 ME=0.5, est betas")
final_table(dataset=hypred_AGAIG_SD0_ME.75[hypred_AGAIG_SD0_ME.75$section=="low",], 
            caption = "Hypred low Tier, AGAIG, SD=0 ME=0.75, est betas")


final_table(dataset=hypred_AGAIG_ME.05[hypred_AGAIG_ME.05$section=="low",], 
            caption = "Hypred low Tier, AGAIG, SD=1 ME=0.05, est betas")
final_table(dataset=hypred_AGAIG_ME.12[hypred_AGAIG_ME.12$section=="low",], 
            caption = "Hypred low Tier, AGAIG, SD=1 ME=0.12, est betas")
final_table(dataset=hypred_AGAIG_ME.25[hypred_AGAIG_ME.25$section=="low",], 
            caption = "Hypred low Tier, AGAIG, SD=1 ME=0.25, est betas")
final_table(dataset=hypred_AGAIG_ME.5[hypred_AGAIG_ME.5$section=="low",], 
            caption = "Hypred low Tier, AGAIG, SD=1 ME=0.5, est betas")

#####################################################################################################

final_table(dataset=hypred_reverse_SD0_ME0[hypred_reverse_SD0_ME0$section=="mid",], 
            caption = "Hypred Mid Tier, Reverse, SD=0 ME=0, est betas")
final_table(dataset=hypred_reverse_SD0_ME.25[hypred_reverse_SD0_ME.25$section=="mid",], 
            caption = "Hypred Mid Tier, Reverse, SD=0 ME=0.25, est betas")
final_table(dataset=hypred_reverse_SD0_ME.5[hypred_reverse_SD0_ME.5$section=="mid",], 
            caption = "Hypred Mid Tier, Reverse, SD=0 ME=0.5, est betas")
final_table(dataset=hypred_reverse_SD0_ME.75[hypred_reverse_SD0_ME.75$section=="mid",], 
            caption = "Hypred Mid Tier, Reverse, SD=0 ME=0.75, est betas")


final_table(dataset=hypred_reverse_ME.05[hypred_reverse_ME.05$section=="mid",], 
            caption = "Hypred Mid Tier, Reverse, SD=1 ME=0.05, est betas")
final_table(dataset=hypred_reverse_ME.12[hypred_reverse_ME.12$section=="mid",], 
            caption = "Hypred Mid Tier, Reverse, SD=1 ME=0.12, est betas")
final_table(dataset=hypred_reverse_ME.25[hypred_reverse_ME.25$section=="mid",], 
            caption = "Hypred Mid Tier, Reverse, SD=1 ME=0.25, est betas")
final_table(dataset=hypred_reverse_ME.5[hypred_reverse_ME.5$section=="mid",], 
            caption = "Hypred Mid Tier, Reverse, SD=1 ME=0.5, est betas")



final_table(dataset=hypred_reverse_SD0_ME0[hypred_reverse_SD0_ME0$section=="low",], 
            caption = "Hypred low Tier, Reverse, SD=0 ME=0, est betas")
final_table(dataset=hypred_reverse_SD0_ME.25[hypred_reverse_SD0_ME.25$section=="low",], 
            caption = "Hypred low Tier, Reverse, SD=0 ME=0.25, est betas")
final_table(dataset=hypred_reverse_SD0_ME.5[hypred_reverse_SD0_ME.5$section=="low",], 
            caption = "Hypred low Tier, Reverse, SD=0 ME=0.5, est betas")
final_table(dataset=hypred_reverse_SD0_ME.75[hypred_reverse_SD0_ME.75$section=="low",], 
            caption = "Hypred low Tier, Reverse, SD=0 ME=0.75, est betas")


final_table(dataset=hypred_reverse_ME.05[hypred_reverse_ME.05$section=="low",], 
            caption = "Hypred low Tier, Reverse, SD=1 ME=0.05, est betas")
final_table(dataset=hypred_reverse_ME.12[hypred_reverse_ME.12$section=="low",], 
            caption = "Hypred low Tier, Reverse, SD=1 ME=0.12, est betas")
final_table(dataset=hypred_reverse_ME.25[hypred_reverse_ME.25$section=="low",], 
            caption = "Hypred low Tier, Reverse, SD=1 ME=0.25, est betas")
final_table(dataset=hypred_reverse_ME.5[hypred_reverse_ME.5$section=="low",], 
            caption = "Hypred low Tier, Reverse, SD=1 ME=0.5, est betas")

            
comb_hypred <- function(dataset, mid_method, low_method){
  temp <- dataset
  
  top <- temp[temp$section=="top",]
  mid <- temp[temp$section=="mid",]
  low <- temp[temp$section=="low",]
  
  if (mid_method=="mss"){
    mid$eprevfail <- mid$mss_eprevfail
    mid$edet <- mid$mss_esens*mid$mss_eprevfail
    mid$t <- mid$mss_t
    mid$rds <- mid$mss_rds
  }
  else if (mid_method=="linreg"){
    mid$eprevfail <- mid$linreg_eprevfail
    mid$edet <- mid$linreg_esens*mid$linreg_eprevfail
    mid$t <- mid$linreg_t
    mid$rds <- mid$linreg_rds
  }
  else if (mid_method=="lrsoe"){
    mid$eprevfail <- mid$lrsoe_eprevfail
    mid$edet <- mid$mss_esens*mid$lrsoe_eprevfail
    mid$t <- mid$lrsoe_t
    mid$rds <- mid$lrsoe_rds
  }
  else if (mid_method=="mincov"){
    mid$eprevfail <- mid$mincov_eprevfail
    mid$edet <- mid$mss_esens*mid$mincov_eprevfail
    mid$t <- mid$mincov_t
    mid$rds <- mid$mincov_rds
  }
  else if (mid_method=="mini"){
    mid$eprevfail <- mid$mini_eprevfail
    mid$edet <- mid$mss_esens*mid$mini_eprevfail
    mid$t <- mid$mini_t
    mid$rds <- mid$mini_rds
  }
  
  if (low_method=="mss"){
    low$eprevfail <- low$mss_eprevfail
    low$edet <- low$mss_esens*low$mss_eprevfail
    low$t <- low$mss_t
    low$rds <- low$mss_rds
  }
  else if (low_method=="linreg"){
    low$eprevfail <- low$linreg_eprevfail
    low$edet <- low$linreg_esens*low$linreg_eprevfail
    low$t <- low$linreg_t
    low$rds <- low$linreg_rds
  }
  else if (low_method=="lrsoe"){
    low$eprevfail <- low$lrsoe_eprevfail
    low$edet <- low$mss_esens*low$lrsoe_eprevfail
    low$t <- low$lrsoe_t
    low$rds <- low$lrsoe_rds
  }
  else if (low_method=="mincov"){
    low$eprevfail <- low$mincov_eprevfail
    low$edet <- low$mss_esens*low$mincov_eprevfail
    low$t <- low$mincov_t
    low$rds <- low$mincov_rds
  }
  else if (low_method=="mini"){
    low$eprevfail <- low$mini_eprevfail
    low$edet <- low$mss_esens*low$mini_eprevfail
    low$t <- low$mini_t
    low$rds <- low$mini_rds
  }
  
  sens_denom <- top$all_eprevfail + sum(low$eprevfail, na.rm=TRUE) + sum(mid$eprevfail, na.rm=TRUE)
  if (is.na(top$all_esens)){
    sens_num <- sum(low$edet, na.rm=TRUE) + sum(mid$edet, na.rm=TRUE)
  }
  else {sens_num <- (top$all_esens*top$all_eprevfail) + sum(low$edet, na.rm=TRUE) + sum(mid$edet, na.rm=TRUE)}
  sens <- round(sens_num/sens_denom, digits=4)
  
  eff_denom <- top$all_t + (length(mid$t)*100) + (length(low$t)*100)
  eff_num <- top$all_t + sum(mid$t) + sum(low$t)
  eff <- (eff_denom - eff_num)/eff_denom
  
  rds_denom <- length(mid$t) + length(low$t) + (top$all_t/100)
  rds_num <- (top$all_t/100) + sum(mid$rds) + sum(low$rds)
  rds <- round(rds_num/rds_denom, digits=1)
  
  return(c(sens, eff, rds))
}


comb_hypred(hypred_AGAIG_SD0_ME0, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_AGAIG_SD0_ME.25, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_AGAIG_SD0_ME.5, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_AGAIG_SD0_ME.75, mid_method="mincov", low_method="mincov")

comb_hypred(hypred_AGAIG_ME.05, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_AGAIG_ME.12, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_AGAIG_ME.25, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_AGAIG_ME.5, mid_method="mincov", low_method="mincov")

debug(comb_hypred)
comb_hypred(hypred_reverse_SD0_ME0, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_reverse_SD0_ME.25, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_reverse_SD0_ME.5, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_reverse_SD0_ME.75, mid_method="mincov", low_method="mincov")

comb_hypred(hypred_reverse_ME.05, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_reverse_ME.12, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_reverse_ME.25, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_reverse_ME.5, mid_method="mincov", low_method="mincov")

comb_hypred(hypred_uganda_ME0, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_uganda_ME.05, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_uganda_ME.12, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_uganda_ME.25, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_uganda_ME.5, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_uganda_ME.75, mid_method="mincov", low_method="mincov")

