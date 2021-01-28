#==============================================================================
# FILENAME: results_anal_hypred.R
# PROJECT: 	Pooled testing in HIV
# PURPOSE: Creates LaTeX code for summary result tables for the HyPred method
# AUTHOR: Adam Brand


# INPUT datasets: the results datasets in the 'Results' folder


# R VERSION: 3.6.1
#==============================================================================
#Notes: this program is different from results_anal in that it creates separate tables for the
#       middle and bottom risk tier groups. Then in the last step this program combines results
#       from all 3 risk tier groups (depending on which method you use for the middle and bottom tier)
#       and outputs the mean of the 3 performance indicators. This is what is reported in tables 1-3
#       for the HyPred method.

# =============================================================================

library(qwraps2)
library(dplyr)
library(rlang)
library(xtable)
library(PropCIs)
library(DescTools)
library(here)

### AGAIG, SD=1
hypred_AGAIG_SD1_ME0 <- readRDS("Results/SimResults/Results_Hypred_AGAIG_SD1_ME0_rand.rds")
hypred_AGAIG_SD1_ME.025 <- readRDS("Results/SimResults/Results_Hypred_AGAIG_SD1_ME.025_rand.rds")
hypred_AGAIG_SD1_ME.05 <- readRDS("Results/SimResults/Results_Hypred_AGAIG_SD1_ME.05_rand.rds")
hypred_AGAIG_SD1_ME.075 <- readRDS("Results/SimResults/Results_Hypred_AGAIG_SD1_ME.075_rand.rds")
hypred_AGAIG_SD1_ME.1 <- readRDS("Results/SimResults/Results_Hypred_AGAIG_SD1_ME.1_rand.rds")
hypred_AGAIG_SD1_ME.125 <- readRDS("Results/SimResults/Results_Hypred_AGAIG_SD1_ME.125_rand.rds")
hypred_AGAIG_SD1_ME.15 <- readRDS("Results/SimResults/Results_Hypred_AGAIG_SD1_ME.15_rand.rds")
hypred_AGAIG_SD1_ME.175 <- readRDS("Results/SimResults/Results_Hypred_AGAIG_SD1_ME.175_rand.rds")
hypred_AGAIG_SD1_ME.2 <- readRDS("Results/SimResults/Results_Hypred_AGAIG_SD1_ME.2_rand.rds")
hypred_AGAIG_SD1_ME.225 <- readRDS("Results/SimResults/Results_Hypred_AGAIG_SD1_ME.225_rand.rds")
hypred_AGAIG_SD1_ME.25 <- readRDS("Results/SimResults/Results_Hypred_AGAIG_SD1_ME.25_rand.rds")
hypred_AGAIG_SD1_ME.5 <- readRDS("Results/SimResults/Results_Hypred_AGAIG_SD1_ME.5_rand.rds")
hypred_AGAIG_SD1_ME.75 <- readRDS("Results/SimResults/Results_Hypred_AGAIG_SD1_ME.75_rand.rds")


# Reverse, SD=1
hypred_reverse_SD1_ME0 <- readRDS("Results/SimResults/Results_Hypred_reverse_SD1_ME0_rand.rds")
hypred_reverse_SD1_ME.025 <- readRDS("Results/SimResults/Results_Hypred_reverse_SD1_ME.025_rand.rds")
hypred_reverse_SD1_ME.05 <- readRDS("Results/SimResults/Results_Hypred_reverse_SD1_ME.05_rand.rds")
hypred_reverse_SD1_ME.075 <- readRDS("Results/SimResults/Results_Hypred_reverse_SD1_ME.075_rand.rds")
hypred_reverse_SD1_ME.1 <- readRDS("Results/SimResults/Results_Hypred_reverse_SD1_ME.1_rand.rds")
hypred_reverse_SD1_ME.125 <- readRDS("Results/SimResults/Results_Hypred_reverse_SD1_ME.125_rand.rds")
hypred_reverse_SD1_ME.15 <- readRDS("Results/SimResults/Results_Hypred_reverse_SD1_ME.15_rand.rds")
hypred_reverse_SD1_ME.175 <- readRDS("Results/SimResults/Results_Hypred_reverse_SD1_ME.175_rand.rds")
hypred_reverse_SD1_ME.2 <- readRDS("Results/SimResults/Results_Hypred_reverse_SD1_ME.2_rand.rds")
hypred_reverse_SD1_ME.225 <- readRDS("Results/SimResults/Results_Hypred_reverse_SD1_ME.225_rand.rds")
hypred_reverse_SD1_ME.25 <- readRDS("Results/SimResults/Results_Hypred_reverse_SD1_ME.25_rand.rds")
hypred_reverse_SD1_ME.5 <- readRDS("Results/SimResults/Results_Hypred_reverse_SD1_ME.5_rand.rds")
hypred_reverse_SD1_ME.75 <- readRDS("Results/SimResults/Results_Hypred_reverse_SD1_ME.75_rand.rds")


# No association, SD=1
hypred_noassoc_SD1_ME0 <- readRDS("Results/SimResults/Results_Hypred_noassoc_SD1_ME0_rand.rds")
hypred_noassoc_SD1_ME.025 <- readRDS("Results/SimResults/Results_Hypred_noassoc_SD1_ME.025_rand.rds")
hypred_noassoc_SD1_ME.05 <- readRDS("Results/SimResults/Results_Hypred_noassoc_SD1_ME.05_rand.rds")
hypred_noassoc_SD1_ME.075 <- readRDS("Results/SimResults/Results_Hypred_noassoc_SD1_ME.075_rand.rds")
hypred_noassoc_SD1_ME.1 <- readRDS("Results/SimResults/Results_Hypred_noassoc_SD1_ME.1_rand.rds")
hypred_noassoc_SD1_ME.125 <- readRDS("Results/SimResults/Results_Hypred_noassoc_SD1_ME.125_rand.rds")
hypred_noassoc_SD1_ME.15 <- readRDS("Results/SimResults/Results_Hypred_noassoc_SD1_ME.15_rand.rds")
hypred_noassoc_SD1_ME.175 <- readRDS("Results/SimResults/Results_Hypred_noassoc_SD1_ME.175_rand.rds")
hypred_noassoc_SD1_ME.2 <- readRDS("Results/SimResults/Results_Hypred_noassoc_SD1_ME.2_rand.rds")
hypred_noassoc_SD1_ME.225 <- readRDS("Results/SimResults/Results_Hypred_noassoc_SD1_ME.225_rand.rds")
hypred_noassoc_SD1_ME.25 <- readRDS("Results/SimResults/Results_Hypred_noassoc_SD1_ME.25_rand.rds")
hypred_noassoc_SD1_ME.5 <- readRDS("Results/SimResults/Results_Hypred_noassoc_SD1_ME.5_rand.rds")
hypred_noassoc_SD1_ME.75 <- readRDS("Results/SimResults/Results_Hypred_noassoc_SD1_ME.75_rand.rds")


# Misspecified, SD=1
hypred_misspec_SD1_ME0 <- readRDS("Results/SimResults/Results_Hypred_misspec_SD1_ME0_rand.rds")
hypred_misspec_SD1_ME.025 <- readRDS("Results/SimResults/Results_Hypred_misspec_SD1_ME.025_rand.rds")
hypred_misspec_SD1_ME.05 <- readRDS("Results/SimResults/Results_Hypred_misspec_SD1_ME.05_rand.rds")
hypred_misspec_SD1_ME.075 <- readRDS("Results/SimResults/Results_Hypred_misspec_SD1_ME.075_rand.rds")
hypred_misspec_SD1_ME.1 <- readRDS("Results/SimResults/Results_Hypred_misspec_SD1_ME.1_rand.rds")
hypred_misspec_SD1_ME.125 <- readRDS("Results/SimResults/Results_Hypred_misspec_SD1_ME.125_rand.rds")
hypred_misspec_SD1_ME.15 <- readRDS("Results/SimResults/Results_Hypred_misspec_SD1_ME.15_rand.rds")
hypred_misspec_SD1_ME.175 <- readRDS("Results/SimResults/Results_Hypred_misspec_SD1_ME.175_rand.rds")
hypred_misspec_SD1_ME.2 <- readRDS("Results/SimResults/Results_Hypred_misspec_SD1_ME.2_rand.rds")
hypred_misspec_SD1_ME.225 <- readRDS("Results/SimResults/Results_Hypred_misspec_SD1_ME.225_rand.rds")
hypred_misspec_SD1_ME.25 <- readRDS("Results/SimResults/Results_Hypred_misspec_SD1_ME.25_rand.rds")
hypred_misspec_SD1_ME.5 <- readRDS("Results/SimResults/Results_Hypred_misspec_SD1_ME.5_rand.rds")
hypred_misspec_SD1_ME.75 <- readRDS("Results/SimResults/Results_Hypred_misspec_SD1_ME.75_rand.rds")


                ##############################

### AGAIG, SD=0
hypred_AGAIG_SD0_ME0 <- readRDS("Results/SimResults/Results_Hypred_AGAIG_SD0_ME0_rand.rds")
hypred_AGAIG_SD0_ME.025 <- readRDS("Results/SimResults/Results_Hypred_AGAIG_SD0_ME.025_rand.rds")
hypred_AGAIG_SD0_ME.05 <- readRDS("Results/SimResults/Results_Hypred_AGAIG_SD0_ME.05_rand.rds")
hypred_AGAIG_SD0_ME.075 <- readRDS("Results/SimResults/Results_Hypred_AGAIG_SD0_ME.075_rand.rds")
hypred_AGAIG_SD0_ME.1 <- readRDS("Results/SimResults/Results_Hypred_AGAIG_SD0_ME.1_rand.rds")
hypred_AGAIG_SD0_ME.125 <- readRDS("Results/SimResults/Results_Hypred_AGAIG_SD0_ME.125_rand.rds")
hypred_AGAIG_SD0_ME.15 <- readRDS("Results/SimResults/Results_Hypred_AGAIG_SD0_ME.15_rand.rds")
hypred_AGAIG_SD0_ME.175 <- readRDS("Results/SimResults/Results_Hypred_AGAIG_SD0_ME.175_rand.rds")
hypred_AGAIG_SD0_ME.2 <- readRDS("Results/SimResults/Results_Hypred_AGAIG_SD0_ME.2_rand.rds")
hypred_AGAIG_SD0_ME.225 <- readRDS("Results/SimResults/Results_Hypred_AGAIG_SD0_ME.225_rand.rds")
hypred_AGAIG_SD0_ME.25 <- readRDS("Results/SimResults/Results_Hypred_AGAIG_SD0_ME.25_rand.rds")
hypred_AGAIG_SD0_ME.5 <- readRDS("Results/SimResults/Results_Hypred_AGAIG_SD0_ME.5_rand.rds")
hypred_AGAIG_SD0_ME.75 <- readRDS("Results/SimResults/Results_Hypred_AGAIG_SD0_ME.75_rand.rds")


# Reverse, SD=0
hypred_reverse_SD0_ME0 <- readRDS("Results/SimResults/Results_Hypred_reverse_SD0_ME0_rand.rds")
hypred_reverse_SD0_ME.025 <- readRDS("Results/SimResults/Results_Hypred_reverse_SD0_ME.025_rand.rds")
hypred_reverse_SD0_ME.05 <- readRDS("Results/SimResults/Results_Hypred_reverse_SD0_ME.05_rand.rds")
hypred_reverse_SD0_ME.075 <- readRDS("Results/SimResults/Results_Hypred_reverse_SD0_ME.075_rand.rds")
hypred_reverse_SD0_ME.1 <- readRDS("Results/SimResults/Results_Hypred_reverse_SD0_ME.1_rand.rds")
hypred_reverse_SD0_ME.125 <- readRDS("Results/SimResults/Results_Hypred_reverse_SD0_ME.125_rand.rds")
hypred_reverse_SD0_ME.15 <- readRDS("Results/SimResults/Results_Hypred_reverse_SD0_ME.15_rand.rds")
hypred_reverse_SD0_ME.175 <- readRDS("Results/SimResults/Results_Hypred_reverse_SD0_ME.175_rand.rds")
hypred_reverse_SD0_ME.2 <- readRDS("Results/SimResults/Results_Hypred_reverse_SD0_ME.2_rand.rds")
hypred_reverse_SD0_ME.225 <- readRDS("Results/SimResults/Results_Hypred_reverse_SD0_ME.225_rand.rds")
hypred_reverse_SD0_ME.25 <- readRDS("Results/SimResults/Results_Hypred_reverse_SD0_ME.25_rand.rds")
hypred_reverse_SD0_ME.5 <- readRDS("Results/SimResults/Results_Hypred_reverse_SD0_ME.5_rand.rds")
hypred_reverse_SD0_ME.75 <- readRDS("Results/SimResults/Results_Hypred_reverse_SD0_ME.75_rand.rds")


# No association, SD=0
hypred_noassoc_SD0_ME0 <- readRDS("Results/SimResults/Results_Hypred_noassoc_SD0_ME0_rand.rds")
hypred_noassoc_SD0_ME.025 <- readRDS("Results/SimResults/Results_Hypred_noassoc_SD0_ME.025_rand.rds")
hypred_noassoc_SD0_ME.05 <- readRDS("Results/SimResults/Results_Hypred_noassoc_SD0_ME.05_rand.rds")
hypred_noassoc_SD0_ME.075 <- readRDS("Results/SimResults/Results_Hypred_noassoc_SD0_ME.075_rand.rds")
hypred_noassoc_SD0_ME.1 <- readRDS("Results/SimResults/Results_Hypred_noassoc_SD0_ME.1_rand.rds")
hypred_noassoc_SD0_ME.125 <- readRDS("Results/SimResults/Results_Hypred_noassoc_SD0_ME.125_rand.rds")
hypred_noassoc_SD0_ME.15 <- readRDS("Results/SimResults/Results_Hypred_noassoc_SD0_ME.15_rand.rds")
hypred_noassoc_SD0_ME.175 <- readRDS("Results/SimResults/Results_Hypred_noassoc_SD0_ME.175_rand.rds")
hypred_noassoc_SD0_ME.2 <- readRDS("Results/SimResults/Results_Hypred_noassoc_SD0_ME.2_rand.rds")
hypred_noassoc_SD0_ME.225 <- readRDS("Results/SimResults/Results_Hypred_noassoc_SD0_ME.225_rand.rds")
hypred_noassoc_SD0_ME.25 <- readRDS("Results/SimResults/Results_Hypred_noassoc_SD0_ME.25_rand.rds")
hypred_noassoc_SD0_ME.5 <- readRDS("Results/SimResults/Results_Hypred_noassoc_SD0_ME.5_rand.rds")
hypred_noassoc_SD0_ME.75 <- readRDS("Results/SimResults/Results_Hypred_noassoc_SD0_ME.75_rand.rds")


# Misspecified, SD=0
hypred_misspec_SD0_ME0 <- readRDS("Results/SimResults/Results_Hypred_misspec_SD0_ME0_rand.rds")
hypred_misspec_SD0_ME.025 <- readRDS("Results/SimResults/Results_Hypred_misspec_SD0_ME.025_rand.rds")
hypred_misspec_SD0_ME.05 <- readRDS("Results/SimResults/Results_Hypred_misspec_SD0_ME.05_rand.rds")
hypred_misspec_SD0_ME.075 <- readRDS("Results/SimResults/Results_Hypred_misspec_SD0_ME.075_rand.rds")
hypred_misspec_SD0_ME.1 <- readRDS("Results/SimResults/Results_Hypred_misspec_SD0_ME.1_rand.rds")
hypred_misspec_SD0_ME.125 <- readRDS("Results/SimResults/Results_Hypred_misspec_SD0_ME.125_rand.rds")
hypred_misspec_SD0_ME.15 <- readRDS("Results/SimResults/Results_Hypred_misspec_SD0_ME.15_rand.rds")
hypred_misspec_SD0_ME.175 <- readRDS("Results/SimResults/Results_Hypred_misspec_SD0_ME.175_rand.rds")
hypred_misspec_SD0_ME.2 <- readRDS("Results/SimResults/Results_Hypred_misspec_SD0_ME.2_rand.rds")
hypred_misspec_SD0_ME.225 <- readRDS("Results/SimResults/Results_Hypred_misspec_SD0_ME.225_rand.rds")
hypred_misspec_SD0_ME.25 <- readRDS("Results/SimResults/Results_Hypred_misspec_SD0_ME.25_rand.rds")
hypred_misspec_SD0_ME.5 <- readRDS("Results/SimResults/Results_Hypred_misspec_SD0_ME.5_rand.rds")
hypred_misspec_SD0_ME.75 <- readRDS("Results/SimResults/Results_Hypred_misspec_SD0_ME.75_rand.rds")



####################################################################################

hypred_uganda_ME0 <- readRDS("Results/UgandaResults/Uganda_hypred_ME0.rds")
hypred_uganda_ME.025 <- readRDS("Results/UgandaResults/Uganda_hypred_ME.025.rds")
hypred_uganda_ME.05 <- readRDS("Results/UgandaResults/Uganda_hypred_ME.05.rds")
hypred_uganda_ME.075 <- readRDS("Results/UgandaResults/Uganda_hypred_ME.075.rds")
hypred_uganda_ME.1 <- readRDS("Results/UgandaResults/Uganda_hypred_ME.1.rds")
hypred_uganda_ME.125 <- readRDS("Results/UgandaResults/Uganda_hypred_ME.125.rds")
hypred_uganda_ME.15 <- readRDS("Results/UgandaResults/Uganda_hypred_ME.15.rds")
hypred_uganda_ME.175 <- readRDS("Results/UgandaResults/Uganda_hypred_ME.175.rds")
hypred_uganda_ME.2 <- readRDS("Results/UgandaResults/Uganda_hypred_ME.2.rds")
hypred_uganda_ME.225 <- readRDS("Results/UgandaResults/Uganda_hypred_ME.225.rds")
hypred_uganda_ME.25 <- readRDS("Results/UgandaResults/Uganda_hypred_ME.25.rds")
hypred_uganda_ME.5 <- readRDS("Results/UgandaResults/Uganda_hypred_ME.5.rds")
hypred_uganda_ME.75 <- readRDS("Results/UgandaResults/Uganda_hypred_ME.75.rds")


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

# computing 95% CIs for the 3 performance criteria
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
    if (pool_method=="mincov"){names(summ)[1] <- "MiniPred"}
    if (pool_method=="mini"){names(summ)[1] <- "Mini+alg"}
    
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


# creates results tables for each of the 3 tiers of risk groups separately for each method within the tier
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

# the below result using final_table outputs latex code for a table of summary results for 
# the mid and low tiers separately. Based on these simulations we decided which method to use for
# each tier on the real Uganda data.
# combined results are below this code.

# NO LONGER INCLUDED IN THE MANUSCRIPT; COMMENTED OUT

# ## AGAIG, SD=0, mid tier
# final_table(dataset=hypred_AGAIG_SD0_ME0[hypred_AGAIG_SD0_ME0$section=="mid",], 
#             caption = "Hypred Mid Tier, AGAIG, SD=0 ME=0, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD0_ME.025[hypred_AGAIG_SD0_ME.025$section=="mid",], 
#             caption = "Hypred Mid Tier, AGAIG, SD=0 ME=0.025, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD0_ME.05[hypred_AGAIG_SD0_ME.05$section=="mid",], 
#             caption = "Hypred Mid Tier, AGAIG, SD=0 ME=0.05, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD0_ME.075[hypred_AGAIG_SD0_ME.075$section=="mid",], 
#             caption = "Hypred Mid Tier, AGAIG, SD=0 ME=0.075, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD0_ME.1[hypred_AGAIG_SD0_ME.1$section=="mid",], 
#             caption = "Hypred Mid Tier, AGAIG, SD=0 ME=0.1, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD0_ME.125[hypred_AGAIG_SD0_ME.125$section=="mid",], 
#             caption = "Hypred Mid Tier, AGAIG, SD=0 ME=0.125, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD0_ME.15[hypred_AGAIG_SD0_ME.15$section=="mid",], 
#             caption = "Hypred Mid Tier, AGAIG, SD=0 ME=0.15, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD0_ME.175[hypred_AGAIG_SD0_ME.175$section=="mid",], 
#             caption = "Hypred Mid Tier, AGAIG, SD=0 ME=0.175, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD0_ME.2[hypred_AGAIG_SD0_ME.2$section=="mid",], 
#             caption = "Hypred Mid Tier, AGAIG, SD=0 ME=0.2, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD0_ME.225[hypred_AGAIG_SD0_ME.225$section=="mid",], 
#             caption = "Hypred Mid Tier, AGAIG, SD=0 ME=0.225, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD0_ME.25[hypred_AGAIG_SD0_ME.25$section=="mid",], 
#             caption = "Hypred Mid Tier, AGAIG, SD=0 ME=0.25, Estimated betas")
# 
# final_table(dataset=hypred_AGAIG_SD0_ME.5[hypred_AGAIG_SD0_ME.5$section=="mid",], 
#             caption = "Hypred Mid Tier, AGAIG, SD=0 ME=0.5, Estimated betas")
# 
# final_table(dataset=hypred_AGAIG_SD0_ME.75[hypred_AGAIG_SD0_ME.75$section=="mid",], 
#             caption = "Hypred Mid Tier, AGAIG, SD=0 ME=0.75, Estimated betas")
# 
# 
# 
# # AGAIG, SD=1, mid tier
# final_table(dataset=hypred_AGAIG_SD1_ME0[hypred_AGAIG_SD1_ME0$section=="mid",], 
#             caption = "Hypred Mid Tier, AGAIG, SD=1 ME=0, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD1_ME.025[hypred_AGAIG_SD1_ME.025$section=="mid",], 
#             caption = "Hypred Mid Tier, AGAIG, SD=1 ME=0.025, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD1_ME.05[hypred_AGAIG_SD1_ME.05$section=="mid",], 
#             caption = "Hypred Mid Tier, AGAIG, SD=1 ME=0.05, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD1_ME.075[hypred_AGAIG_SD1_ME.075$section=="mid",], 
#             caption = "Hypred Mid Tier, AGAIG, SD=1 ME=0.075, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD1_ME.1[hypred_AGAIG_SD1_ME.1$section=="mid",], 
#             caption = "Hypred Mid Tier, AGAIG, SD=1 ME=0.1, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD1_ME.125[hypred_AGAIG_SD1_ME.125$section=="mid",], 
#             caption = "Hypred Mid Tier, AGAIG, SD=1 ME=0.125, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD1_ME.15[hypred_AGAIG_SD1_ME.15$section=="mid",], 
#             caption = "Hypred Mid Tier, AGAIG, SD=1 ME=0.15, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD1_ME.175[hypred_AGAIG_SD1_ME.175$section=="mid",], 
#             caption = "Hypred Mid Tier, AGAIG, SD=1 ME=0.175, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD1_ME.2[hypred_AGAIG_SD1_ME.2$section=="mid",], 
#             caption = "Hypred Mid Tier, AGAIG, SD=1 ME=0.2, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD1_ME.225[hypred_AGAIG_SD1_ME.225$section=="mid",], 
#             caption = "Hypred Mid Tier, AGAIG, SD=1 ME=0.225, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD1_ME.25[hypred_AGAIG_SD1_ME.25$section=="mid",], 
#             caption = "Hypred Mid Tier, AGAIG, SD=1 ME=0.25, Estimated betas")
# 
# final_table(dataset=hypred_AGAIG_SD1_ME.5[hypred_AGAIG_SD1_ME.5$section=="mid",], 
#             caption = "Hypred Mid Tier, AGAIG, SD=1 ME=0.5, Estimated betas")
# 
# final_table(dataset=hypred_AGAIG_SD1_ME.75[hypred_AGAIG_SD1_ME.75$section=="mid",], 
#             caption = "Hypred Mid Tier, AGAIG, SD=1 ME=0.75, Estimated betas")
# 
# 
# ## AGAIG, SD=0, low tier
# final_table(dataset=hypred_AGAIG_SD0_ME0[hypred_AGAIG_SD0_ME0$section=="low",], 
#             caption = "Hypred Low Tier, AGAIG, SD=0 ME=0, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD0_ME.025[hypred_AGAIG_SD0_ME.025$section=="low",], 
#             caption = "Hypred Low Tier, AGAIG, SD=0 ME=0.025, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD0_ME.05[hypred_AGAIG_SD0_ME.05$section=="low",], 
#             caption = "Hypred Low Tier, AGAIG, SD=0 ME=0.05, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD0_ME.075[hypred_AGAIG_SD0_ME.075$section=="low",], 
#             caption = "Hypred Low Tier, AGAIG, SD=0 ME=0.075, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD0_ME.1[hypred_AGAIG_SD0_ME.1$section=="low",], 
#             caption = "Hypred Low Tier, AGAIG, SD=0 ME=0.1, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD0_ME.125[hypred_AGAIG_SD0_ME.125$section=="low",], 
#             caption = "Hypred Low Tier, AGAIG, SD=0 ME=0.125, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD0_ME.15[hypred_AGAIG_SD0_ME.15$section=="low",], 
#             caption = "Hypred Low Tier, AGAIG, SD=0 ME=0.15, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD0_ME.175[hypred_AGAIG_SD0_ME.175$section=="low",], 
#             caption = "Hypred Low Tier, AGAIG, SD=0 ME=0.175, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD0_ME.2[hypred_AGAIG_SD0_ME.2$section=="low",], 
#             caption = "Hypred Low Tier, AGAIG, SD=0 ME=0.2, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD0_ME.225[hypred_AGAIG_SD0_ME.225$section=="low",], 
#             caption = "Hypred Low Tier, AGAIG, SD=0 ME=0.225, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD0_ME.25[hypred_AGAIG_SD0_ME.25$section=="low",], 
#             caption = "Hypred Low Tier, AGAIG, SD=0 ME=0.25, Estimated betas")
# 
# final_table(dataset=hypred_AGAIG_SD0_ME.5[hypred_AGAIG_SD0_ME.5$section=="low",], 
#             caption = "Hypred Low Tier, AGAIG, SD=0 ME=0.5, Estimated betas")
# 
# final_table(dataset=hypred_AGAIG_SD0_ME.75[hypred_AGAIG_SD0_ME.75$section=="low",], 
#             caption = "Hypred Low Tier, AGAIG, SD=0 ME=0.75, Estimated betas")
# 
# 
# ## AGAIG, SD=1, low tier
# final_table(dataset=hypred_AGAIG_SD1_ME0[hypred_AGAIG_SD1_ME0$section=="low",], 
#             caption = "Hypred Low Tier, AGAIG, SD=1 ME=0, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD1_ME.025[hypred_AGAIG_SD1_ME.025$section=="low",], 
#             caption = "Hypred Low Tier, AGAIG, SD=1 ME=0.025, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD1_ME.05[hypred_AGAIG_SD1_ME.05$section=="low",], 
#             caption = "Hypred Low Tier, AGAIG, SD=1 ME=0.05, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD1_ME.075[hypred_AGAIG_SD1_ME.075$section=="low",], 
#             caption = "Hypred Low Tier, AGAIG, SD=1 ME=0.075, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD1_ME.1[hypred_AGAIG_SD1_ME.1$section=="low",], 
#             caption = "Hypred Low Tier, AGAIG, SD=1 ME=0.1, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD1_ME.125[hypred_AGAIG_SD1_ME.125$section=="low",], 
#             caption = "Hypred Low Tier, AGAIG, SD=1 ME=0.125, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD1_ME.15[hypred_AGAIG_SD1_ME.15$section=="low",], 
#             caption = "Hypred Low Tier, AGAIG, SD=1 ME=0.15, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD1_ME.175[hypred_AGAIG_SD1_ME.175$section=="low",], 
#             caption = "Hypred Low Tier, AGAIG, SD=1 ME=0.175, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD1_ME.2[hypred_AGAIG_SD1_ME.2$section=="low",], 
#             caption = "Hypred Low Tier, AGAIG, SD=1 ME=0.2, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD1_ME.225[hypred_AGAIG_SD1_ME.225$section=="low",], 
#             caption = "Hypred Low Tier, AGAIG, SD=1 ME=0.225, Estimated Betas")
# 
# final_table(dataset=hypred_AGAIG_SD1_ME.25[hypred_AGAIG_SD1_ME.25$section=="low",], 
#             caption = "Hypred Low Tier, AGAIG, SD=1 ME=0.25, Estimated betas")
# 
# final_table(dataset=hypred_AGAIG_SD1_ME.5[hypred_AGAIG_SD1_ME.5$section=="low",], 
#             caption = "Hypred Low Tier, AGAIG, SD=1 ME=0.5, Estimated betas")
# 
# final_table(dataset=hypred_AGAIG_SD1_ME.75[hypred_AGAIG_SD1_ME.75$section=="low",], 
#             caption = "Hypred Low Tier, AGAIG, SD=1 ME=0.75, Estimated betas")
# 
# #####################################################################################################
# 
# ## Reverse, SD=0, mid tier
# final_table(dataset=hypred_reverse_SD0_ME0[hypred_reverse_SD0_ME0$section=="mid",], 
#             caption = "Hypred Mid Tier, Reverse, SD=0 ME=0, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD0_ME.025[hypred_reverse_SD0_ME.025$section=="mid",], 
#             caption = "Hypred Mid Tier, Reverse, SD=0 ME=0.025, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD0_ME.05[hypred_reverse_SD0_ME.05$section=="mid",], 
#             caption = "Hypred Mid Tier, Reverse, SD=0 ME=0.05, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD0_ME.075[hypred_reverse_SD0_ME.075$section=="mid",], 
#             caption = "Hypred Mid Tier, Reverse, SD=0 ME=0.075, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD0_ME.1[hypred_reverse_SD0_ME.1$section=="mid",], 
#             caption = "Hypred Mid Tier, Reverse, SD=0 ME=0.1, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD0_ME.125[hypred_reverse_SD0_ME.125$section=="mid",], 
#             caption = "Hypred Mid Tier, Reverse, SD=0 ME=0.125, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD0_ME.15[hypred_reverse_SD0_ME.15$section=="mid",], 
#             caption = "Hypred Mid Tier, Reverse, SD=0 ME=0.15, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD0_ME.175[hypred_reverse_SD0_ME.175$section=="mid",], 
#             caption = "Hypred Mid Tier, Reverse, SD=0 ME=0.175, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD0_ME.2[hypred_reverse_SD0_ME.2$section=="mid",], 
#             caption = "Hypred Mid Tier, Reverse, SD=0 ME=0.2, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD0_ME.225[hypred_reverse_SD0_ME.225$section=="mid",], 
#             caption = "Hypred Mid Tier, Reverse, SD=0 ME=0.225, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD0_ME.25[hypred_reverse_SD0_ME.25$section=="mid",], 
#             caption = "Hypred Mid Tier, Reverse, SD=0 ME=0.25, Estimated betas")
# 
# final_table(dataset=hypred_reverse_SD0_ME.5[hypred_reverse_SD0_ME.5$section=="mid",], 
#             caption = "Hypred Mid Tier, Reverse, SD=0 ME=0.5, Estimated betas")
# 
# final_table(dataset=hypred_reverse_SD0_ME.75[hypred_reverse_SD0_ME.75$section=="mid",], 
#             caption = "Hypred Mid Tier, Reverse, SD=0 ME=0.75, Estimated betas")
# 
# 
# 
# # Reverse, SD=1, mid tier
# final_table(dataset=hypred_reverse_SD1_ME0[hypred_reverse_SD1_ME0$section=="mid",], 
#             caption = "Hypred Mid Tier, Reverse, SD=1 ME=0, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD1_ME.025[hypred_reverse_SD1_ME.025$section=="mid",], 
#             caption = "Hypred Mid Tier, Reverse, SD=1 ME=0.025, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD1_ME.05[hypred_reverse_SD1_ME.05$section=="mid",], 
#             caption = "Hypred Mid Tier, Reverse, SD=1 ME=0.05, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD1_ME.075[hypred_reverse_SD1_ME.075$section=="mid",], 
#             caption = "Hypred Mid Tier, Reverse, SD=1 ME=0.075, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD1_ME.1[hypred_reverse_SD1_ME.1$section=="mid",], 
#             caption = "Hypred Mid Tier, Reverse, SD=1 ME=0.1, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD1_ME.125[hypred_reverse_SD1_ME.125$section=="mid",], 
#             caption = "Hypred Mid Tier, Reverse, SD=1 ME=0.125, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD1_ME.15[hypred_reverse_SD1_ME.15$section=="mid",], 
#             caption = "Hypred Mid Tier, Reverse, SD=1 ME=0.15, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD1_ME.175[hypred_reverse_SD1_ME.175$section=="mid",], 
#             caption = "Hypred Mid Tier, Reverse, SD=1 ME=0.175, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD1_ME.2[hypred_reverse_SD1_ME.2$section=="mid",], 
#             caption = "Hypred Mid Tier, Reverse, SD=1 ME=0.2, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD1_ME.225[hypred_reverse_SD1_ME.225$section=="mid",], 
#             caption = "Hypred Mid Tier, Reverse, SD=1 ME=0.225, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD1_ME.25[hypred_reverse_SD1_ME.25$section=="mid",], 
#             caption = "Hypred Mid Tier, Reverse, SD=1 ME=0.25, Estimated betas")
# 
# final_table(dataset=hypred_reverse_SD1_ME.5[hypred_reverse_SD1_ME.5$section=="mid",], 
#             caption = "Hypred Mid Tier, Reverse, SD=1 ME=0.5, Estimated betas")
# 
# final_table(dataset=hypred_reverse_SD1_ME.75[hypred_reverse_SD1_ME.75$section=="mid",], 
#             caption = "Hypred Mid Tier, Reverse, SD=1 ME=0.75, Estimated betas")
# 
# 
# ## Reverse, SD=0, low tier
# final_table(dataset=hypred_reverse_SD0_ME0[hypred_reverse_SD0_ME0$section=="low",], 
#             caption = "Hypred Low Tier, Reverse, SD=0 ME=0, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD0_ME.025[hypred_reverse_SD0_ME.025$section=="low",], 
#             caption = "Hypred Low Tier, Reverse, SD=0 ME=0.025, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD0_ME.05[hypred_reverse_SD0_ME.05$section=="low",], 
#             caption = "Hypred Low Tier, Reverse, SD=0 ME=0.05, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD0_ME.075[hypred_reverse_SD0_ME.075$section=="low",], 
#             caption = "Hypred Low Tier, Reverse, SD=0 ME=0.075, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD0_ME.1[hypred_reverse_SD0_ME.1$section=="low",], 
#             caption = "Hypred Low Tier, Reverse, SD=0 ME=0.1, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD0_ME.125[hypred_reverse_SD0_ME.125$section=="low",], 
#             caption = "Hypred Low Tier, Reverse, SD=0 ME=0.125, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD0_ME.15[hypred_reverse_SD0_ME.15$section=="low",], 
#             caption = "Hypred Low Tier, Reverse, SD=0 ME=0.15, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD0_ME.175[hypred_reverse_SD0_ME.175$section=="low",], 
#             caption = "Hypred Low Tier, Reverse, SD=0 ME=0.175, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD0_ME.2[hypred_reverse_SD0_ME.2$section=="low",], 
#             caption = "Hypred Low Tier, Reverse, SD=0 ME=0.2, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD0_ME.225[hypred_reverse_SD0_ME.225$section=="low",], 
#             caption = "Hypred Low Tier, Reverse, SD=0 ME=0.225, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD0_ME.25[hypred_reverse_SD0_ME.25$section=="low",], 
#             caption = "Hypred Low Tier, Reverse, SD=0 ME=0.25, Estimated betas")
# 
# final_table(dataset=hypred_reverse_SD0_ME.5[hypred_reverse_SD0_ME.5$section=="low",], 
#             caption = "Hypred Low Tier, Reverse, SD=0 ME=0.5, Estimated betas")
# 
# final_table(dataset=hypred_reverse_SD0_ME.75[hypred_reverse_SD0_ME.75$section=="low",], 
#             caption = "Hypred Low Tier, Reverse, SD=0 ME=0.75, Estimated betas")
# 
# 
# ## Reverse, SD=1, low tier
# final_table(dataset=hypred_reverse_SD1_ME0[hypred_reverse_SD1_ME0$section=="low",], 
#             caption = "Hypred Low Tier, Reverse, SD=1 ME=0, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD1_ME.025[hypred_reverse_SD1_ME.025$section=="low",], 
#             caption = "Hypred Low Tier, Reverse, SD=1 ME=0.025, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD1_ME.05[hypred_reverse_SD1_ME.05$section=="low",], 
#             caption = "Hypred Low Tier, Reverse, SD=1 ME=0.05, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD1_ME.075[hypred_reverse_SD1_ME.075$section=="low",], 
#             caption = "Hypred Low Tier, Reverse, SD=1 ME=0.075, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD1_ME.1[hypred_reverse_SD1_ME.1$section=="low",], 
#             caption = "Hypred Low Tier, Reverse, SD=1 ME=0.1, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD1_ME.125[hypred_reverse_SD1_ME.125$section=="low",], 
#             caption = "Hypred Low Tier, Reverse, SD=1 ME=0.125, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD1_ME.15[hypred_reverse_SD1_ME.15$section=="low",], 
#             caption = "Hypred Low Tier, Reverse, SD=1 ME=0.15, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD1_ME.175[hypred_reverse_SD1_ME.175$section=="low",], 
#             caption = "Hypred Low Tier, Reverse, SD=1 ME=0.175, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD1_ME.2[hypred_reverse_SD1_ME.2$section=="low",], 
#             caption = "Hypred Low Tier, Reverse, SD=1 ME=0.2, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD1_ME.225[hypred_reverse_SD1_ME.225$section=="low",], 
#             caption = "Hypred Low Tier, Reverse, SD=1 ME=0.225, Estimated Betas")
# 
# final_table(dataset=hypred_reverse_SD1_ME.25[hypred_reverse_SD1_ME.25$section=="low",], 
#             caption = "Hypred Low Tier, Reverse, SD=1 ME=0.25, Estimated betas")
# 
# final_table(dataset=hypred_reverse_SD1_ME.5[hypred_reverse_SD1_ME.5$section=="low",], 
#             caption = "Hypred Low Tier, Reverse, SD=1 ME=0.5, Estimated betas")
# 
# final_table(dataset=hypred_reverse_SD1_ME.75[hypred_reverse_SD1_ME.75$section=="low",], 
#             caption = "Hypred Low Tier, Reverse, SD=1 ME=0.75, Estimated betas")
# 
# 
# ###########################################################################################
# 
# ## No association, SD=0, mid tier
# final_table(dataset=hypred_noassoc_SD0_ME0[hypred_noassoc_SD0_ME0$section=="mid",], 
#             caption = "Hypred Mid Tier, No Association, SD=0 ME=0, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD0_ME.025[hypred_noassoc_SD0_ME.025$section=="mid",], 
#             caption = "Hypred Mid Tier, No Association, SD=0 ME=0.025, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD0_ME.05[hypred_noassoc_SD0_ME.05$section=="mid",], 
#             caption = "Hypred Mid Tier, No Association, SD=0 ME=0.05, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD0_ME.075[hypred_noassoc_SD0_ME.075$section=="mid",], 
#             caption = "Hypred Mid Tier, No Association, SD=0 ME=0.075, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD0_ME.1[hypred_noassoc_SD0_ME.1$section=="mid",], 
#             caption = "Hypred Mid Tier, No Association, SD=0 ME=0.1, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD0_ME.125[hypred_noassoc_SD0_ME.125$section=="mid",], 
#             caption = "Hypred Mid Tier, No Association, SD=0 ME=0.125, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD0_ME.15[hypred_noassoc_SD0_ME.15$section=="mid",], 
#             caption = "Hypred Mid Tier, No Association, SD=0 ME=0.15, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD0_ME.175[hypred_noassoc_SD0_ME.175$section=="mid",], 
#             caption = "Hypred Mid Tier, No Association, SD=0 ME=0.175, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD0_ME.2[hypred_noassoc_SD0_ME.2$section=="mid",], 
#             caption = "Hypred Mid Tier, No Association, SD=0 ME=0.2, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD0_ME.225[hypred_noassoc_SD0_ME.225$section=="mid",], 
#             caption = "Hypred Mid Tier, No Association, SD=0 ME=0.225, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD0_ME.25[hypred_noassoc_SD0_ME.25$section=="mid",], 
#             caption = "Hypred Mid Tier, No Association, SD=0 ME=0.25, Estimated betas")
# 
# final_table(dataset=hypred_noassoc_SD0_ME.5[hypred_noassoc_SD0_ME.5$section=="mid",], 
#             caption = "Hypred Mid Tier, No Association, SD=0 ME=0.5, Estimated betas")
# 
# final_table(dataset=hypred_noassoc_SD0_ME.75[hypred_noassoc_SD0_ME.75$section=="mid",], 
#             caption = "Hypred Mid Tier, No Association, SD=0 ME=0.75, Estimated betas")
# 
# 
# 
# # No Association, SD=1, mid tier
# final_table(dataset=hypred_noassoc_SD1_ME0[hypred_noassoc_SD1_ME0$section=="mid",], 
#             caption = "Hypred Mid Tier, No Association, SD=1 ME=0, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD1_ME.025[hypred_noassoc_SD1_ME.025$section=="mid",], 
#             caption = "Hypred Mid Tier, No Association, SD=1 ME=0.025, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD1_ME.05[hypred_noassoc_SD1_ME.05$section=="mid",], 
#             caption = "Hypred Mid Tier, No Association, SD=1 ME=0.05, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD1_ME.075[hypred_noassoc_SD1_ME.075$section=="mid",], 
#             caption = "Hypred Mid Tier, No Association, SD=1 ME=0.075, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD1_ME.1[hypred_noassoc_SD1_ME.1$section=="mid",], 
#             caption = "Hypred Mid Tier, No Association, SD=1 ME=0.1, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD1_ME.125[hypred_noassoc_SD1_ME.125$section=="mid",], 
#             caption = "Hypred Mid Tier, No Association, SD=1 ME=0.125, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD1_ME.15[hypred_noassoc_SD1_ME.15$section=="mid",], 
#             caption = "Hypred Mid Tier, No Association, SD=1 ME=0.15, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD1_ME.175[hypred_noassoc_SD1_ME.175$section=="mid",], 
#             caption = "Hypred Mid Tier, No Association, SD=1 ME=0.175, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD1_ME.2[hypred_noassoc_SD1_ME.2$section=="mid",], 
#             caption = "Hypred Mid Tier, No Association, SD=1 ME=0.2, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD1_ME.225[hypred_noassoc_SD1_ME.225$section=="mid",], 
#             caption = "Hypred Mid Tier, No Association, SD=1 ME=0.225, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD1_ME.25[hypred_noassoc_SD1_ME.25$section=="mid",], 
#             caption = "Hypred Mid Tier, No Association, SD=1 ME=0.25, Estimated betas")
# 
# final_table(dataset=hypred_noassoc_SD1_ME.5[hypred_noassoc_SD1_ME.5$section=="mid",], 
#             caption = "Hypred Mid Tier, No Association, SD=1 ME=0.5, Estimated betas")
# 
# final_table(dataset=hypred_noassoc_SD1_ME.75[hypred_noassoc_SD1_ME.75$section=="mid",], 
#             caption = "Hypred Mid Tier, No Association, SD=1 ME=0.75, Estimated betas")
# 
# 
# ## No Association, SD=0, low tier
# final_table(dataset=hypred_noassoc_SD0_ME0[hypred_noassoc_SD0_ME0$section=="low",], 
#             caption = "Hypred Low Tier, No Association, SD=0 ME=0, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD0_ME.025[hypred_noassoc_SD0_ME.025$section=="low",], 
#             caption = "Hypred Low Tier, No Association, SD=0 ME=0.025, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD0_ME.05[hypred_noassoc_SD0_ME.05$section=="low",], 
#             caption = "Hypred Low Tier, No Association, SD=0 ME=0.05, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD0_ME.075[hypred_noassoc_SD0_ME.075$section=="low",], 
#             caption = "Hypred Low Tier, No Association, SD=0 ME=0.075, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD0_ME.1[hypred_noassoc_SD0_ME.1$section=="low",], 
#             caption = "Hypred Low Tier, No Association, SD=0 ME=0.1, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD0_ME.125[hypred_noassoc_SD0_ME.125$section=="low",], 
#             caption = "Hypred Low Tier, No Association, SD=0 ME=0.125, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD0_ME.15[hypred_noassoc_SD0_ME.15$section=="low",], 
#             caption = "Hypred Low Tier, No Association, SD=0 ME=0.15, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD0_ME.175[hypred_noassoc_SD0_ME.175$section=="low",], 
#             caption = "Hypred Low Tier, No Association, SD=0 ME=0.175, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD0_ME.2[hypred_noassoc_SD0_ME.2$section=="low",], 
#             caption = "Hypred Low Tier, No Association, SD=0 ME=0.2, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD0_ME.225[hypred_noassoc_SD0_ME.225$section=="low",], 
#             caption = "Hypred Low Tier, No Association, SD=0 ME=0.225, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD0_ME.25[hypred_noassoc_SD0_ME.25$section=="low",], 
#             caption = "Hypred Low Tier, No Association, SD=0 ME=0.25, Estimated betas")
# 
# final_table(dataset=hypred_noassoc_SD0_ME.5[hypred_noassoc_SD0_ME.5$section=="low",], 
#             caption = "Hypred Low Tier, No Association, SD=0 ME=0.5, Estimated betas")
# 
# final_table(dataset=hypred_noassoc_SD0_ME.75[hypred_noassoc_SD0_ME.75$section=="low",], 
#             caption = "Hypred Low Tier, No Association, SD=0 ME=0.75, Estimated betas")
# 
# 
# ## No Association, SD=1, low tier
# final_table(dataset=hypred_noassoc_SD1_ME0[hypred_noassoc_SD1_ME0$section=="low",], 
#             caption = "Hypred Low Tier, No Association, SD=1 ME=0, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD1_ME.025[hypred_noassoc_SD1_ME.025$section=="low",], 
#             caption = "Hypred Low Tier, No Association, SD=1 ME=0.025, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD1_ME.05[hypred_noassoc_SD1_ME.05$section=="low",], 
#             caption = "Hypred Low Tier, No Association, SD=1 ME=0.05, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD1_ME.075[hypred_noassoc_SD1_ME.075$section=="low",], 
#             caption = "Hypred Low Tier, No Association, SD=1 ME=0.075, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD1_ME.1[hypred_noassoc_SD1_ME.1$section=="low",], 
#             caption = "Hypred Low Tier, No Association, SD=1 ME=0.1, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD1_ME.125[hypred_noassoc_SD1_ME.125$section=="low",], 
#             caption = "Hypred Low Tier, No Association, SD=1 ME=0.125, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD1_ME.15[hypred_noassoc_SD1_ME.15$section=="low",], 
#             caption = "Hypred Low Tier, No Association, SD=1 ME=0.15, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD1_ME.175[hypred_noassoc_SD1_ME.175$section=="low",], 
#             caption = "Hypred Low Tier, No Association, SD=1 ME=0.175, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD1_ME.2[hypred_noassoc_SD1_ME.2$section=="low",], 
#             caption = "Hypred Low Tier, No Association, SD=1 ME=0.2, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD1_ME.225[hypred_noassoc_SD1_ME.225$section=="low",], 
#             caption = "Hypred Low Tier, No Association, SD=1 ME=0.225, Estimated Betas")
# 
# final_table(dataset=hypred_noassoc_SD1_ME.25[hypred_noassoc_SD1_ME.25$section=="low",], 
#             caption = "Hypred Low Tier, No Association, SD=1 ME=0.25, Estimated betas")
# 
# final_table(dataset=hypred_noassoc_SD1_ME.5[hypred_noassoc_SD1_ME.5$section=="low",], 
#             caption = "Hypred Low Tier, No Association, SD=1 ME=0.5, Estimated betas")
# 
# final_table(dataset=hypred_noassoc_SD1_ME.75[hypred_noassoc_SD1_ME.75$section=="low",], 
#             caption = "Hypred Low Tier, No Association, SD=1 ME=0.75, Estimated betas")
# #####################################################################################################
# 
# 
# ## Misspecified, SD=0, mid tier
# final_table(dataset=hypred_misspec_SD0_ME0[hypred_misspec_SD0_ME0$section=="mid",], 
#             caption = "Hypred Mid Tier, Misspecified, SD=0 ME=0, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD0_ME.025[hypred_misspec_SD0_ME.025$section=="mid",], 
#             caption = "Hypred Mid Tier, Misspecified, SD=0 ME=0.025, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD0_ME.05[hypred_misspec_SD0_ME.05$section=="mid",], 
#             caption = "Hypred Mid Tier, Misspecified, SD=0 ME=0.05, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD0_ME.075[hypred_misspec_SD0_ME.075$section=="mid",], 
#             caption = "Hypred Mid Tier, Misspecified, SD=0 ME=0.075, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD0_ME.1[hypred_misspec_SD0_ME.1$section=="mid",], 
#             caption = "Hypred Mid Tier, Misspecified, SD=0 ME=0.1, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD0_ME.125[hypred_misspec_SD0_ME.125$section=="mid",], 
#             caption = "Hypred Mid Tier, Misspecified, SD=0 ME=0.125, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD0_ME.15[hypred_misspec_SD0_ME.15$section=="mid",], 
#             caption = "Hypred Mid Tier, Misspecified, SD=0 ME=0.15, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD0_ME.175[hypred_misspec_SD0_ME.175$section=="mid",], 
#             caption = "Hypred Mid Tier, Misspecified, SD=0 ME=0.175, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD0_ME.2[hypred_misspec_SD0_ME.2$section=="mid",], 
#             caption = "Hypred Mid Tier, Misspecified, SD=0 ME=0.2, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD0_ME.225[hypred_misspec_SD0_ME.225$section=="mid",], 
#             caption = "Hypred Mid Tier, Misspecified, SD=0 ME=0.225, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD0_ME.25[hypred_misspec_SD0_ME.25$section=="mid",], 
#             caption = "Hypred Mid Tier, Misspecified, SD=0 ME=0.25, Estimated betas")
# 
# final_table(dataset=hypred_misspec_SD0_ME.5[hypred_misspec_SD0_ME.5$section=="mid",], 
#             caption = "Hypred Mid Tier, Misspecified, SD=0 ME=0.5, Estimated betas")
# 
# final_table(dataset=hypred_misspec_SD0_ME.75[hypred_misspec_SD0_ME.75$section=="mid",], 
#             caption = "Hypred Mid Tier, Misspecified, SD=0 ME=0.75, Estimated betas")
# 
# 
# 
# ## Misspecified, SD=0, low tier
# final_table(dataset=hypred_misspec_SD0_ME0[hypred_misspec_SD0_ME0$section=="low",], 
#             caption = "Hypred Low Tier, Misspecified, SD=0 ME=0, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD0_ME.025[hypred_misspec_SD0_ME.025$section=="low",], 
#             caption = "Hypred Low Tier, Misspecified, SD=0 ME=0.025, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD0_ME.05[hypred_misspec_SD0_ME.05$section=="low",], 
#             caption = "Hypred Low Tier, Misspecified, SD=0 ME=0.05, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD0_ME.075[hypred_misspec_SD0_ME.075$section=="low",], 
#             caption = "Hypred Low Tier, Misspecified, SD=0 ME=0.075, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD0_ME.1[hypred_misspec_SD0_ME.1$section=="low",], 
#             caption = "Hypred Low Tier, Misspecified, SD=0 ME=0.1, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD0_ME.125[hypred_misspec_SD0_ME.125$section=="low",], 
#             caption = "Hypred Low Tier, Misspecified, SD=0 ME=0.125, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD0_ME.15[hypred_misspec_SD0_ME.15$section=="low",], 
#             caption = "Hypred Low Tier, Misspecified, SD=0 ME=0.15, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD0_ME.175[hypred_misspec_SD0_ME.175$section=="low",], 
#             caption = "Hypred Low Tier, Misspecified, SD=0 ME=0.175, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD0_ME.2[hypred_misspec_SD0_ME.2$section=="low",], 
#             caption = "Hypred Low Tier, Misspecified, SD=0 ME=0.2, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD0_ME.225[hypred_misspec_SD0_ME.225$section=="low",], 
#             caption = "Hypred Low Tier, Misspecified, SD=0 ME=0.225, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD0_ME.25[hypred_misspec_SD0_ME.25$section=="low",], 
#             caption = "Hypred Low Tier, Misspecified, SD=0 ME=0.25, Estimated betas")
# 
# final_table(dataset=hypred_misspec_SD0_ME.5[hypred_misspec_SD0_ME.5$section=="low",], 
#             caption = "Hypred Low Tier, Misspecified, SD=0 ME=0.5, Estimated betas")
# 
# final_table(dataset=hypred_misspec_SD0_ME.75[hypred_misspec_SD0_ME.75$section=="low",], 
#             caption = "Hypred Low Tier, Misspecified, SD=0 ME=0.75, Estimated betas")
# 
# 
# 
# 
# # Misspecified, SD=1, mid tier
# final_table(dataset=hypred_misspec_SD1_ME0[hypred_misspec_SD1_ME0$section=="mid",], 
#             caption = "Hypred Mid Tier, Misspecified, SD=1 ME=0, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD1_ME.025[hypred_misspec_SD1_ME.025$section=="mid",], 
#             caption = "Hypred Mid Tier, Misspecified, SD=1 ME=0.025, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD1_ME.05[hypred_misspec_SD1_ME.05$section=="mid",], 
#             caption = "Hypred Mid Tier, Misspecified, SD=1 ME=0.05, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD1_ME.075[hypred_misspec_SD1_ME.075$section=="mid",], 
#             caption = "Hypred Mid Tier, Misspecified, SD=1 ME=0.075, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD1_ME.1[hypred_misspec_SD1_ME.1$section=="mid",], 
#             caption = "Hypred Mid Tier, Misspecified, SD=1 ME=0.1, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD1_ME.125[hypred_misspec_SD1_ME.125$section=="mid",], 
#             caption = "Hypred Mid Tier, Misspecified, SD=1 ME=0.125, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD1_ME.15[hypred_misspec_SD1_ME.15$section=="mid",], 
#             caption = "Hypred Mid Tier, Misspecified, SD=1 ME=0.15, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD1_ME.175[hypred_misspec_SD1_ME.175$section=="mid",], 
#             caption = "Hypred Mid Tier, Misspecified, SD=1 ME=0.175, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD1_ME.2[hypred_misspec_SD1_ME.2$section=="mid",], 
#             caption = "Hypred Mid Tier, Misspecified, SD=1 ME=0.2, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD1_ME.225[hypred_misspec_SD1_ME.225$section=="mid",], 
#             caption = "Hypred Mid Tier, Misspecified, SD=1 ME=0.225, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD1_ME.25[hypred_misspec_SD1_ME.25$section=="mid",], 
#             caption = "Hypred Mid Tier, Misspecified, SD=1 ME=0.25, Estimated betas")
# 
# final_table(dataset=hypred_misspec_SD1_ME.5[hypred_misspec_SD1_ME.5$section=="mid",], 
#             caption = "Hypred Mid Tier, Misspecified, SD=1 ME=0.5, Estimated betas")
# 
# final_table(dataset=hypred_misspec_SD1_ME.75[hypred_misspec_SD1_ME.75$section=="mid",], 
#             caption = "Hypred Mid Tier, Misspecified, SD=1 ME=0.75, Estimated betas")
# 
# 
# 
# ## Misspecified, SD=1, low tier
# final_table(dataset=hypred_misspec_SD1_ME0[hypred_misspec_SD1_ME0$section=="low",], 
#             caption = "Hypred Low Tier, Misspecified, SD=1 ME=0, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD1_ME.025[hypred_misspec_SD1_ME.025$section=="low",], 
#             caption = "Hypred Low Tier, Misspecified, SD=1 ME=0.025, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD1_ME.05[hypred_misspec_SD1_ME.05$section=="low",], 
#             caption = "Hypred Low Tier, Misspecified, SD=1 ME=0.05, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD1_ME.075[hypred_misspec_SD1_ME.075$section=="low",], 
#             caption = "Hypred Low Tier, Misspecified, SD=1 ME=0.075, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD1_ME.1[hypred_misspec_SD1_ME.1$section=="low",], 
#             caption = "Hypred Low Tier, Misspecified, SD=1 ME=0.1, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD1_ME.125[hypred_misspec_SD1_ME.125$section=="low",], 
#             caption = "Hypred Low Tier, Misspecified, SD=1 ME=0.125, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD1_ME.15[hypred_misspec_SD1_ME.15$section=="low",], 
#             caption = "Hypred Low Tier, Misspecified, SD=1 ME=0.15, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD1_ME.175[hypred_misspec_SD1_ME.175$section=="low",], 
#             caption = "Hypred Low Tier, Misspecified, SD=1 ME=0.175, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD1_ME.2[hypred_misspec_SD1_ME.2$section=="low",], 
#             caption = "Hypred Low Tier, Misspecified, SD=1 ME=0.2, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD1_ME.225[hypred_misspec_SD1_ME.225$section=="low",], 
#             caption = "Hypred Low Tier, Misspecified, SD=1 ME=0.225, Estimated Betas")
# 
# final_table(dataset=hypred_misspec_SD1_ME.25[hypred_misspec_SD1_ME.25$section=="low",], 
#             caption = "Hypred Low Tier, Misspecified, SD=1 ME=0.25, Estimated betas")
# 
# final_table(dataset=hypred_misspec_SD1_ME.5[hypred_misspec_SD1_ME.5$section=="low",], 
#             caption = "Hypred Low Tier, Misspecified, SD=1 ME=0.5, Estimated betas")
# 
# final_table(dataset=hypred_misspec_SD1_ME.75[hypred_misspec_SD1_ME.75$section=="low",], 
#             caption = "Hypred Low Tier, Misspecified, SD=1 ME=0.75, Estimated betas")
# 
# 
# 
# ####################################  Hypred Uganda results #####################################
# 
# ## Hypred mid tier
# 
# final_table(dataset=hypred_uganda_ME0[hypred_uganda_ME0$section=="mid",], 
#             caption="Uganda Results - Hypred Mid Tier: ME=0")
# 
# final_table(dataset=hypred_uganda_ME.025[hypred_uganda_ME.025$section=="mid",], 
#             caption="Uganda Results - Hypred Mid Tier: ME=0.025")
# 
# final_table(dataset=hypred_uganda_ME.05[hypred_uganda_ME.05$section=="mid",], 
#             caption="Uganda Results - Hypred Mid Tier: ME=0.05")
# 
# final_table(dataset=hypred_uganda_ME.075[hypred_uganda_ME.075$section=="mid",], 
#             caption="Uganda Results - Hypred Mid Tier: ME=0.075")
# 
# final_table(dataset=hypred_uganda_ME.1[hypred_uganda_ME.1$section=="mid",], 
#             caption="Uganda Results - Hypred Mid Tier: ME=0.1")
# 
# final_table(dataset=hypred_uganda_ME.125[hypred_uganda_ME.125$section=="mid",], 
#             caption="Uganda Results - Hypred Mid Tier: ME=0.125")
# 
# final_table(dataset=hypred_uganda_ME.15[hypred_uganda_ME.15$section=="mid",], 
#             caption="Uganda Results - Hypred Mid Tier: ME=0.15")
# 
# final_table(dataset=hypred_uganda_ME.175[hypred_uganda_ME.175$section=="mid",], 
#             caption="Uganda Results - Hypred Mid Tier: ME=0.175")
# 
# final_table(dataset=hypred_uganda_ME.2[hypred_uganda_ME.2$section=="mid",], 
#             caption="Uganda Results - Hypred Mid Tier: ME=0.2")
# 
# final_table(dataset=hypred_uganda_ME.225[hypred_uganda_ME.225$section=="mid",], 
#             caption="Uganda Results - Hypred Mid Tier: ME=0.225")
# 
# final_table(dataset=hypred_uganda_ME.25[hypred_uganda_ME.25$section=="mid",], 
#             caption="Uganda Results - Hypred Mid Tier: ME=0.25")
# 
# final_table(dataset=hypred_uganda_ME.5[hypred_uganda_ME.5$section=="mid",], 
#             caption="Uganda Results - Hypred Mid Tier: ME=0.5")
# 
# final_table(dataset=hypred_uganda_ME.75[hypred_uganda_ME.75$section=="mid",], 
#             caption="Uganda Results - Hypred Mid Tier: ME=0.75")
# 
# 
# ## Low tier
# 
# final_table(dataset=hypred_uganda_ME0[hypred_uganda_ME0$section=="low",], 
#             caption="Uganda Results - Hypred Low Tier: ME=0")
# 
# final_table(dataset=hypred_uganda_ME.025[hypred_uganda_ME.025$section=="low",], 
#             caption="Uganda Results - Hypred Low Tier: ME=0.025")
# 
# final_table(dataset=hypred_uganda_ME.05[hypred_uganda_ME.05$section=="low",], 
#             caption="Uganda Results - Hypred Low Tier: ME=0.05")
# 
# final_table(dataset=hypred_uganda_ME.075[hypred_uganda_ME.075$section=="low",], 
#             caption="Uganda Results - Hypred Low Tier: ME=0.075")
# 
# final_table(dataset=hypred_uganda_ME.1[hypred_uganda_ME.1$section=="low",], 
#             caption="Uganda Results - Hypred Low Tier: ME=0.1")
# 
# final_table(dataset=hypred_uganda_ME.125[hypred_uganda_ME.125$section=="low",], 
#             caption="Uganda Results - Hypred Low Tier: ME=0.125")
# 
# final_table(dataset=hypred_uganda_ME.15[hypred_uganda_ME.15$section=="low",], 
#             caption="Uganda Results - Hypred Low Tier: ME=0.15")
# 
# final_table(dataset=hypred_uganda_ME.175[hypred_uganda_ME.175$section=="low",], 
#             caption="Uganda Results - Hypred Low Tier: ME=0.175")
# 
# final_table(dataset=hypred_uganda_ME.2[hypred_uganda_ME.2$section=="low",], 
#             caption="Uganda Results - Hypred Low Tier: ME=0.2")
# 
# final_table(dataset=hypred_uganda_ME.225[hypred_uganda_ME.225$section=="low",], 
#             caption="Uganda Results - Hypred Low Tier: ME=0.225")
# 
# final_table(dataset=hypred_uganda_ME.25[hypred_uganda_ME.25$section=="low",], 
#             caption="Uganda Results - Hypred Low Tier: ME=0.25")
# 
# final_table(dataset=hypred_uganda_ME.5[hypred_uganda_ME.5$section=="low",], 
#             caption="Uganda Results - Hypred Low Tier: ME=0.5")
# 
# final_table(dataset=hypred_uganda_ME.75[hypred_uganda_ME.75$section=="low",], 
#             caption="Uganda Results - Hypred Low Tier: ME=0.75")
# 


# combines results for the 3 tier risk groups for the hypred method; user must define which method
# to use for the mid and low tier. High tier is assumed individual testing
# function only outputs means, no CIs
# this is what is reported in Tables 1-3 in the paper
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
    mid$edet <- mid$lrsoe_esens*mid$lrsoe_eprevfail
    mid$t <- mid$lrsoe_t
    mid$rds <- mid$lrsoe_rds
  }
  else if (mid_method=="mincov"){
    mid$eprevfail <- mid$mincov_eprevfail
    mid$edet <- mid$mincov_esens*mid$mincov_eprevfail
    mid$t <- mid$mincov_t
    mid$rds <- mid$mincov_rds
  }
  else if (mid_method=="mini"){
    mid$eprevfail <- mid$mini_eprevfail
    mid$edet <- mid$mini_esens*mid$mini_eprevfail
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
    low$edet <- low$lrsoe_esens*low$lrsoe_eprevfail
    low$t <- low$lrsoe_t
    low$rds <- low$lrsoe_rds
  }
  else if (low_method=="mincov"){
    low$eprevfail <- low$mincov_eprevfail
    low$edet <- low$mincov_esens*low$mincov_eprevfail
    low$t <- low$mincov_t
    low$rds <- low$mincov_rds
  }
  else if (low_method=="mini"){
    low$eprevfail <- low$mini_eprevfail
    low$edet <- low$mini_esens*low$mini_eprevfail
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
  rds_num <- (top$all_t/20) + sum(mid$rds) + sum(low$rds)
  rds <- round(rds_num/rds_denom, digits=1)
  
  return(c(sens, eff, rds))
}



# AGAIG, SD=0
comb_hypred(hypred_AGAIG_SD0_ME0, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_AGAIG_SD0_ME.025, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_AGAIG_SD0_ME.05, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_AGAIG_SD0_ME.075, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_AGAIG_SD0_ME.1, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_AGAIG_SD0_ME.125, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_AGAIG_SD0_ME.15, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_AGAIG_SD0_ME.175, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_AGAIG_SD0_ME.2, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_AGAIG_SD0_ME.225, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_AGAIG_SD0_ME.25, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_AGAIG_SD0_ME.5, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_AGAIG_SD0_ME.75, mid_method="mincov", low_method="mincov")

# AGAIG, SD=1
comb_hypred(hypred_AGAIG_SD1_ME0, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_AGAIG_SD1_ME.025, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_AGAIG_SD1_ME.05, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_AGAIG_SD1_ME.075, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_AGAIG_SD1_ME.1, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_AGAIG_SD1_ME.125, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_AGAIG_SD1_ME.15, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_AGAIG_SD1_ME.175, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_AGAIG_SD1_ME.2, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_AGAIG_SD1_ME.225, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_AGAIG_SD1_ME.25, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_AGAIG_SD1_ME.5, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_AGAIG_SD1_ME.75, mid_method="mincov", low_method="mincov")

# reverse, SD=0
comb_hypred(hypred_reverse_SD0_ME0, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_reverse_SD0_ME.025, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_reverse_SD0_ME.05, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_reverse_SD0_ME.075, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_reverse_SD0_ME.1, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_reverse_SD0_ME.125, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_reverse_SD0_ME.15, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_reverse_SD0_ME.175, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_reverse_SD0_ME.2, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_reverse_SD0_ME.225, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_reverse_SD0_ME.25, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_reverse_SD0_ME.5, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_reverse_SD0_ME.75, mid_method="mincov", low_method="mincov")

# reverse, SD=1
comb_hypred(hypred_reverse_SD1_ME0, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_reverse_SD1_ME.025, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_reverse_SD1_ME.05, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_reverse_SD1_ME.075, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_reverse_SD1_ME.1, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_reverse_SD1_ME.125, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_reverse_SD1_ME.15, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_reverse_SD1_ME.175, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_reverse_SD1_ME.2, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_reverse_SD1_ME.225, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_reverse_SD1_ME.25, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_reverse_SD1_ME.5, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_reverse_SD1_ME.75, mid_method="mincov", low_method="mincov")

# uganda real data
comb_hypred(hypred_uganda_ME0, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_uganda_ME.025, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_uganda_ME.05, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_uganda_ME.075, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_uganda_ME.1, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_uganda_ME.125, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_uganda_ME.15, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_uganda_ME.175, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_uganda_ME.2, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_uganda_ME.225, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_uganda_ME.25, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_uganda_ME.5, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_uganda_ME.75, mid_method="mincov", low_method="mincov")


## No assoc Sd=0
comb_hypred(hypred_noassoc_SD0_ME0, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_noassoc_SD0_ME.025, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_noassoc_SD0_ME.05, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_noassoc_SD0_ME.075, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_noassoc_SD0_ME.1, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_noassoc_SD0_ME.125, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_noassoc_SD0_ME.15, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_noassoc_SD0_ME.175, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_noassoc_SD0_ME.2, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_noassoc_SD0_ME.225, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_noassoc_SD0_ME.25, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_noassoc_SD0_ME.5, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_noassoc_SD0_ME.75, mid_method="mincov", low_method="mincov")


## No assoc Sd=1
comb_hypred(hypred_noassoc_SD1_ME0, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_noassoc_SD1_ME.025, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_noassoc_SD1_ME.05, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_noassoc_SD1_ME.075, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_noassoc_SD1_ME.1, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_noassoc_SD1_ME.125, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_noassoc_SD1_ME.15, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_noassoc_SD1_ME.175, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_noassoc_SD1_ME.2, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_noassoc_SD1_ME.225, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_noassoc_SD1_ME.25, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_noassoc_SD1_ME.5, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_noassoc_SD1_ME.75, mid_method="mincov", low_method="mincov")


## Misspec Sd=0
comb_hypred(hypred_misspec_SD0_ME0, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_misspec_SD0_ME.025, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_misspec_SD0_ME.05, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_misspec_SD0_ME.075, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_misspec_SD0_ME.1, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_misspec_SD0_ME.125, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_misspec_SD0_ME.15, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_misspec_SD0_ME.175, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_misspec_SD0_ME.2, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_misspec_SD0_ME.225, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_misspec_SD0_ME.25, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_misspec_SD0_ME.5, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_misspec_SD0_ME.75, mid_method="mincov", low_method="mincov")


## Misspec Sd=1
comb_hypred(hypred_misspec_SD1_ME0, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_misspec_SD1_ME.025, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_misspec_SD1_ME.05, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_misspec_SD1_ME.075, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_misspec_SD1_ME.1, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_misspec_SD1_ME.125, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_misspec_SD1_ME.15, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_misspec_SD1_ME.175, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_misspec_SD1_ME.2, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_misspec_SD1_ME.225, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_misspec_SD1_ME.25, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_misspec_SD1_ME.5, mid_method="mincov", low_method="mincov")
comb_hypred(hypred_misspec_SD1_ME.75, mid_method="mincov", low_method="mincov")


