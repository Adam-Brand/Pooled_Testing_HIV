#==============================================================================
# FILENAME: results_anal_uganda.R
# PROJECT: 	Pooled testing in HIV
# PURPOSE: compile the method evaluation results using the real uganda data into performance summaries
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
library(here)

# reading in results from methods except HyPred
uganda_ME0 <- readRDS("Results/UgandaResults/Uganda_ME0.rds")
uganda_ME.025 <- readRDS("Results/UgandaResults/Uganda_ME.025.rds")
uganda_ME.05 <- readRDS("Results/UgandaResults/Uganda_ME.05.rds")
uganda_ME.075 <- readRDS("Results/UgandaResults/Uganda_ME.075.rds")
uganda_ME.1 <- readRDS("Results/UgandaResults/Uganda_ME.1.rds")
uganda_ME.125 <- readRDS("Results/UgandaResults/Uganda_ME.125.rds")
uganda_ME.15 <- readRDS("Results/UgandaResults/Uganda_ME.15.rds")
uganda_ME.175 <- readRDS("Results/UgandaResults/Uganda_ME.175.rds")
uganda_ME.2 <- readRDS("Results/UgandaResults/Uganda_ME.2.rds")
uganda_ME.225 <- readRDS("Results/UgandaResults/Uganda_ME.225.rds")
uganda_ME.25 <- readRDS("Results/UgandaResults/Uganda_ME.25.rds")
uganda_ME.5 <- readRDS("Results/UgandaResults/Uganda_ME.5.rds")
uganda_ME.75 <- readRDS("Results/UgandaResults/Uganda_ME.75.rds")




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


###################### results for methods except HyPred #####################################

final_table(dataset=uganda_ME0, 
            caption="Uganda Results: ME=0")

final_table(dataset=uganda_ME.025, 
            caption="Uganda Results: ME=0.025")

final_table(dataset=uganda_ME.05, 
            caption="Uganda Results: ME=0.05")

final_table(dataset=uganda_ME.075, 
            caption="Uganda Results: ME=0.075")

final_table(dataset=uganda_ME.1, 
            caption="Uganda Results: ME=0.1")

final_table(dataset=uganda_ME.125, 
            caption="Uganda Results: ME=0.125")

final_table(dataset=uganda_ME.15, 
            caption="Uganda Results: ME=0.15")

final_table(dataset=uganda_ME.175, 
            caption="Uganda Results: ME=0.175")

final_table(dataset=uganda_ME.2, 
            caption="Uganda Results: ME=0.2")

final_table(dataset=uganda_ME.225, 
            caption="Uganda Results: ME=0.225")

final_table(dataset=uganda_ME.25, 
            caption="Uganda Results: ME=0.25")

final_table(dataset=uganda_ME.5, 
            caption="Uganda Results: ME=0.5")

final_table(dataset=uganda_ME.75, 
            caption="Uganda Results: ME=0.75")

