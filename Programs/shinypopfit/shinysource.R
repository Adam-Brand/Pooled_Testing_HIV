#==============================================================================
# FILENAME: shinysource.R
# PROJECT: 	Pooled testing in HIV
# PURPOSE: this is the source file for the shiny app called shinypopfit
#          
#          
# AUTHOR: Adam Brand


# INPUT datasets: None

# OUTPUT: None, this is only a source file


# R VERSION: 3.6.1
#==============================================================================
#Notes: as a source file, it needs to be in the same folder as the server.R file, but does not need to be run





# =============================================================================



# Finally, also in a separate R file, call this one shinysource.R

Under500VL <- function(n){
  a <- round(.85*n)
  b <- round(.05*n)
  c <- n - a - b
  y <- c(runif(a,0,50), runif(b,50,100),runif(c,100,500))
  return(y)
}

Over500VL <- function(n){
  a <- round(.93*n)
  b <- n - a
  y <- c(replicate(a, 10^(2.7 + rgamma(1,shape = 1.6, scale = 0.5))),
         replicate(b, 10^(2.7 + rgamma(1,shape = 3.18, scale = 0.5))))
  return(y)
}

prevover500 <- function(prevfail, cutoff){
  if (50 <= cutoff & cutoff < 100){
    prevover500 = (prevfail - .10 - (.05*(100-cutoff)/50))/(1 - .1 - (.05*(100-cutoff)/50))}
  else if (100 <= cutoff & cutoff < 500){
    prevover500 = (prevfail - (.10*(500-cutoff)/400))/(1 - (.10*(500-cutoff)/400))}
  else if (cutoff == 500){
    prevover500 = prevfail}
  else if (cutoff > 500){
    prevover500 = prevfail/((.93*(1 - pgamma(log10(cutoff)-2.7, shape = 1.6, scale = 0.5)) +
                               .07*(1 - pgamma(log10(cutoff)-2.7, shape = 3.18, scale = 0.5))))}
  else {prevover500 = "error"}
  return(prevover500)
}

createpop <- function(n, prevover500, cutoff, prevfail=TRUE){
  if (prevfail==FALSE){
    prevover500 = prevover500}
  else {prevover500 <- prevover500(prevover500, cutoff)}
  a <- round(n*prevover500)
  b <- n - a
  over500 <- Over500VL(a)
  under500 <- Under500VL(b)
  pop <- c(over500, under500)
  return(pop)
}

# Simulates a dataframe with adherence, PFS and viral load for n patients
# based on the covariate model, described in detai later
VL.data.adh <- function(n, shape1, shape2, b0, b1, b2, b3, pf.prev, sd){
  pf <- rbinom(n=n, size=1, prob=pf.prev)
  adhere <- (100 - 100*rbeta(n=n, shape1 = shape1, shape2=shape2))
  var <- rnorm(n=n, mean=0, sd = sd)
  log.VL <- b0 + b1*adhere + b2*pf + b3*adhere*pf + var
  VL <- 10^log.VL
  data <- cbind(VL, adhere, pf)
  return(data.frame(data))
}

# Plots the histograms
plot.pop <- function(n, prev.over.cutoff, cutoff, shape1, shape2, pf.prev, b0, b1, b2, b3, sd){
  data <- createpop(n, prev.over.cutoff, cutoff, prevfail=TRUE)
  Logdata <- log10(data)
  modeldata <- VL.data.adh(n=n, shape1=shape1, shape2=shape2, b0=b0, b1=b1, b2=b2, b3=b3, pf.prev=pf.prev, sd=sd)
  modelVL <- log10(modeldata$VL)
  hist(Logdata, freq=FALSE, col=rgb(1,0,0,0.5), xlim = c(0,10), breaks=21, 
       xlab="Log10(VL); Red is Pop VL, Blue is Model Fit", main="Log10 VL")
  hist(modelVL, col=rgb(0,0,1,0.5), freq=FALSE, add=TRUE, breaks=21)
  box()
}
