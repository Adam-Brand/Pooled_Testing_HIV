##### RA Pooling Source ######


###### simulating VL levels based on
###### distributions described in "Pooled Nucleic Acid Testing to Identify Antiretroviral Treatment Failure During HIV Infection"

# This function simulated viral load values under 500 copies/mL
Under500VL <- function(n){
  a <- round(.85*n)
  b <- round(.05*n)
  c <- n - a - b
  y <- c(runif(a,0,50), runif(b,50,100),runif(c,100,500))
  return(y)
}

# responses for over 500 copies/mL
#### mixed gamma distr.

Over500VL <- function(n){
  a <- round(.93*n)
  b <- n - a
  y <- c(replicate(a, 10^(2.7 + rgamma(1,shape = 1.6, scale = 0.5))),
         replicate(b, 10^(2.7 + rgamma(1,shape = 3.18, scale = 0.5))))
  return(y)
}


#### Function which sets the prev of failure based on the prevover500 & the cutoff
#### minimum cutoff is 50
prevfail <- function(prevover500, cutoff){
  if (50 <= cutoff & cutoff < 100){
    prevfail = prevover500 + .10 + .05((100-cutoff)/50)}
  else if (100 <= cutoff & cutoff < 500){
    prevfail = prevover500 + .10*((500-cutoff)/400)}
  else if (cutoff == 500){
    prevfail = prevover500}
  else if (cutoff > 500){
    prevfail = prevover500*((.93*(1 - pgamma(log10(cutoff)-2.7, shape = 1.6, scale = 0.5)) +
                               .07*(1 - pgamma(log10(cutoff)-2.7, shape = 3.18, scale = 0.5))))}
  else {prevfail = "error"}
  return(prevfail)
}

##### This function sets the prevover500 based on the prev of failure
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

##### Function to create population of HIV patients based on prev over 500
##### and the cutoff

createpop <- function(n, prevover500, cutoff, prevfail=FALSE){
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


#### function which makes one testable matrix from the population
one.matrix <- function(data, matsize, SE){
  m <- matrix(nrow=matsize, ncol=matsize)
  e <- matrix(nrow=matsize, ncol=matsize)
  for (i in 1:matsize){
    for(j in 1:matsize){
      index <- sample.int(length(data$VL),1,replace=FALSE)
      m[i,j] = data[index, "VL"]
      e[i,j] = data[index,"VLobs"]
      data <- data[-c(index),]
    }}
  x <- data.frame(cbind(m,e))
  rowst <- replicate(matsize, 0)
  for (i in 1:matsize){rowst[i] = sum(x[i,1:matsize]/matsize)}
  rows <- rowst*(10^(rnorm(1, sd=SE)))
  colst <- replicate(matsize, 0)
  for (i in 1:matsize){colst[i] = sum(x[1:matsize,i]/matsize)}
  cols <- colst*(10^(rnorm(1, sd=SE)))
  return(cbind(x,rows,cols))
}


#### this zeros out all rows and columns for the true values 
##   if the row and column totals are less than cutoff/n
zerocolrow <- function(y, matsize, cutoff, lowlimit){
  
  ###checks each row to see if its below the cutoff/matsize
  for(i in 1:matsize){
    if (y$rows[i] < cutoff/matsize){
      ### this chunk zeros out the row total 
      num <- y$rows[i]
      y$rows[i] = 0
      ### alters col totals based on assigning equal VL value to each box in row which
      ### has not already been tested
      for (k in 1:matsize){
        if (y[i, (matsize+k)] == 0){}
        else{
          y$cols[k] <- y$cols[k] - num/matsize
          y[i,k] = 0
          y[[i, (matsize+k)]] = 0
        }
      }
    }
    else{}}
  for(i in 1:matsize){
    if (y$cols[i] < cutoff/matsize){
      num <- y$cols[i]
      y$cols[i] = 0
      for (k in 1:matsize){
        if (y[k, (matsize+i)] == 0){}
        else{
          y$rows[k] <- y$rows[k] - num/matsize
          y[k,i] = 0
          y[[k, (matsize+i)]] = 0
        }
      }
    }
    else{}}
  return(y)
}


#### this function checks to see if all the elements in a row or column are all zero.  If so, it sets that column or row value to zero
checkrowcol <- function(y, n){
  for (i in 1:n){
    if(sum(y[i,1:n]) == 0){y$rows[i] = 0}
    else{}
    if(sum(y[1:n,i]) == 0){y$cols[i] = 0}
    else{}
  }
  return(y)
}


#### This function checks to see if the listed index is in a matrix
#### returns True or False
### ind needs to be input as a vector e.g. c(i,j)
check.indices <- function(mat, ind){
  any(apply(mat, 1, function(x, want) isTRUE(all.equal(x, want)), ind))
}

##### changes data frame into vector
df2vec <- function(dataframe, left, right){
  z <- NULL
  y <- NULL
  for (i in (left):(right)){
    y <- as.vector(dataframe[,i])
    z <- c(z,y)}
  return(z)}

#### calculates the prevalences of the original matrix
#### returns true prev over 500, true prevfail, error prev over 500, error prevfail
prev <- function(x, matsize, cutoff){
  probs = c(.9,.95,.99,1)
  xt <- x[,1:matsize]
  xe <- x[,(matsize+1):(matsize*2)]
  xtv <- df2vec(x, 1, matsize)
  xev <- df2vec(x, matsize+1, matsize*2)
  tquant <- quantile(xtv, probs=probs)
  equant <- quantile(xev, probs=probs)
  tprev500 <- length(xt[xt>=500])/(matsize^2)
  eprev500 <- length(xe[xe>=500])/(matsize^2)
  tprevfail <- length(xt[xt>=cutoff])
  eprevfail <- length(xe[xe>=cutoff])
 
  
  result <- data.frame(cbind(tprevfail, eprevfail))
  
  return(result)
}


#### creates data frame of all possible sums of row and column totals
### indexed by the ith row and jth column
### if an element is zero, it stays zero
### x is a matrix
sums <- function(x, matsize){
  s <- matrix(nrow = matsize, ncol = matsize)
  sums <- data.frame(s)
  for (i in 1:matsize){
    for (j in 1:matsize){
      if (x[i,j]==0){
        sums[i,j]=0}
      else{sums[i,j] <- x$rows[i] + x$cols[j]}}}
  return(sums)
}

##### function which compares a reduced matrix to the original 
#### and calculates the sensitivity
#### this returns prop of true failues we detect, prop of error fails we detect
sens <- function(y, z, matsize, cutoff){
  a <- 0
  b <- 0
  n <- matsize +1
  m <- matsize*2
  ## this loop checks for each true failure, how many do we conclude are failures
  for (i in 1:matsize){
    for (j in 1:matsize){
      if (z[i,j]>cutoff){
        if (y[i,j]==0 & y[[i,(matsize + j)]]>=cutoff){
          a <- a+1}
        else{a <- a}}
      else{}}}
  ## this loops checks for each error failure, how many do we conclude are failures
  for (i in 1:matsize){
    for (j in 1:matsize){
      if (z[[i,(matsize + j)]]>cutoff){
        if (y[i,j]==0 & y[[i,(matsize + j)]]>=cutoff){
          b <- b+1}
        else{b <- b}}
      else{}}}
  truth <- z[,1:matsize]
  error <- z[,n:m]
  true.sens <- a/length(truth[truth>=cutoff])
  error.sens <- b/length(error[error>=cutoff])
  result <- data.frame(cbind(true.sens, error.sens))
  return (result)
}

### this function sets the accuracy bar for true sensitivity; adjusts
### for measurement error (with measurement error, true sensitivity may not be 100%)
accuracy <- function(n, prevover500, cutoff, SE, prevfail=FALSE){
  pop <- createpop(n, prevover500, cutoff, prevfail=prevfail)
  poperror <- pop*(10^(rnorm(1, sd=SE)))
  data <- data.frame(cbind(pop,poperror))
  fail <- length(pop[pop>cutoff])
  x <- 0
  for (i in 1:n){
    if (pop[i] > cutoff & poperror[i] > cutoff){x <- x+1}}
  accuracy <- x/fail
  return(accuracy=accuracy)}


#################### Star of individual methods ##############################################




#### this function takes a one.matrix and adjusts
### the row and column values so that the sum of rows = sum of cols
soe.matrix <- function(x, matsize){
  tot <- (sum(x$rows) + sum(x$cols))/2
  rowdif <- (tot - sum(x$rows))/sum(x$rows)
  coldif <- (tot - sum(x$cols))/sum(x$cols)
  x$rows <- x$rows*(1+rowdif)
  x$cols <- x$cols*(1+coldif)
  return(x) 
}


##### Introduction to Covariates ####

##### this function takes in a population of any size and a specified correlation
##### it then generates x-values following the model x = log10(y) + error
##### it returns the std dev of the error term which most closely estimates the 
##### desired correlation
true.corr <- function(pop, corr, b0, b1){
  popx <- NULL
  z <- seq(from = .01, to = 2.0, by = .01)
  result <- NULL
  diff <- NULL
  for(j in 1:length(z)){
    for (i in 1:length(pop)){
      popx[i] <- (log10(pop[i]) - b0)/b1 + rnorm(1, sd=z[j])}
    result[j] <- cor(log10(pop), popx)}
  diff <- abs(result - corr)
  return(which.min(diff))
}

#### replicates true.corr to get a better estimate
true.corr.est <- function(pop, corr, b0, b1, reps){
  y <- NULL
  y <- replicate(reps, true.corr(pop, corr, b0, b1))
  return(round(mean(y))/100)
}

diag.mat.fill <- function(VLdata, matsize){
  m <- matrix(nrow=matsize, ncol=matsize)
  t <- 0
  g <- 0
  for (j in 1:(2*matsize-1)){
    if (isEven(j)){
      k <- j-g
      for(i in 1:(matsize - g - 1)){
        m[k,i] = VLdata[1]
        VLdata <- tail(VLdata, -1)
        k <- k+1}
      g <- g+1}
    else{
      k <- j-t
      for(i in 1:(matsize - t)){
        m[i,k] = VLdata[1]
        VLdata <- tail(VLdata, -1)
        k <- k+1}
      t <- t+1}}
  return(m)}

rnd.mat.fill <- function(data, matsize, SE, Uganda=FALSE){
  m <- matrix(nrow=matsize, ncol=matsize)
#  e <- matrix(nrow=matsize, ncol=matsize)
  g <- matrix(nrow=matsize, ncol=matsize)
  
  for (i in 1:matsize){
    for(j in 1:matsize){
      index <- sample.int(length(data$VL),1,replace=FALSE)
      m[i,j] = data[index,"VL"]
     # e[i,j] = data[index,"VLobs"]
      g[i,j] = data[index,"predVL"]
      data <- data[-index,]
    }}
  e <- data.frame(m)
if (Uganda==FALSE){
  for (i in 1:matsize){
    for(j in 1:matsize){
      e[i,j] = e[i,j]*(10^(rnorm(1, sd=SE)))}}
}
  x <- data.frame(cbind(m,e,g))
  return(x)
}



corner.mat.fill <- function(data, matsize, SE){
  m <- matrix(nrow=matsize, ncol=matsize)
#  e <- matrix(nrow=matsize, ncol=matsize)
  g <- matrix(nrow=matsize, ncol=matsize)
  for (z in 1:matsize){
    i <- 1
    for (j in z:1){
      m[i,j] = data[1,"VL"]
    #  e[i,j] = data[1,"VLobs"]
      g[i,j] = data[1,"predVL"]
      data <- data[-1,]
      i <- i+1}}
  
  for (z in 2:matsize){
    i <- z
    for (j in matsize:z){
      m[i,j] = data[1,"VL"]
     # e[i,j] = data[1,"VLobs"]
      g[i,j] = data[1,"predVL"]
      data <- data[-1,]
      i <- i+1}}
  
  e <- data.frame(m)
  for (i in 1:matsize){
    for(j in 1:matsize){
      e[i,j] = e[i,j]*(10^(rnorm(1, sd=SE)))}}
  x <- data.frame(cbind(m,e,g))
  return(x)
}

row.mat.fill <- function(VLdata, matsize){
  m <- matrix(nrow=matsize, ncol=matsize)
  for (i in 1:matsize){
    for (j in 1:matsize){
      m[i,j] = VLdata[1]
      VLdata <- tail(VLdata, -1)}}
  return(m)
}

border.mat.fill <- function(VLdata, matsize){
  m <- matrix(nrow=matsize, ncol=matsize)
  z <- matsize
  t <- 0
  for(z in matsize:(matsize-1)){
    t <- t + 1
    for (j in t:z){
      m[t,j] = VLdata[1]
      VLdata <- tail(VLdata, -1)}
    for (i in (t+1):(z-1)){
      m[i,z] = VLdata[1]
      VLdata <- tail(VLdata, -1)}
    for (j in z:t){
      m[z,j] = VLdata[1]
      VLdata <- tail(VLdata, -1)}
    for (i in (z-1):(t+1)){
      m[i,t] = VLdata[1]
      VLdata <- tail(VLdata, -1)}}
  for (i in 3:(matsize-2)){
    for (j in 3:(matsize-2)){
      m[i,j] = sample(VLdata,1,replace=FALSE)}}
  return(m)
}

#### this matrix builder can build a matrix based on a covariate with given correlation
#### to VL levels using different methods. Ordering and filling of the matrices is
### based solely off of the randomly generated covariate values with given correlation
## filltype's are: rnd - doesn't use covariate values
## diag - puts highest values on largest diagonals
## corner - puts highest values in top left corner
## row - puts highest values in lowest row number from left to right
## border - puts highest values around the border (2 outer borders only)
one.cov.matrix <- function(data, matsize, SE, filltype, Uganda=FALSE){
  
  if (filltype=="rnd"){
    m <- rnd.mat.fill(data, matsize, SE, Uganda)}
  else if (filltype=="diag"){
    data <- data[order(data[,"predVL"], decreasing=TRUE),]
    m <- diag.mat.fill(VL, matsize)}
  else if (filltype=="corner"){
    data <- data[order(data[,"predVL"], decreasing=TRUE),]
    m <- corner.mat.fill(data, matsize, SE)}
  else if (filltype=="row"){
    data <- data[order(data[,"predVL"], decreasing=TRUE),]
    m <- row.mat.fill(VL, matsize)}
  else if (filltype=="border"){
    data <- data[order(data[,"predVL"], decreasing=TRUE),]
    m <- border.mat.fill(VL, matsize)}
  else {m <- rnd.mat.fill(data, matsize)}
  ### creates data frame to reduce from our matrix m
  x <- data.frame(m[,1:(matsize*2)])
  g <- data.frame(m[,(matsize*2+1):(matsize*3)])
  rowst <- replicate(matsize, 0)
  for (i in 1:matsize){rowst[i] = sum(x[i,1:matsize]/matsize)}
  rows <- rowst*(10^(rnorm(1, sd=SE)))
  colst <- replicate(matsize, 0)
  for (i in 1:matsize){colst[i] = sum(x[1:matsize,i]/matsize)}
  cols <- colst*(10^(rnorm(1, sd=SE)))
  return(cbind(x,rows,cols,g))
}



### Sums Smarter Method

#### This function tests the matrix squares corresponding to the highest row and column totals
#### making sure not to test the same square twice by zeroing out the tested square
#### this function returns the reduced matrix with a column for how many tests were conducted
### the reduced matrix shows which boxes were tested but does not zero out any columns or rows
### this reduced matrix will be used to compare against the original to compute sensitivity
### This function also accepts input tstperrd which indicates the number of boxes tested per round
### It needs the above function, check.indices to work
###### Modified Simple Search Method (Smt sums)

reduce.mat.smt <- function(x, y, matsize, cutoff, tstperd, lowlimit){
  # tracks number of tests
  t <- matsize*2
  # tracks number of rounds
  rounds <- 1
  # diagnostic; irrelevant
  check <- matrix(nrow=0, ncol=2)
  # Calculates prevaence metrics for diagnostics
  prev <- prev(x, matsize, cutoff)
  # function that checks every row and column average.  If less than t/n, 
  # this function classifies all patients in that row/column as not experiencing
  # treatment failure and zeros out that row/column
  x <- zerocolrow(x,matsize,cutoff, lowlimit)
  
  # Start of the while loop which determines when to stop testing
  while (max(x$cols) > cutoff/matsize & max(x$rows) > cutoff/matsize){
    # tracks how many rows and columns can contain failures
    num.rows.fail <- length(x$rows[x$rows > cutoff/matsize])
    num.cols.fail <- length(x$cols[x$cols > cutoff/matsize])
    # calculates the minimum of the above two numbers
    # useful in choosing how many tests to conduct
    min.num.fails <- min(num.rows.fail, num.cols.fail)
    # function that creates data frame of all possible sums of row and column totals
    # indexed by the ith row and jth column corresponding to the ith, jth sample.
    # If a patient has been classified, their sum will be zero regardless of their
    # corresponding row and column averages
    sums <- sums(x,matsize)
    z <- NULL
    d <- matrix(nrow=0, ncol=2)
    # orders the sums in a list in descending order
    z <- order(sums, decreasing=TRUE)
    # extracts the indicies corresponding to each sum total
    # only extracts indices of patients not yet classified
    for (k in 1:length(z)){
      j <- ceiling(z[k]/matsize)
      if (z[k] %% matsize == 0){i <- matsize}
      else{i <- z[k] %% matsize}
      if (x[i,j] == 0){}
      else{d <- rbind(d, c(i,j))}
    }
    # d is a dataframe of the indices of all unclassified patients from highest sum
    # to lowest
    d <- data.frame(d)
    # the number of unclassfied patients remaining
    max.num.fails <- length(d[,1])
    
    # The next three if/else statements choose the number of tests to conduct this
    # round
    
    # This determines if the number of unclassified patients is less than
    # or equal to our maximum number of tests per round
    if (max.num.fails <= tstperd){
      for (i in 1:length(d[,1])){
        # extracts the indices from data frame d
        a <- d[i,1]
        b <- d[i,2]
        # extracts the tested value (includes measurement error)
        xerror <- x[[a, (b+matsize)]]/matsize
        # Subtracts the appropriate value for the corresponding row and column
        x$rows[a] <- x$rows[a] - xerror
        x$cols[b] <- x$cols[b] - xerror
        # zeros out the true value and tested value from the x matrix.
        # This allows us to track which patients have been classified.
        x[[a,(matsize+b)]] = 0
        x[[a,b]] = 0
        # zeros out the true value in the y matrix to indicate that we tested
        # this patient.  Useful for calculating sensitivity in a later function.
        y[[a,b]] = 0
        # increases our number of tests by 1.
        t <- t+1
        # re-declares d as a data frame.  Important for when d has only 1 row.
        d <- data.frame(d)
      }
    }
    # If we pass the above check, this determines if the minimum of the number of 
    # rows/columns which could contain a failure is less than the maximum number of
    # tests divided by 2.  If so, we relax the restriction that all patients
    # tested be in distinct rows and columns.
    else if (min.num.fails < tstperd/2){
      for (i in 1:tstperd){
        # extracts the indices from data frame d
        a <- d[i,1]
        b <- d[i,2]
        # extracts the tested value (includes measurement error)
        xerror <- x[[a, (b+matsize)]]/matsize
        # Subtracts the appropriate value for the corresponding row and column
        x$rows[a] <- x$rows[a] - xerror
        x$cols[b] <- x$cols[b] - xerror
        # zeros out the true value and tested value from the x matrix.
        # This allows us to track which patients have been classified.
        x[[a,(matsize+b)]] = 0
        x[[a,b]] = 0
        # zeros out the true value in the y matrix to indicate that we tested
        # this patient.  Useful for calculating sensitivity in a later function.
        y[[a,b]] = 0
        # increases our number of tests by 1.
        t <- t+1
      }
    }
    # If we pass the above two checks, we will test the maximum number of tests
    # ensuring that the samples tested are in distinct rows and columns, starting
    # with the highest sum
    else{
      for (i in 1:tstperd){
        if (length(d[,1]) == 0){}
        else{
          # extracts the indices from data frame d
          a <- d[1,1]
          b <- d[1,2]
          # extracts the tested value (includes measurement error)
          xerror <- x[[a, (b+matsize)]]/matsize
          # Subtracts the appropriate value for the corresponding row and column
          x$rows[a] <- x$rows[a] - xerror
          x$cols[b] <- x$cols[b] - xerror
          # zeros out the true value and tested value from the x matrix.
          # This allows us to track which patients have been classified.
          x[[a,(matsize+b)]] = 0
          x[[a,b]] = 0
          # zeros out the true value in the y matrix to indicate that we tested
          # this patient.  Useful for calculating sensitivity in a later function.
          y[[a,b]] = 0
          # This removes all sample indices in the same row or column as the
          # sample neing tested.
          d <- data.frame(d[d[,1] != a & d[,2] != b,])
          # increases our number of tests by 1.
          t <- t+1
        }
      }
    }
    # Increases our number of testing rounds by 1
    rounds = rounds+1
    # function that checks every row and column average.  If less than t/n, 
    # this function classifies all patients in that row/column as not experiencing
    # treatment failure and zeros out that row/column
    x <- zerocolrow(x,matsize,cutoff,lowlimit)
    # Function that checks if we have classified all patients in a row/column
    # and if so, zeros out that row/column
    x <- checkrowcol(x,matsize)
  }
  # this returns matrix y with the true values zeroed out for every patient
  # we actually tested (not just classified).  A later function will compare these to 
  # the observable values to calculate sensitivity.  It also attachs extra columns onto 
  # matrix y as a way of also returning number of tests,rounds and prev metrics
  result <- data.frame(cbind(y,t, rounds, prev))
  
  return(result)
}

#### This function takes a random matrix and returns the number of tests conducted
### as well as the true and error sensitivity
test.one.smt <- function(data, matsize, cutoff, SE, tstperd, lowlimit, Uganda=FALSE){
  ## creates a matrix of random entries from the population
  x <- one.cov.matrix(data, matsize, SE, filltype="rnd", Uganda)
  ## 2 exact copies of the original matrix with rows and cols columns
  y <- x
  z <- x
  ## returns a matrix with zero'ed out boxes which were tested and 
  ## number of tests as the last column (matsize + 3)
  y <- reduce.mat.smt(x, y, matsize, cutoff, tstperd, lowlimit)
  t <- y[[1, "t"]]
  rounds <- y[[1, "rounds"]]
  tprevfail <- y[[1, "tprevfail"]]
  eprevfail <- y[[1, "eprevfail"]]
  x <- sens(y, z, matsize, cutoff)
  
  result <- data.frame(cbind(x,t, rounds, tprevfail, eprevfail))
  
  result <- result %>% 
    rename(
      mss_tsens = true.sens,
      mss_esens = error.sens,
      mss_t = t,
      mss_rds = rounds,
      mss_tprevfail = tprevfail,
      mss_eprevfail = eprevfail
    )
  return(result) 
}


################ this matrix takes in a testing matrix (can be partially tested)
################ and returns an estimated matrix where the predictions are fitted to the observed row and column
################ averages

mat.linreg2 <- function(data, mat, guess, matsize, prec, precrd){
  guess1 <- guess
  rows <- NULL
  cols <- NULL
  rowpercent <- matrix(nrow=matsize, ncol=matsize)
  colpercent <- matrix(nrow=matsize, ncol=matsize)
  totpercent <- matrix(nrow=matsize, ncol=matsize)
  for(i in 1:matsize){
    for (j in 1:matsize){
      if (mat[i,j]==0){guess1[i,j]=0}
      else{
        guess1[i,j] <- guess[i,j]
      }}}
  for (i in 1:matsize){
    rows[i] <- sum(guess1[i,])
    cols[i] <- sum(guess1[,i])
  }
  for (i in 1:matsize){
    for (j in 1:matsize){
      rowpercent[i,j] <- guess1[i,j]/rows[i]
      colpercent[i,j] <- guess1[i,j]/cols[j]
      totpercent[i,j] <- guess1[i,j]/sum(guess1[,1:matsize])
    }
  }
  for (i in 1:matsize){
    for (j in 1:matsize){
      if (guess1[i,j]==0){guess1[i,j]=0}
      else{
        guess1[i,j] <- (rowpercent[i,j]*matsize*mat$rows[i] + colpercent[i,j]*matsize*mat$cols[j]
                       + totpercent[i,j]*matsize*sum(mat$rows))/3
      }}}
  for (i in 1:matsize){
    rows[i] <- sum(guess1[i,])/matsize
    cols[i] <- sum(guess1[,i])/matsize
  }
  counter <- 1
  while (counter <= precrd & (max(abs(rows - mat$rows)) > prec | max(abs(cols - mat$cols)) > prec)){
    
    for(i in 1:matsize){
      for(j in 1:matsize){
        if (guess1[i,j]==0){guess1[i,j]=0}
        else{
          guess1[i,j] <- guess1[i,j]*(mat$cols[j]/cols[j])}}}
    
    for (i in 1:matsize){
      rows[i] <- sum(guess1[i,])/matsize
      cols[i] <- sum(guess1[,i])/matsize
    }
    
    for(i in 1:matsize){
      for(j in 1:matsize){
        if (guess1[i,j]==0){guess1[i,j]=0}
        else{
          guess1[i,j] <- guess1[i,j]*(mat$rows[i]/rows[i])}}}
    
    for (i in 1:matsize){
      rows[i] <- sum(guess1[i,])/matsize
      cols[i] <- sum(guess1[,i])/matsize
    }
    counter <- counter+1
  }
  prec.max <- max(max(abs(rows - mat$rows)), max(abs(cols - mat$cols)))
  return(data.frame(cbind(guess1, prec.max)))
}

############## Start of the linreg method ##############################

reduce.mat.linreg <- function(data, x, y, guess, matsize, prec, precrd, 
                              cutoff, tstperd, lowlimit){
  ## tracks number of trials
  t <- matsize*2
  # tracks number of rounds
  rounds <- 1
  prev <- prev(x, matsize, cutoff)
  #### zeros out the column or row total and all corresponding elements if 
  ### row/column total < cutoff/matsize
  x <- zerocolrow(x,matsize,cutoff, lowlimit)
  prec.max <- NULL
  counter <- 1
  #### preserving the first predicted guesses for VL that we pass in with 'guess'
  guess_first <- guess
  while (max(x$cols) > cutoff/matsize & max(x$rows) > cutoff/matsize){
    # alters original matrix rows and cols for the soe method
    z <- soe.matrix(x)
    xerror = NULL
    tstperd <- min(tstperd, min(length(x$rows[x$rows>0]), length(x$cols[x$cols>0])))
    ### method which chooses the boxes to test
    guess <- mat.linreg2(data, z, guess_first, matsize, prec, precrd)
    prec.max[counter] <- guess[[1, (matsize+1)]]
    guess <- data.matrix(guess[, 1:matsize])
    counter <- counter + 1
    d <- matrix(nrow=0, ncol=2)
    for (i in 1:tstperd){
      num <- which.max(guess)
      if (num %% 10 == 0){k <- 10}
      else {k <- num %% 10}
      j <- ceiling(num/10)
      d <- rbind(d, c(k,j))
      guess[k,j] = 0
    }
    for (i in 1:length(d[,1])){
      xerror[i] <- (x[[d[i,1],(matsize+d[i,2])]])/matsize
      a <- d[i,1]
      x$rows[a] <- x$rows[a] - xerror[i]
      b <- d[i,2]
      x$cols[b] <- x$cols[b] - xerror[i]
      x[[a,(matsize+b)]] = 0
      x[[a,b]] = 0
      y[[a,b]] = 0}
    #### number of rounds of testing needed, includes the 2 to start
    rounds = rounds+1
    #### zeros out the column or row total and all corresponding elements if 
    ### row/column total < cutoff/matsize
    x <- zerocolrow(x,matsize,cutoff, lowlimit)
    ### this tests whether all elements in a row or column have been tested
    #### if so, it zeros out the row/column total
    x <- checkrowcol(x,matsize)
    t <- t + length(d[,1])}
    prec.maxmax <- max(prec.max)
  
  result <- data.frame(cbind(y,t, rounds, prev))
  
  return(result)
}

test.one.linreg <- function(data, matsize, prec, precrd, cutoff, tstperd, SE, lowlimit, filltyp, Uganda=FALSE){
  ## creates a matrix of random entries from the population
  x <- one.cov.matrix(data, matsize, SE, filltyp, Uganda)
  ## 2 exact copies of the original matrix with rows and cols columns
  y <- x[,1:((matsize*2)+2)]
  z <- x[,1:((matsize*2)+2)]
  
  guess <- x[,((matsize*2)+3):((matsize*3)+2)]
  ## returns a matrix with zero'ed out boxes which were tested and 
  ## number of tests as the last column (matsize + 3)
  y <- reduce.mat.linreg(data, x, y, guess, matsize, prec, precrd, cutoff, tstperd, lowlimit)
  t <- y[[1, "t"]]
  rounds <- y[[1, "rounds"]]
  tprevfail <- y[[1, "tprevfail"]]
  eprevfail <- y[[1, "eprevfail"]]
  x <- sens(y, z, matsize, cutoff)
  
  result <- data.frame(cbind(x,t, rounds, tprevfail, eprevfail))
  
  result <- result %>% 
    rename(
      linreg_tsens = true.sens,
      linreg_esens = error.sens,
      linreg_t = t,
      linreg_rds = rounds,
      linreg_tprevfail = tprevfail,
      linreg_eprevfail = eprevfail
    )
  return(result)  
}



##### solve soe matrix

### This function takes in an soe matrix and index list and solves it for the fewest
### possible failures

soe.solve.weights <- function(z, d, matsize, lowlimit){
  ### this counter tracks the indicies of the rows of d
  ### this code will iterate over the top rows of d until we have a solution to
  ### the matrix
  counter2 <- 1
  ###### solves the matrix using weights and pooled totals to assign VL values
  ## the weights are used to decide which boxes to put the row/column totals in
  
  ## fixed is a matrix which tracks the indices of boxes already solved
  fixed <- matrix(nrow=0, ncol=2)
  #### adds all zeroed out boxes to 'fixed' unless they're already there
  for (i in 1:matsize){
    for (j in 1:matsize){
      if (z[i,j] == 0){
        if (check.indices(fixed, c(i, j))){}
        else {fixed <- rbind(fixed, c(i,j))}
      }
      else {}
    }
  }
  
  while (length(fixed[,1]) < (matsize^2)){
    ### tracks how many boxes were 'zeroed' on a given interation
    num.boxes.zero <- 0
    
    ### checks the first box to see if the row or col total for that box is greater
    if ((z$rows[d[counter2,1]]) > (z$cols[d[counter2,2]])){
      
      ### if the row is higher, we use the col total
      num <- z$cols[d[counter2,2]]
      ### this applies the col total to the appropriate box and zeros out the other boxes
      for (i in 1:matsize){
        ### This ensures we don't zero a box we've already fixed in our solution
        if(check.indices(fixed, c(i, d[counter2,2]))){}
        else{
          ### this zeros out all the boxes not fixed in the col
          z[i, d[counter2,2]] = lowlimit
          ### tracks number of boxes zeroed
          num.boxes.zero <- num.boxes.zero + 1
          ### alters each row total by our lower VL limit/matsize for each box zeroed
          z$rows[i] <- z$rows[i] - lowlimit/matsize
          fixed <- rbind(fixed, c(i,  d[counter2,2]))}
      }
      
      box.num <- (num*matsize) - (lowlimit*(num.boxes.zero - 1))
      #### applies our solution VL number to the appropriate box
      z[[d[counter2,1], d[counter2,2]]]= box.num
      ##### also zeros out the corresponding 'error values' (Probably not relevant)
      z[[d[counter2,1],(matsize+d[counter2,2])]] = 0
      #### zeros out the col total
      z$cols[d[counter2,2]] <- 0
      ### subtracts the num from the row total and adds back in lowlimit/count.rowz that was subtracted earlier
      z$rows[d[counter2,1]] <- z$rows[d[counter2,1]] - (box.num/matsize) + lowlimit/matsize
      #### clears all coordinates in matrix d with the same col coordinate
      d <- subset(d, d[,2] != d[counter2,2])
      ### reorders matrix d by weight
      d <- d[order(d[,3], decreasing = TRUE),]
      ### ensures d is a matrix
      if (is.vector(d)){
        d <- matrix(d, nrow=1, ncol=3)}
      else{}
    }
    
    ### the same comments apply to when col total > row total
    else{
      
      #### extracts the row average
      num <- z$rows[d[counter2,1]]
      for (i in 1:matsize){
        #### ensrues we do not touch boxes where we have fixed values of 0 or otherwise
        if(check.indices(fixed, c(d[counter2,1], i))){}
        else{
          ### assigns fixed value zero to others in row
          z[d[counter2,1],i] = lowlimit
          ### tracks number of boxes zeroed
          num.boxes.zero <- num.boxes.zero + 1
          ### alters each row total by our lower VL limit/matsize for each box zeroed
          z$cols[i] <- z$cols[i] - lowlimit/matsize
          fixed <- rbind(fixed, c(d[counter2,1], i))}
      }
      box.num <- (num*matsize) - (lowlimit*(num.boxes.zero - 1))
      z[[d[counter2,1], d[counter2,2]]]= box.num
      z[[d[counter2,1],(matsize+d[counter2,2])]] = 0
      z$rows[d[counter2,1]] <- 0
      z$cols[d[counter2,2]] <- z$cols[d[counter2,2]] - (box.num/matsize) + lowlimit/matsize
      d <- subset(d, d[,1] != d[counter2,1])
      d <- d[order(d[,3], decreasing = TRUE),]
      if (is.vector(d)){
        d <- matrix(d, nrow=1, ncol=3)}
      else{}
      
    }
    #### This ends our code for solving our z matrix
  } 
  return(z)
}

##### This is the linreg soe method.  It is the reverse of the soelinreg method
##### Here we assign weights based on linreg first
##### then we solve the soe matrix using those weights to choose solution

##### SOE Full with Linreg Method w/testing round limit

reduce.mat.linregsoe <- function(data, x, y, guess, matsize, cutoff,
                                 prec, precrd, lowlimit){
  
  
  ## tracks number of trials
  t <- matsize*2
  # tracks number of rounds
  rounds <- 1
  prev <- prev(x, matsize, cutoff)
  x <- zerocolrow(x,matsize,cutoff, lowlimit)
  prec.max <- NULL
  counter <- 1
  check <- matrix(nrow=0, ncol=2)
  ##### preserving the initial VL guess
  guess_first <- guess
  while (max(x$cols) > cutoff/matsize & max(x$rows) > cutoff/matsize){
    # alters original matrix rows and cols for the soe method
    z <- soe.matrix(x)
    d <- matrix(nrow=0, ncol=3)
    xerror = NULL
    ##### This chunk solves the linreg matrix and assigns weights
    guess <- mat.linreg2(data, z, guess_first, matsize, prec, precrd)
    ### tracks the distance btwn the guess matrix row/col total and observed values
    prec.max[counter] <- guess[[1, (matsize+1)]]
    guess <- data.matrix(guess[, 1:matsize])
    ### this counter is for prec.max tracking
    counter <- counter + 1
    ### sum of all guessed VL values
    total <- sum(guess)
    #### this converts VL values to weights and tracks the data in matrix d
    for (i in 1:matsize){
      for (j in 1:matsize){
        guess[i,j] <- guess[i,j]/total
        d <- rbind(d, c(i,j,guess[i,j]))
      }
    }
    ### this orders matrix d by highest weight first
    d <- d[order(d[,3], decreasing = TRUE),]
    ### This code ensures that d is a matrix, even if only 1 row
    if (is.vector(d)){
      d <- matrix(d, nrow=1, ncol=3)}
    else{}
    #### this sollves our soe matrix based on the weights
    z <- soe.solve.weights(z, d, matsize, lowlimit)
    
    ### start of our code to test VL values
    #### extracts only our solution matrix used to decide which boxes to test
    testmat <- data.frame(z[,1:matsize])
    ### tracks number of boxes > 0 in our solution (testable boxes)
    over.limit <- testmat[testmat>lowlimit]
    index <- NULL
    # loops over the number of tests this round where max tests in number in our solution
    for (k in 1:length(over.limit)){
      ## extracts the list coordinates of the highest box and converts to matrix coords
      index <- head(order(testmat, decreasing=TRUE))
      ### this code converst the list coordinates into matrix coords
      j <- ceiling(index[1]/matsize)
      if (index[1] %% 10 == 0){i <- 10}
      else{i <- index[1] %% 10}
      ### zeros our our testmat box to ensure we do not test twice
      testmat[i,j] <- 0
      ### diagnostic code, keeps running total of all boxes checked
      check <- rbind(check, c(i,j))
      ### extracts the error value for the coordinate and divides by the matsize
      xerror <- (x[[i,(matsize+j)]])/matsize
      x$rows[i] <- x$rows[i] - xerror
      x$cols[j] <- x$cols[j] - xerror
      x[[i,(matsize+j)]] = 0
      x[[i,j]] = 0
      y[[i,j]] = 0
    }
    #### number of rounds of testing needed, includes the 2 to start
    rounds = rounds+1
    #### zeros out the column or row total and all corresponding elements if 
    ### row/column total < cutoff/matsize
    x <- zerocolrow(x,matsize,cutoff, lowlimit)
    ### this tests whether all elements in a row or column have been tested
    #### if so, it zeros out the row/column total
    x <- checkrowcol(x,matsize)
    t <- t + length(over.limit)}
  prec.maxmax <- max(prec.max)
  
  result <- data.frame(cbind(y,t, rounds, prev))
  
  return(result)
}

test.one.linregsoe <- function(data, matsize, prec, precrd,
                               cutoff, SE, lowlimit, filltyp, Uganda=FALSE){
  
  ## creates a matrix of random entries from the population
  x <- one.cov.matrix(data, matsize, SE, filltyp, Uganda)
  ## 2 exact copies of the original matrix with rows and cols columns
  y <- x[,1:((matsize*2)+2)]
  z <- x[,1:((matsize*2)+2)]
  
  guess <- x[,((matsize*2)+3):((matsize*3)+2)]
  ## returns a matrix with zero'ed out boxes which were tested and 
  ## number of tests as the last column (matsize + 3)
  y <- reduce.mat.linregsoe(data, x, y, guess, matsize, cutoff,
                            prec, precrd, lowlimit)
  t <- y[[1, "t"]]
  rounds <- y[[1, "rounds"]]
  tprevfail <- y[[1, "tprevfail"]]
  eprevfail <- y[[1, "eprevfail"]]
  x <- sens(y, z, matsize, cutoff)
  
  result <- data.frame(cbind(x,t, rounds, tprevfail, eprevfail))
  
  result <- result %>% 
    rename(
      lrsoe_tsens = true.sens,
      lrsoe_esens = error.sens,
      lrsoe_t = t,
      lrsoe_rds = rounds,
      lrsoe_tprevfail = tprevfail,
      lrsoe_eprevfail = eprevfail
    )
  return(result)  
}



#### This is the minipool method, using a full size matrix, but only
#### testing rows individually
#### this allows for 1 test per row per round, reducing the overall number of rounds
#### this compares apples to apples for comparing to matrix methods

##### Minipool w/covariates

zero.row.mini <- function(x, matsize, cutoff){
  for (i in 1:matsize){
    if (x$rows[i] < cutoff/matsize){
      x$rows[i] = 0
      for (j in 1:matsize){
        x[i,j] = 0
      }
    }
  }
  return(x)
}

check.row.mini <- function(x, matsize){
  for (i in 1:matsize){
    if(sum(x[i,1:matsize]) == 0){x$rows[i] = 0}
    else{}
  }
  return(x)
}


mat.linreg.mini2 <- function(data, mat, guess, matsize){
  guess1 <- guess
  
  for(i in 1:matsize){
    for (j in 1:matsize){
      if (mat[i,j]==0){guess1[i,j]=0}
      else{
        guess1[i,j] <- guess[i,j]
      }}}
  
  return(guess1)
}

reduce.mat.mini.cov <- function(data, x, y, guess, matsize, cutoff){
  ### tracks number of tests
  t <- matsize
  rounds <- 1
  check <- matrix(nrow=0, ncol=2)
  prev <- prev(x, matsize, cutoff)
  x <- zero.row.mini(x, matsize, cutoff)
  ###### preserving the initial VL guess
  guess_first <- guess
  while (max(x$rows) > (cutoff/matsize)){
    xerror <- NULL
    guess <- mat.linreg.mini2(data, x, guess_first, matsize)
    
    for (i in 1:matsize){
      if (x$rows[i] == 0){}
      else{
        index <- which.max(guess[i,])
        xerror <- x[[i,(index+matsize)]]/matsize
        x$rows[i] <- x$rows[i] - xerror
        x[[i,(index+matsize)]] <- 0
        x[i, index] = 0
        y[i, index] = 0
        t <- t+1
        check <- rbind(check, c(i,index))
      }
    }
    rounds <- rounds+1
    x <- zero.row.mini(x, matsize, cutoff)
    x <- check.row.mini(x, matsize)
  }
  result <- data.frame(cbind(y,t, rounds, prev))
  
  return(result)
}


test.one.mini.cov <- function(data, matsize, cutoff, SE, filltyp, Uganda=FALSE){
  #### creates a matrix of VL values with row and colum totals
  x <- one.cov.matrix(data, matsize, SE, filltyp, Uganda)
  #### this strips off the column totals
  y <- x[,1:((matsize*2)+2)]
  z <- x[,1:((matsize*2)+2)]
  guess <- x[,((matsize*2)+3):((matsize*3)+2)]
  
  y <- reduce.mat.mini.cov(data, x, y, guess, matsize, cutoff)
  t <- y[[1, "t"]]
  rounds <- y[[1, "rounds"]]
  tprevfail <- y[[1, "tprevfail"]]
  eprevfail <- y[[1, "eprevfail"]]
  x <- sens(y, z, matsize, cutoff)
  
  result <- data.frame(cbind(x,t, rounds, tprevfail, eprevfail))
  
  result <- result %>% 
    rename(
      mincov_tsens = true.sens,
      mincov_esens = error.sens,
      mincov_t = t,
      mincov_rds = rounds,
      mincov_tprevfail = tprevfail,
      mincov_eprevfail = eprevfail
    )
  return(result)  
}

reduce.mat.mini <- function(data, x, y, matsize, cutoff){
  ### tracks number of tests
  t <- matsize
  rounds <- 1
  check <- matrix(nrow=0, ncol=2)
  prev <- prev(x, matsize, cutoff)
  x <- zero.row.mini(x, matsize, cutoff)
  
  while (max(x$rows) > (cutoff/matsize)){
    xerror <- NULL
    
    for (i in 1:matsize){
      if (x$rows[i] == 0){}
      else{
        index <- 1
        while (x[i, index] == 0){
          index <- index + 1
        }
        xerror <- x[[i,(index+matsize)]]/matsize
        x$rows[i] <- x$rows[i] - xerror
        x[[i,(index+matsize)]] <- 0
        x[i, index] = 0
        y[i, index] = 0
        t <- t+1
        check <- rbind(check, c(i,index))
      }
    }
    rounds <- rounds+1
    x <- zero.row.mini(x, matsize, cutoff)
    x <- check.row.mini(x, matsize)
  }
  result <- data.frame(cbind(y,t, rounds, prev))
  
  return(result)
}


test.one.mini <- function(data, matsize, cutoff, SE, Uganda=FALSE){
  #### creates a matrix of VL values with row and colum totals
  x <- one.cov.matrix(data, matsize, SE, filltype="rnd", Uganda)
  #### this strips off the column totals
  y <- x
  z <- x
  y <- reduce.mat.mini(data, x, y, matsize, cutoff)
  
  t <- y[[1, "t"]]
  rounds <- y[[1, "rounds"]]
  tprevfail <- y[[1, "tprevfail"]]
  eprevfail <- y[[1, "eprevfail"]]
  x <- sens(y, z, matsize, cutoff)
  
  result <- data.frame(cbind(x,t, rounds, tprevfail, eprevfail))
  
  result <- result %>% 
    rename(
      mini_tsens = true.sens,
      mini_esens = error.sens,
      mini_t = t,
      mini_rds = rounds,
      mini_tprevfail = tprevfail,
      mini_eprevfail = eprevfail
    )
  return(result)   
}


####################### The testing functionsare below ############################

pool.alg.cov <- function(reps, data, matsize, prec, precrd,
                               cutoff, SE, tstperd, lowlimit, filltyp, Uganda=FALSE){
  tstperrd <- tstperd
  
  linreg <- NULL
  smart <- NULL
  linregsoe <- NULL
  mini.cov <- NULL
  mini.pool <- NULL
  
  #data$VLobs <- data$VL*10^(rnorm(length(data$VL), sd=SE))
  
  for (i in 1:(reps)){
    check <- NULL
    index <- sample.int(n=length(data$VL), size=(matsize^2), replace=FALSE)
    data1 <- data[c(index),]
    set.seed(1)
    linreg1 <- test.one.linreg(data1, matsize, prec, precrd, cutoff, tstperd, SE, lowlimit, filltyp, Uganda)
    
    set.seed(1)
    smart1 <- test.one.smt(data1, matsize, cutoff, SE, tstperrd, lowlimit, Uganda)
    
    set.seed(1)
    linregsoe1 <- test.one.linregsoe(data1, matsize, prec, precrd, cutoff, SE, lowlimit, filltyp, Uganda)
    
    set.seed(1)
    mini.cov1 <- test.one.mini.cov(data1, matsize, cutoff, SE, filltyp, Uganda)
    
    set.seed(1)
    mini.pool1 <- test.one.mini(data1, matsize, cutoff, SE, Uganda)
    
    ### removes tested subjects from population    
    data <- data[-c(index),]
    
    linreg <- rbind(linreg, linreg1)
    smart <- rbind(smart, smart1)
    linregsoe <- rbind(linregsoe, linregsoe1)
    mini.cov <- rbind(mini.cov, mini.cov1)
    mini.pool <- rbind(mini.pool, mini.pool1)
    
    # check <- data.frame(cbind(linreg, smart, linregsoe, mini.cov, mini.pool))
    # write.table(check, file="C:/Users/Barny/Dropbox/KI_Project_4/Results/temp//Results_reverse_ME.5_rand.R")
  
    }
  
  # check <- NULL
  # index <- sample.int(n=length(data$VL), size=(matsize^2), replace=FALSE)
  # data1 <- data[c(index),]
  # set.seed(1)
  # linreg1 <- test.one.linreg(data1, matsize, prec, precrd, cutoff, tstperd, SE, lowlimit, filltyp)
  # 
  # set.seed(1)
  # smart1 <- test.one.smt(data1, matsize, cutoff, SE, tstperrd, lowlimit)
  # 
  # set.seed(1)
  # linregsoe1 <- test.one.linregsoe(data1, matsize, prec, precrd, cutoff, SE, lowlimit, filltyp)
  # 
  # set.seed(1)
  # mini.cov1 <- test.one.mini.cov(data1, matsize, cutoff, SE)
  # 
  # set.seed(1)
  # mini.pool1 <- test.one.mini(data1, matsize, cutoff, SE)
  # 
  # ### removes tested subjects from population    
  # data <- data[-c(index),]

  
  result <- data.frame(cbind(linreg, smart, 
                             linregsoe, 
                             mini.cov, mini.pool))
  
  return(result)
  
}


####### this function predicts VL based on the simulation model

predictVL <- function(data, b0star, b1star, b2star, b3star){
  data$predVL <- 10^(b0star + b1star*data$adhere + b2star*data$pf + b3star*data$adhere*data$pf)
  return(data)
}


hypred <- function(reps, data, matsize, prec, precrd,
                   cutoff, SE, tstperd, lowlimit, filltyp, top_percent, bot_percent){
  
  samples <- reps*(matsize^2)
  data <- data[1:samples,]
  
  tstperrd <- tstperd
  
  linreg <- NULL
  smart <- NULL
  linregsoe <- NULL
  mini.cov <- NULL
  mini.pool <- NULL
  
  linreg_low <- NULL
  smart_low <- NULL
  linregsoe_low <- NULL
  mini.cov_low <- NULL
  mini.pool_low <- NULL
  
  all <- NULL
  
  data <- data[order(data[,"predVL"], decreasing=TRUE),]
  
  a <- length(data$VL)*top_percent
  b <- length(data$VL)*(1-bot_percent)
  
  top_tier <- data[1:a,]
  mid_tier <- data[(a+1):b,]
  low_tier <- data[(b+1):length(data$VL),]
  
  mid_tier_reps <- length(mid_tier$VL)/(matsize^2)
  low_tier_reps <- length(low_tier$VL)/(matsize^2)
  
  top_tier$error <- top_tier$VL*(10^(rnorm(length(top_tier$VL), sd=SE)))
  
  all_tprevfail <- sum(top_tier$VL >=1000)
  all_eprevfail <- sum(top_tier$error >=1000)
  all_tsens <- sum(top_tier$VL >= 1000 & top_tier$error >=1000)/all_tprevfail
  all_esens <- sum(top_tier$error >=1000)/all_eprevfail
  all_t <- length(top_tier$VL)
  all_rds <- 1
  section <- "top"
  
  linreg_tsens <- NA
  linreg_esens <- NA
  linreg_t <- NA
  linreg_rds <- NA
  linreg_tprevfail <- NA
  linreg_eprevfail <- NA
  
  mss_tsens <- NA
  mss_esens <- NA
  mss_t <- NA
  mss_rds <- NA
  mss_tprevfail <- NA
  mss_eprevfail <- NA
  
  lrsoe_tsens <- NA
  lrsoe_esens <- NA
  lrsoe_t <- NA
  lrsoe_rds <- NA
  lrsoe_tprevfail <- NA
  lrsoe_eprevfail <- NA
  
  mincov_tsens <- NA
  mincov_esens <- NA
  mincov_t <- NA
  mincov_rds <- NA
  mincov_tprevfail <- NA
  mincov_eprevfail <- NA
  
  mini_tsens <- NA
  mini_esens <- NA
  mini_t <- NA
  mini_rds <- NA
  mini_tprevfail <- NA
  mini_eprevfail <- NA
  
  temp <- data.frame(cbind(linreg_tsens, linreg_esens, linreg_t, linreg_rds, linreg_tprevfail, linreg_eprevfail,
                                      mss_tsens, mss_esens, mss_t, mss_rds, mss_tprevfail, mss_eprevfail,
                                      lrsoe_tsens, lrsoe_esens, lrsoe_t, lrsoe_rds, lrsoe_tprevfail, lrsoe_eprevfail,
                                      mincov_tsens, mincov_esens, mincov_t, mincov_rds, mincov_tprevfail, mincov_eprevfail,
                                      mini_tsens, mini_esens, mini_t, mini_rds, mini_tprevfail, mini_eprevfail))
  
  temp <- data.frame(lapply(temp, function(x) as.numeric(x)))
  
  
  result_top <- data.frame(cbind(temp,
                                all_tsens, all_esens, all_t, all_rds, all_tprevfail, all_eprevfail, section))
  
  for (i in 1:mid_tier_reps){
    check <- NULL
    index <- sample.int(n=length(mid_tier$VL), size=(matsize^2), replace=FALSE)
    data1 <- mid_tier[c(index),]
    set.seed(1)
    linreg1 <- test.one.linreg(data1, matsize, prec, precrd, cutoff, tstperd, SE, lowlimit, filltyp)
    
    set.seed(1)
    smart1 <- test.one.smt(data1, matsize, cutoff, SE, tstperrd, lowlimit)
    
    set.seed(1)
    linregsoe1 <- test.one.linregsoe(data1, matsize, prec, precrd, cutoff, SE, lowlimit, filltyp)
    
    set.seed(1)
    mini.cov1 <- test.one.mini.cov(data1, matsize, cutoff, SE)
    
    set.seed(1)
    mini.pool1 <- test.one.mini(data1, matsize, cutoff, SE)
    
    ### removes tested subjects from population    
    
    mid_tier <- mid_tier[-c(index),]
    
    linreg <- rbind(linreg, linreg1)
    smart <- rbind(smart, smart1)
    linregsoe <- rbind(linregsoe, linregsoe1)
    mini.cov <- rbind(mini.cov, mini.cov1)
    mini.pool <- rbind(mini.pool, mini.pool1)
    
  }
  
  section <- "mid"
  
  all_tsens <- NA
  all_esens <- NA
  all_t <- NA
  all_rds <- NA
  all_tprevfail <- NA
  all_eprevfail <- NA
  
  temp <- data.frame(cbind(all_tsens, all_esens, all_t, all_rds, all_tprevfail, all_eprevfail))
  
  temp <- data.frame(lapply(temp, function(x) as.numeric(x)))
  
  result_mid <- data.frame(cbind(linreg, smart, linregsoe, mini.cov, mini.pool, 
                                 temp,
                                 section))
  
  for (i in 1:low_tier_reps){
    check <- NULL
    index <- sample.int(n=length(low_tier$VL), size=(matsize^2), replace=FALSE)
    data1 <- low_tier[c(index),]
    set.seed(1)
    linreg1 <- test.one.linreg(data1, matsize, prec, precrd, cutoff, tstperd, SE, lowlimit, filltyp)
    
    set.seed(1)
    smart1 <- test.one.smt(data1, matsize, cutoff, SE, tstperrd, lowlimit)
    
    set.seed(1)
    linregsoe1 <- test.one.linregsoe(data1, matsize, prec, precrd, cutoff, SE, lowlimit, filltyp)
    
    set.seed(1)
    mini.cov1 <- test.one.mini.cov(data1, matsize, cutoff, SE)
    
    set.seed(1)
    mini.pool1 <- test.one.mini(data1, matsize, cutoff, SE)
    
    ### removes tested subjects from population    
    
    low_tier <- low_tier[-c(index),]
    
    linreg_low <- rbind(linreg_low, linreg1)
    smart_low <- rbind(smart_low, smart1)
    linregsoe_low <- rbind(linregsoe_low, linregsoe1)
    mini.cov_low <- rbind(mini.cov_low, mini.cov1)
    mini.pool_low <- rbind(mini.pool_low, mini.pool1)
    
  }
  
  section <- "low"
  
  result_low <- data.frame(cbind(linreg_low, smart_low, linregsoe_low, mini.cov_low, mini.pool_low, 
                                 temp,
                                 section))
  
  result <- data.frame(rbind(result_top, result_mid, result_low))
  
  return(result)
  
}




hypred_uganda <- function(reps, data, matsize, prec, precrd,
                   cutoff, SE, tstperd, lowlimit, filltyp, Uganda=FALSE){
  
  
  tstperrd <- tstperd
  
  linreg <- NULL
  smart <- NULL
  linregsoe <- NULL
  mini.cov <- NULL
  mini.pool <- NULL
  
  linreg_low <- NULL
  smart_low <- NULL
  linregsoe_low <- NULL
  mini.cov_low <- NULL
  mini.pool_low <- NULL
  
  all <- NULL
  
  data <- data[order(data[,"predVL"], decreasing=TRUE),]
  
  top_tier <- data[1:300,]
  mid_tier <- data[301:1300,]
  low_tier <- data[1301:3600,]
  
  mid_tier_reps <- floor(length(mid_tier$VL)/(matsize^2))
  low_tier_reps <- floor(length(low_tier$VL)/(matsize^2))
  
  top_tier$error <- top_tier$VL
  
  all_tprevfail <- sum(top_tier$VL >=1000)
  all_eprevfail <- sum(top_tier$error >=1000)
  all_tsens <- sum(top_tier$VL >= 1000 & top_tier$error >=1000)/all_tprevfail
  all_esens <- sum(top_tier$error >=1000)/all_eprevfail
  all_t <- length(top_tier$VL)
  all_rds <- 1
  section <- "top"
  
  linreg_tsens <- NA
  linreg_esens <- NA
  linreg_t <- NA
  linreg_rds <- NA
  linreg_tprevfail <- NA
  linreg_eprevfail <- NA
  
  mss_tsens <- NA
  mss_esens <- NA
  mss_t <- NA
  mss_rds <- NA
  mss_tprevfail <- NA
  mss_eprevfail <- NA
  
  lrsoe_tsens <- NA
  lrsoe_esens <- NA
  lrsoe_t <- NA
  lrsoe_rds <- NA
  lrsoe_tprevfail <- NA
  lrsoe_eprevfail <- NA
  
  mincov_tsens <- NA
  mincov_esens <- NA
  mincov_t <- NA
  mincov_rds <- NA
  mincov_tprevfail <- NA
  mincov_eprevfail <- NA
  
  mini_tsens <- NA
  mini_esens <- NA
  mini_t <- NA
  mini_rds <- NA
  mini_tprevfail <- NA
  mini_eprevfail <- NA
  
  temp <- data.frame(cbind(linreg_tsens, linreg_esens, linreg_t, linreg_rds, linreg_tprevfail, linreg_eprevfail,
                           mss_tsens, mss_esens, mss_t, mss_rds, mss_tprevfail, mss_eprevfail,
                           lrsoe_tsens, lrsoe_esens, lrsoe_t, lrsoe_rds, lrsoe_tprevfail, lrsoe_eprevfail,
                           mincov_tsens, mincov_esens, mincov_t, mincov_rds, mincov_tprevfail, mincov_eprevfail,
                           mini_tsens, mini_esens, mini_t, mini_rds, mini_tprevfail, mini_eprevfail))
  
  temp <- data.frame(lapply(temp, function(x) as.numeric(x)))
  
  
  result_top <- data.frame(cbind(temp,
                                 all_tsens, all_esens, all_t, all_rds, all_tprevfail, all_eprevfail, section))
  
  for (i in 1:mid_tier_reps){
    check <- NULL
    index <- sample.int(n=length(mid_tier$VL), size=(matsize^2), replace=FALSE)
    data1 <- mid_tier[c(index),]
    set.seed(1)
    linreg1 <- test.one.linreg(data1, matsize, prec, precrd, cutoff, tstperd, SE, lowlimit, filltyp, Uganda)
    
    set.seed(1)
    smart1 <- test.one.smt(data1, matsize, cutoff, SE, tstperrd, lowlimit, Uganda)
    
    set.seed(1)
    linregsoe1 <- test.one.linregsoe(data1, matsize, prec, precrd, cutoff, SE, lowlimit, filltyp, Uganda)
    
    set.seed(1)
    mini.cov1 <- test.one.mini.cov(data1, matsize, cutoff, SE, filltyp, Uganda)
    
    set.seed(1)
    mini.pool1 <- test.one.mini(data1, matsize, cutoff, SE, Uganda)
    
    ### removes tested subjects from population    
    
    mid_tier <- mid_tier[-c(index),]
    
    linreg <- rbind(linreg, linreg1)
    smart <- rbind(smart, smart1)
    linregsoe <- rbind(linregsoe, linregsoe1)
    mini.cov <- rbind(mini.cov, mini.cov1)
    mini.pool <- rbind(mini.pool, mini.pool1)
    
  }
  
  section <- "mid"
  
  all_tsens <- NA
  all_esens <- NA
  all_t <- NA
  all_rds <- NA
  all_tprevfail <- NA
  all_eprevfail <- NA
  
  temp <- data.frame(cbind(all_tsens, all_esens, all_t, all_rds, all_tprevfail, all_eprevfail))
  
  temp <- data.frame(lapply(temp, function(x) as.numeric(x)))
  
  result_mid <- data.frame(cbind(linreg, smart, linregsoe, mini.cov, mini.pool, 
                                 temp,
                                 section))
  
  for (i in 1:low_tier_reps){
    check <- NULL
    index <- sample.int(n=length(low_tier$VL), size=(matsize^2), replace=FALSE)
    data1 <- low_tier[c(index),]
    set.seed(1)
    linreg1 <- test.one.linreg(data1, matsize, prec, precrd, cutoff, tstperd, SE, lowlimit, filltyp, Uganda)
    
    set.seed(1)
    smart1 <- test.one.smt(data1, matsize, cutoff, SE, tstperrd, lowlimit, Uganda)
    
    set.seed(1)
    linregsoe1 <- test.one.linregsoe(data1, matsize, prec, precrd, cutoff, SE, lowlimit, filltyp, Uganda)
    
    set.seed(1)
    mini.cov1 <- test.one.mini.cov(data1, matsize, cutoff, SE, filltyp, Uganda)
    
    set.seed(1)
    mini.pool1 <- test.one.mini(data1, matsize, cutoff, SE, Uganda)
    
    ### removes tested subjects from population    
    
    low_tier <- low_tier[-c(index),]
    
    linreg_low <- rbind(linreg_low, linreg1)
    smart_low <- rbind(smart_low, smart1)
    linregsoe_low <- rbind(linregsoe_low, linregsoe1)
    mini.cov_low <- rbind(mini.cov_low, mini.cov1)
    mini.pool_low <- rbind(mini.pool_low, mini.pool1)
    
  }
  
  section <- "low"
  
  result_low <- data.frame(cbind(linreg_low, smart_low, linregsoe_low, mini.cov_low, mini.pool_low, 
                                 temp,
                                 section))
  
  result <- data.frame(rbind(result_top, result_mid, result_low))
  
  return(result)
  
}




