
# based on https://www.ncbi.nlm.nih.gov/pubmed/12005759
# Permutation entropy: a natural complexity measure for time series.


logistic.map <- function(r, x, N, M){
  # https://magesblog.com/post/2012-03-17-logistic-map-feigenbaum-diagram/
  ## r: bifurcation parameter
  ## x: initial value
  ## N: number of iteration
  ## M: number of iteration points to be returned
  z <- 1:N
  z[1] <- x
  for(i in c(1:(N-1))){
    z[i+1] <- r *z[i]  * (1 - z[i])
  }
  ## Return the last M iterations 
  z[c((N-M):N)]
}
## Set scanning range for bifurcation parameter r
my.r <- seq(2.5, 4, by=0.003)
system.time(Orbit <- sapply(my.r, logistic.map,  x=0.1, N=1000, M=300))
#user  system elapsed 
#0.649   0.007   0.657 
Orbit <- as.vector(Orbit)
r <- sort(rep(my.r, 301))
#plot(Orbit ~ r, pch='.')


randomWalk <- function(n) {
  # first we have some sequence
  seqvec <- vector(length = n)
  for (i in 1:n){
    if (runif(n=1) > 0.5) {
      seqvec[i] <- rnorm(mean=5,n=1)
    }
    else {
      seqvec[i] <- rnorm(mean=0,n=1)
    }
  }  
  return(seqvec)
}

symbol <- function(x) {
  # x is a short vector
  s <- c()
  for (i in 1:(length(x)-1)) {
    s <- c(s, as.character(as.numeric((x[i] > x[i+1]))))
  }
  paste(s,collapse="")
}

symbolTransform <- function(x, ord) {
  # takes a numeric vector
  # and returns a sequence of symbols
  # ord is the markov order
  symbList <- c()
  for (i in 1:(length(x)-ord)) {
    symbList <- c(symbList, symbol(x[i:(i+ord)]))    
  }
  symbList
}

symbEntropy <- function(x,ord,norm=F) {
  # takes a sequence of symbols returns entropy
  # symbs are the symbols
  # x is the time series
  # ord is the markov order
  # norm indicates that the entropy should be normalized per symbol of order n
  o <- ord-1 # to get two time steps, we need this step : step+(ord-1)
  symbs <- symbolTransform(x, ord=o)
  symbTable <- table(symbs)
  denom <- length(symbs)-o+1
  ent <- 0
  for (i in 1:length(symbTable)){
    ent <- ent - ((symbTable[i]/denom)*log2((symbTable[i]/denom)))
  }
  if (norm) {
    return(ent/(ord-1))
  } else {
    return(ent)    
  }
}

#x <- randomWalk(100)
x <- c(4,7,9,10,6,11,3)
#x <- Orbit # see logistic.map above
e <- symbEntropy(x, ord=2) # order 2 means symbols of two time steps.


# WOW, really resilient to noise.
enoise <- c()
x <- c(4,7,9,10,6,11,3)
for (ni in seq(from=0,to=1,by=0.1)) {
  xnoise <- x + rnorm(mean=0, sd=0.5, n=length(x))
  enoise <- c(enoise, symbEntropy(xnoise, ord=1))
}


###################
# Symbolic Transfer Entropy
# following paper: https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.100.158101

buildJoint <- function(xsym, ysym, n, l) {
  jointCounts <- table(xsym[1:(n-l+1)], ysym[1:(n-l+1)], xsym[l:n])
  jointProb <- jointCounts/sum(jointCounts)
  jointProb
}

safediv <- function(a,b) {
  # divide a/b
  # if b != 0
  # a can be a vector
  if (b == 0) {
    return(rep(0, length(a)))
  } else {
    a/b
  }
}

buildCond1 <- function(xsym, ysym, n, l){
  a <- xsym[1:(n-l+1)]
  b <- ysym[1:(n-l+1)]
  c <- xsym[l:n]
  countxy <- table(a,b)  # given x_i and y_i
  countxyz <- table(a,b,c)  # prob of x_(i+l)
  for (si in unique(a)){
    for (sj in unique(b)) {
      countxyz[si,sj,] <- safediv(countxyz[si,sj,], countxy[si,sj])
    }
  }
  countxyz
}

buildCond2 <- function(xsym, n, l){
  # want P(x_l | x)
  a <- xsym[1:(n-l+1)]
  b <- xsym[l:n]
  countxx <- table(a,b)
  countx  <- rowSums(countxx) 
  for (si in unique(a)) {
    countxx[si,] <- safediv(countxx[si,], countx[si])
  }
  countxx
}


ste <- function(y,x,ord,l) {
  # x,y are time series
  # ord is the markov order or embedding dimension
  # here we assume the direction is y -> x (y influencing x)
  # l is the time delay (for y -> time delay l -> x)
  # predicting x out into the future using current x and y.
  xsym <- symbolTransform(x,ord)
  ysym <- symbolTransform(y,ord)
  n <- length(xsym)
  jointProb <- buildJoint(xsym,ysym,n,l)  
  condProbTop <- buildCond1(xsym,ysym,n,l) 
  condProbBot <- buildCond2(xsym,n,l)
  a <- xsym[1:(n-l+1)]
  b <- ysym[1:(n-l+1)]
  c <- xsym[l:n]
  te <- 0
  for (i in 1:(n-l+1)) {
    #     sum of  joint(a,b,c) * log2  P(c | a,b) /  P(c | a)
    te <- te + jointProb[a[i],b[i],c[i]] * log2( safediv(condProbTop[a[i],b[i],c[i]], condProbBot[a[i],c[i]]) )
  }
  te  
}
