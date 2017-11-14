
source("symb_ent.R")

x <- runif(100)
y <- runif(100)
ste(x,y,2,2)
#system.time( ste(y,x,2,2) )
#user  system elapsed 
#0.005   0.000   0.005 

# simulation from: https://www.ncbi.nlm.nih.gov/pubmed/22500692
# Transfer entropy estimation and directional coupling change detection in biomedical time series.

library(NormalLaplace)
y1 <- rnorm(mean=10, sd=1, n=100)
x1 <- sapply(3:100, function(i) (0.7*y[i-2])^2)
y <- y1 + rnl(n=100)
x <- x1 + rnl(n=98)
y <- y[3:100]
ste(y,x,2,2)
ste(x,y,2,2)

# trying different time lags
ste(y,x,2,1) ### think that lag 1 is really like lag 0 somehow... bug here?
#[1] 0
ste(y,x,2,2)  ### this is the 'real' time lag
#[1] 1.344819
ste(y,x,2,3)
#[1] 1.232327
ste(y,x,2,4)
#[1] 0.5231706
ste(y,x,2,5)
#[1] 0.69572

# example from: http://www.sciencedirect.com/science/article/pii/S0167278902004323
#Information transfer in continuous processes

a <- 0.5
b <- 0.5
g <- 0.5
x <- rep(0.5, 10000)
for (i in 1:99) {
  x[i+1] <- a*x[i] + rnorm(n=1, sd=0.1)
}
y <- rep(0.6, 10000)
for (i in 1:999) {
  y[i+1] <- b*y[i] + g*x[i] + rnorm(n=1, sd=0.1)
}
# here x is influencing y
# so expect that ste(x,y) greater  [x -> y]
# because we are predicting y into the future using current x,y

ste(y,x,3,2)
#[1] 0.009677885

ste(x,y,3,2)
#[1] 0.1796151

# and it is... usually.

