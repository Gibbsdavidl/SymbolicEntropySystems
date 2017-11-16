
# coupled lorenz oscillators

# two systems
# each has x,y,z
# and is directionally 
# coupled x2 -> x1 
# 

#http://www2.gsu.edu/~matixb/belykh_node_balance.pdf

last <- function(a) {
  a[length(a)]
}

dx1 <- function(sigva,x1,x2,y,ti,eps) {
  ti*(sigva*(y-x1) + eps*(x2-x1))
}

dx2 <- function(sigva,x,y,ti) {
  ti*(sigva*(y-x))
}
  
dy <- function(rho,x,y,z,ti) {
  ti*(rho*x - z*x - y)
}

dz <- function(beta, x, y, z, ti) {
  ti*(x*y - beta*z)
}


compStep <- function(vl,sigva,rho,beta,eps1,ti) {
  x1 <- last(vl[[1]]); y1 <- last(vl[[2]]); z1 <- last(vl[[3]]);
  x2 <- last(vl[[4]]); y2 <- last(vl[[5]]); z2 <- last(vl[[6]]); 
  
  dx1 <- dx1(sigva,x1,x2,y1,ti,eps1)
  dx2 <- dx2(sigva,x2,y2,ti)
  dy1 <- dy(rho,x1,y1,z1,ti)
  dy2 <- dy(rho,x2,y2,z2,ti)
  dz1 <- dz(beta,x1,y1,z1,ti)
  dz2 <- dz(beta,x2,y2,z2,ti)

  x1 <- c(vl[[1]], x1+dx1); y1 <- c(vl[[2]], y1+dy1); z1 <- c(vl[[3]], z1+dz1);
  x2 <- c(vl[[4]], x2+dx2); y2 <- c(vl[[5]], y2+dy2); z2 <- c(vl[[6]], z2+dz2);
  return(list(x1=x1,y1=y1,z1=z1,x2=x2,y2=y2,z2=z2))  
}

sim <- function(n=50, eps1=1/100) {
  # common params
  sigva <- 10
  rho <- 28
  beta <- 8/3
  timestep <- 0.005
  x1 <- c(runif(1))
  y1 <- c(runif(1))
  z1 <- c(runif(1))
  x2 <- c(runif(1))
  y2 <- c(runif(1))
  z2 <- c(runif(1))
  vl <- list(x1=x1,y1=y1,z1=z1,x2=x2,y2=y2,z2=z2)
  for (tidx in 1:n) {
    vl <- compStep(vl,sigva,rho,beta,eps1,timestep)
  }
  vl
}

# synchronized  
res0 <- sim(10000, eps1=1/2)
plot(x=1:length(res0$x1), y=res0$x1, type='l', xlab='time', ylab='x1', col='blue')
points(x=1:length(res0$x1), y=res0$x2, type='l', xlab='time', ylab='x1', col="red")

# unsynchronized
res1 <- sim(10000, eps1=0/2)
plot(x=1:length(res0$x1), y=res0$x1, type='l', xlab='time', ylab='x1', col='blue')
points(x=1:length(res0$x1), y=res0$x2, type='l', xlab='time', ylab='x1', col="red")

plot(res0$x1, res0$y1, type='l', xlab='x1', ylab='y1')
#plot(res0$z1, res0$z2, type='l', xlab='y1', ylab='y2')
#plot(res0$z1, res0$z2, type='l', xlab='z1', ylab='z2')

source("symb_ent.R")

ste(res0$x1,res0$x2,4,2)
#[1] 1.716397
ste(res0$x2,res0$x1,4,2)
#[1] -0.02907778
ste(res1$x1,res1$x2,4,2)
#[1] -0.5383814
ste(res1$x2,res1$x1,4,2)
#[1] 0.2039921


