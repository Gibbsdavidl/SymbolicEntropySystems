
# coupled lorenz oscillators

# two systems
# each has x,y,z
# and is coupled 
# by an external variable w


last <- function(a) {
  a[length(a)]
}

dx <- function(sigva,x,y,w,ti,eps,mu) {
  ti*(sigva*(y-x) + eps*mu*w)
}

dy <- function(rho,x,y,z,ti) {
  ti*((rho-z)*x - y)
}

dz <- function(beta, x, y, z, ti) {
  ti*(x*y - beta*z)
}

dw <- function(kappa, y, eps2, mu1, mu2, x1, x2, ti) {
  ti*(-kappa*y - (eps2/2)*(mu1*x1 + mu2*x2))
}

compStep <- function(vl,sigva,rho,beta,kappa,mu1,mu2,eps1,eps2,ti) {
  x1 <- last(vl[[1]]); y1 <- last(vl[[2]]); z1 <- last(vl[[3]]);
  x2 <- last(vl[[4]]); y2 <- last(vl[[5]]); z2 <- last(vl[[6]]); w <- last(vl[[7]]);
  
  dx1 <- dx(sigva,x1,y1,w,ti,eps1,mu1)
  dx2 <- dx(sigva,x2,y2,w,ti,eps1,mu2)
  dy1 <- dy(rho,x1,y1,z1,ti)
  dy2 <- dy(rho,x2,y2,z2,ti)
  dz1 <- dz(beta,x1,y1,z1,ti)
  dz2 <- dz(beta,x2,y2,z2,ti)
  dw  <- dw(kappa,w,eps2,mu1,mu2,x1,x2,ti)
  
  x1 <- c(vl[[1]], x1+dx1); y1 <- c(vl[[2]], y1+dy1); z1 <- c(vl[[3]], z1+dz1);
  x2 <- c(vl[[4]], x2+dx2); y2 <- c(vl[[5]], y2+dy2); z2 <- c(vl[[6]], z2+dz2);
  w <- c(vl[[7]], w + dw)
  return(list(x1=x1,y1=y1,z1=z1,x2=x2,y2=y2,z2=z2,w=w))  
}

sim <- function(n=50, eps1=1/100, eps2=1/100) {
  # common params
  sigva <- 10
  rho <- 28
  beta <- 8/3
  mu1 <- -1
  mu2 <- 1
  kappa <- 1/1
  timestep <- 0.005
  x1 <- c(runif(1))
  y1 <- c(runif(1))
  z1 <- c(runif(1))
  x2 <- c(runif(1))
  y2 <- c(runif(1))
  z2 <- c(runif(1))
  w  <- c(0.1)
  vl <- list(x1=x1,y1=y1,z1=z1,x2=x2,y2=y2,z2=z2,w=w)
  for (tidx in 1:n) {
    vl <- compStep(vl,sigva,rho,beta,kappa,mu1,mu2,eps1,eps2,timestep)
  }
  vl
}
  
res0 <- sim(10000, eps1=1/2, eps2=1/1)
plot(x=1:length(res0$x1), y=res0$x1, type='l', xlab='time', ylab='x1', col='blue')
points(x=1:length(res0$x1), y=res0$x2, type='l', xlab='time', ylab='x1', col="red")


#plot(res0$x1, res0$y1, type='l', xlab='x1', ylab='y1')
#plot(res0$z1, res0$z2, type='l', xlab='y1', ylab='y2')
#plot(res0$z1, res0$z2, type='l', xlab='z1', ylab='z2')
