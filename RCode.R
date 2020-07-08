#R Code

library(MASS)

#Gamma Density 
x <- seq(0, 7, by=.001) 
plot(x,  dgamma(x, 2.3, 2.7), type="l", ylim=c(0,2), ylab="Density", main="Gamma Density", sub="Gamma Distribution with Alpha=2.3 and Beta=2.7") 

# Beta Density   
x <- seq(-2, 2, length=1000)
y = dbeta(x,shape1=0.3,shape2=1)
plot(x, y, type="l", lty=1, xlab="x value", ylab="Density",ylim=c(0,2), main = "Beta Density",sub="Beta Density with alpha=0.3, beta=0.7")

#Standard Normal Distribution
x <- seq(-4, 4, length=1000)
hx <- dnorm(x)
plot(x,hx, type="l", lty=1, xlab="x value",
     ylab="Density", main="Standard Normal Distribution",sub = "Standard Normal Density with mean=0 and variance=1")

# BiVariate Normal Distribution
bivn <- mvrnorm(1000, mu = c(0, 0), Sigma = matrix(c(1, .5, .5, 1), 2))
bivn.kde <- kde2d(bivn[,1], bivn[,2], n = 50)
contour(bivn.kde)
image(bivn.kde)
persp(bivn.kde, phi = 45, theta = 30)


# Importance Sampling
w <- function(x) dunif(x,0,1)/dbeta(x,0.3,1)
f <- function(x) x^(2)
X <- rbeta(1000000,0.3,1)
Y <- w(X)*f(X)
c(mean(Y),var(Y))

# Metropolis Algorithm

norm<-function (n, alpha) 
{
  vec <- vector("numeric", n)
  x <- 0
  vec[1] <- x
  for (i in 2:n) {
    rv <- runif(1,-alpha,alpha)
    candidate <- x + rv
    accept <- min(1, dnorm(candidate)/dnorm(x))
    u <- runif(1)
    if (u < accept) 
      x <- candidate
    vec[i] <- x
  }
  vec
}
normvec<-norm(10000,1)
c(mean(normvec),var(normvec))
par(mfrow=c(2,1))
plot(ts(normvec))
hist(normvec,30)
par(mfrow=c(1,1))

# Metropolis Hastings Algorithm

gamm<-function (n, a, b) 
{
  mu <- a/b
  sig <- sqrt(a/(b * b))
  vec <- vector("numeric", n)
  x <- a/b
  vec[1] <- x
  for (i in 2:n) {
    candidate <- rnorm(1, mu, sig)
    accept <- min(1, (dgamma(candidate,a, b)/dgamma(x, a, b))/(dnorm(candidate,mu, sig)/dnorm(x, mu, sig)))
    u <- runif(1)
    if (u < accept) 
      x <- candidate
    vec[i] <- x
  }
  vec
}
vec<-gamm(100000,2.3,2.7)
c(mean(vec),var(vec))
par(mfrow=c(2,1))
plot(ts(vec))
hist(vec,30)
par(mfrow=c(1,1))

# Gibbs Sampling
gibbs<-function (n, rho) 
{
  mat <- matrix(ncol = 2, nrow = n)
  x <- 0
  y <- 0
  mat[1, ] <- c(x, y)
  for (i in 2:n) {
    x <- rnorm(1, rho * y, sqrt(1 - rho^2))
    y <- rnorm(1, rho * x, sqrt(1 - rho^2))
    mat[i, ] <- c(x, y)
  }
  mat
}
bvn<-gibbs(10000,0.97)
mean(bvn)
var(bvn)
par(mfrow=c(3,2))
plot(bvn,col=1:10000)

plot(bvn,type="l")
plot(ts(bvn[,1]), xlab='x')
plot(ts(bvn[,2]), xlab='y')
hist(bvn[,1],40, xlab='Samples', main= 'Histogram of x')
hist(bvn[,2],40, xlab= 'Samples', main= 'Histogram of y')
par(mfrow=c(1,1))

# Gibbs Sampler 

mu <- 5
sigma2 <- 2
phi <- 1/sigma2

#Draw a Sample of size 53 From a N(mu,sigma2)
x <- rnorm(53,mu,sqrt(sigma2))

#Write the data to a data file
write(x,file="./sampledata.dat",sep=",")

#Read the data file and save it as an object
mydata <- scan("./sampledata.dat",sep=",")

#Prior Values
m <- 0
s2 <- 10
a <- 2
b <- 3
n <- length(mydata)

#Allocate Space to save Gibbs sampling draws
totdraws <- 5000
mudraws <- rep(0,5000) #Initial Value of mu is 0
sigma2draws <- rep(1,5000) #Initial Value of s2 is 1

#Begin the Gibbs sampler
for(i in 2:totdraws){
  
  #Generate a Draw for mu
  mstar <- ((1/s2)*m + (n/sigma2draws[i-1])*mean(mydata))/((1/s2)+(n/sigma2draws[i-1]))
  s2star <- 1/((1/s2)+(n/sigma2draws[i-1]))
  mudraws[i] <- rnorm(1,mstar,sqrt(s2star))
  
  #Generate a Draw for Sigma2
  astar <- a + (n/2)
  bstar <- (sum((mydata-mudraws[i])^2)/2) + b
  phitemp <- rgamma(1,astar,bstar)
  sigma2draws[i] <- 1/phitemp
  
}

#Discard the first 500 Draws as Burn-in
mudraws <- mudraws[501:length(mudraws)]
sigma2draws <- sigma2draws[501:length(mudraws)]

#Partition plot into 4 subplots
par(mfrow=c(2,2))

#Plot Trace Plots
plot(mudraws,type="l")
plot(sigma2draws,type="l")

#Plot Marginal Posterior Densities
plot(density(mudraws),col="red",main=expression(paste("Marginal Posterior Density of ",mu)))
plot(density(sigma2draws),col="blue",main=expression(paste("Marginal Posterior Density of ",sigma^2)))

#Shooting Percentage using Metropolis
x = rnorm(1000,3.24,0.33);
T <- 10000; B <- 1000; u.mu <- runif(T); u.sig <- runif(T)
mu <- sig <- numeric(T); mu[1] <- 3; sig[1] <- 5 
REJmu <- 0; REJsig <- 0 
logpost = function(x,mu,sig){ 
  loglike = sum(dnorm(x,mu,sig,log=TRUE))
  return(loglike - log(sig))}
for (t in 2:T) { mut <- mu[t-1]; sigt <- sig[t-1] 
mucand <- mut + runif(1,-0.5,0.5) 
sigcand <- abs(sigt + runif(1,-0.5,0.5)) 
log.alph.mu = logpost(x,mucand,sigt)-logpost(x,mut,sigt) 
if (log(u.mu[t]) < log.alph.mu) mu[t] <- mucand 
else { mu[t] <- mut; REJmu <- REJmu+1 } 
log.alph.sig = logpost(x,mu[t],sigcand)-logpost(x,mu[t],sigt) 
if (log(u.sig[t]) < log.alph.sig) sig[t] <- sigcand 
else { sig[t] <- sigt; REJsig <- REJsig+1 }} 
REJratemu <- REJmu/T; REJratesig <- REJsig/T; 
samp <- 1001:10000; musamp <- mu[samp]; sigsamp <- sig[samp] 
hist(musamp,50, xlab='Shooting Percentages', main='Historgram of Shooting Percentages') 
mu.den <- density(musamp)
plot(mu.den)  
plot(musamp,sigsamp) 



