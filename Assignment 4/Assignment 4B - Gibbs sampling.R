
# comparison between parameters
# "effective sample size"

#-- JENS Code:::

# install.packages("profvis")
install.packages("sna")
library(profvis)
library(microbenchmark)
library(ggplot2)
library(Rfast) # for colMins
library(magic)
library(sna) # for plot.matrix -> sociomatrix

# setwd("C:/Users/Peti/Desktop/KU/Computational statistics/CSwR-master/CSwR-master/data")
setwd("C:/Users/Jens/Desktop/KU/Computational Statistics")
pigment <- read.table("pigment.txt",header=TRUE)


ybar_ <- tapply(pigment$Moisture, pigment$Batch, mean)
J_ <- table(pigment$Batch)
y_ <- pigment$Moisture


# 
betaa2 <- colSums(matrix(betaa0,nrow=2))

# profiling
profvis({
gibbsstep <- function(mu0, alpha0, betaa0, gamma0, eta0, ybar, J, y,
                      sigmaeps=1, sigmaalph=sqrt(86),sigmabeta=sqrt(58)) {
  p <- length(alpha0)
  n <- sum(J)
  q <- length(betaa0)
  V1 <- 1 / (J / sigmaeps^2 + 1 / sigmaalph^2)
  V2 <- 1/  (1 / sigmaeps^2 + 1 / sigmabeta^2)
  group <- rep(1:15, each=2)
  betaa2 <- as.vector(rowsum(betaa0, group))
  
  mu <- rnorm(1, mean(y) - mean(alpha0) - mean(betaa0), sigmaeps/sqrt(n) )
  alpha <- rnorm(p, J*V1/(sigmaeps^2) * ((ybar - mu)-1/J*(betaa2)), sqrt(V1))
  betaa <- rnorm(q,   V2/(sigmaeps^2) * (y-mu-rep(alpha,each=2)), sqrt(V2))
  eta <- mu + rep(alpha,each=2) + betaa
  gammaa <- mu + alpha
    
  list(mu=mu, alpha=alpha, betaa=betaa, eta=eta, gammaa=gammaa) 
}
for(i in 1:1e5) g <- gibbsstep(0,rep(1,15),rep(3,30), ybar_,J_,y_)
})

gibbs <- function(mu0, alpha0, betaa0, data, t, burnin=0, ...) { 
  J <- table(data$Batch)
  y <- data$Moisture
  ybar <- tapply(data$Moisture, data$Batch, mean) #group mean
  
  tmp <- list(mu=mu0, alpha=alpha0, betaa=betaa0)
  for(i in 1:burnin)
    tmp <- gibbsstep(tmp$mu, tmp$alpha, tmp$betaa, ybar, J,y, ...)
  
  mu <- numeric(t)
  alpha <- matrix(0, length(alpha0), t)
  betaa <- matrix(0, length(betaa0), t)
  gammaa <- alpha
  eta <- betaa
  for(i in 1:t) {
    tmp <- gibbsstep(tmp$mu, tmp$alpha, tmp$betaa, ybar, J,y, ...)
    mu[i] <- tmp$mu
    alpha[,i] <- tmp$alpha
    betaa[,i]<- tmp$betaa
    gammaa[,i]<- tmp$gammaa
    eta[,i]<- tmp$eta
  }
  list(mu=mu, alpha=alpha, betaa=betaa, gammaa=gammaa, eta=eta, x=seq(1,t))
}
#==========================================================================================
# end of function

#gibbs(0,rep(0, 15),rep(0, 30),pigment,t=10000) 
#})

tmax <- 1000
set.seed(1234567)
output <- gibbs(0, rep(0,15), rep(0,30), pigment, t=tmax, burnin=3) 

# output

# checking that sections of simulated streams are similar -> stable distribution
summary(output$mu[(tmax/3):(tmax*2/3)])
summary(output$mu[(tmax*2/3):tmax])

b1 <- cor(t(output$betaa[,(tmax*.2):(tmax*.6)]),t(output$betaa[,(tmax*.2):(tmax*.6)]))
b2 <- cor(t(output$betaa[,(tmax*.6):(tmax*1 )]),t(output$betaa[,(tmax*.6):(tmax*1)]))



#plotting the sample paths of the different parameters
plot(seq(1,tmax,1),output$mu,type="l")
plot(seq(1,10000),output$mu[1:10000],type="l")
plot(seq(50000,50100),output$mu[50000:50100],type="l")

plot(seq(1,tmax,1),output$betaa[1,],type = "l")
points(seq(1,tmax,1), output$betaa[1,], type='l',col="Blue")
points(seq(1,tmax,1), output$betaa[4,], type='l',col="Red")

#alpha values
points(seq(1,tmax,1), output$alpha[1,], type='l',col="Blue")
points(seq(1,tmax,1), output$alpha[4,], type='l',col="Red")
points(seq(1,tmax,1), output$alpha[5,], type='l',col="Green")

#autocorrelation
acf(output$mu, lwd = 2, lag.max=1000) #lwd is just 
acf(output$alpha[3,], lwd = 2,lag.max = 1000) #most of the correlations are not insignificant

acf(output$betaa[2,], lwd = 2,lag.max = 1000)
#what can be seen from these autocorrelation functions??

acf(cbind(output$mu,output$alpha[1,]), lwd = 2,lag.max = 1000)

# correlation matrix
cormat <- function(output,corfrom=1,cornum=1000){
  corx <- rbind(output$mu,output$alpha,output$betaa)[,corfrom:(corfrom+cornum-1)]
  cor(t(corx))
}
c1000 <- cormat(output,corfrom=1000)
c10000 <- cormat(output,corfrom=10000)
c100000 <- cormat(output,corfrom=100000)
plot(density(abs(c(c1000))))
lines(density(abs(c(c10000))))
lines(density(abs(c(c100000))))

c1 <- cormat(output)
c100 <- cormat(output,corfrom=100)
c1000 <- cormat(output,corfrom=1000)
c10000 <- cormat(output,corfrom=10000)
c100000 <- cormat(output,corfrom=100000)

# plot(density(abs(c(c1-c100))))
# lines(density(abs(c(c100-c1000))))
plot(density(abs(c(c1000-c10000))))
lines(density(abs(c(c10000-c100000))))

summary(abs(c(c1-c100)))
summary(abs(c(c100-c1000)))
summary(abs(c(c1000-c10000)))
summary(abs(c(c10000-c100000)))

plot(density(abs(c(c1))))
plot(density(abs(c(c100))))
lines(density(abs(c(c1000))))
lines(density(abs(c(c10000))))
lines(density(abs(c(c100000))))

plot(density(abs(c(c1000))))
lines(density(abs(c(c10000))))
lines(density(abs(c(c100000))))

corfrom <- 1
cornum <- 100
corx <- rbind(output$mu,output$alpha,output$betaa)[,corfrom:(corfrom+cornum-1)]
dim(corx)
cormat <- cor(t(corx))

u <- cormat[upper.tri(cormat)]
dim(u)
?upper



#                                    Part 2

#########-----Hierarchical Gibbs---------------------------------------------------------------------------------------


gibbsstepHir <- function(mu0, alpha0, beta0, gammaa0, eta0, J,y,
                          sigmaeps=1, sigmaalph=sqrt(86), sigmabeta=sqrt(58)) {
  I <- length(gammaa0)
  V2 <- 1/ (1/ sigmaeps^2 + 1/sigmabeta^2)
  V3 <- 1/ (J/sigmabeta^2 + 1/sigmaalph^2)
  eta0sum <- colSums(matrix(eta0,nrow=2)) # pairwise sum

  mu     <- rnorm(1, 1/I*sum(gammaa0), sigmaalph/sqrt(I))
  gammaa <- rnorm(I, V3*( (1/sigmabeta^2)*eta0sum + mu/sigmaalph^2), sqrt(V3)) # found it!!!
  eta    <- rnorm(length(eta0), V2*(y/sigmaeps^2 +rep(gammaa,each=2)/sigmabeta^2), sqrt(V2))
  
  alpha <- gammaa - mu
  beta <- eta - rep(gammaa, each=2)
  
  list(mu = mu, gammaa=gammaa, eta=eta, alpha=alpha, beta=beta)
}

h <- gibbsstepHir(0, rep(0,15), rep(0,30), J_,y_)


#-------------------------------------------------------------------------------------

gibbsHir <- function(mu0, gammaa0, eta0, data, t, burnin=0, ...) {
  
  J <- table(data$Batch) # add our own data here, that how many samples are in a batch
  y <- data$Moisture # the individual values

  tmp <- list(mu=mu0, gammaa=gammaa0, eta=eta0)
  for(i in 1:burnin)
    tmp <- gibbsstepHir(tmp$mu, 0,0, tmp$gammaa, tmp$eta, y,J, ...)

  mu <- numeric(t) # making a t length vector of mu
  gammaa <- matrix(0, length(gammaa0), t)
  eta <- matrix(0, length(eta0), t)
  alpha <- matrix(0, length(gammaa0), t)
  beta <- matrix(0, length(eta0), t)

  for(i in 1:t) {
    tmp <- gibbsstepHir(tmp$mu, 0,0, tmp$gammaa, tmp$eta, y,J, ...)
    mu[i] <- tmp$mu
    gammaa[,i] <- tmp$gammaa
    eta[,i] <- tmp$eta
    alpha[,i] <- tmp$alpha
    beta[,i] <- tmp$beta
  }
  
  list(mu=mu, gammaa=gammaa, eta=eta, alpha=alpha, beta=beta, x=seq(1,t))
}

#-----------------------------------------------------------------------------------------------------------------------------
# END of function

tmax <- 1000000
set.seed(1234567)
out <- gibbs(0, rep(0,15), rep(0,30), pigment, t=tmax, burnin=5)
set.seed(1234567)
outHir <- gibbsHir(0, rep(0,15), rep(0,30), pigment, t=tmax, burnin=5)

plot(y=out$mu[5e5:1e6],x=out$x[5e5:1e6], type="l")
plot(y=out$mu[9e5:1e6],x=out$x[9e5:1e6], type="l")
plot(y=out$mu[99e4:1e6],x=out$x[99e4:1e6], type="l")
plot(y=out$mu[999e3:1e6],x=out$x[999e3:1e6], type="l")
plot(y=out$alpha[1,999e3:1e6],x=out$x[999e3:1e6], type="l")
plot(y=out$beta[1,999e3:1e6],x=out$x[999e3:1e6], type="l")

out$my <- out$mu
out$my[999e3] <- -5
plot(y=out$my[999e3:1e6],x=out$x[999e3:1e6], type="l", col="red", main="Gibbs sampler: Sample paths\n(Mu is red, alpha is green, beta is blue)", ylab="Realizations")
lines(y=out$alpha[1,999e3:1e6],x=out$x[999e3:1e6], type="l", col="green")
lines(y=out$beta[1,999e3:1e6],x=out$x[999e3:1e6], type="l", col="blue")

outHir$my <- outHir$mu
outHir$my[999e3] <- -5
plot(y=outHir$my[999e3:1e6],x=outHir$x[999e3:1e6], type="l", col="red", main="Hierarchical Gibbs sampler: Sample paths\n(Mu is red, alpha is green, beta is blue)", ylab="Realizations")
lines(y=outHir$alpha[1,999e3:1e6],x=outHir$x[999e3:1e6], type="l", col="green")
lines(y=outHir$beta[1,999e3:1e6],x=outHir$x[999e3:1e6], type="l", col="blue")

?plot
plot(y=outHir$mu[999e3:1e6],x=outHir$x[999e3:1e6], type="l")
plot(y=outHir$alpha[1,999e3:1e6],x=outHir$x[999e3:1e6], type="l")
plot(y=outHir$beta[1,999e3:1e6],x=outHir$x[999e3:1e6], type="l")


acf(out$mu[5e5:1e6], lwd = 2,lag.max = 1000)
acf(out$mu[5e5:1e6], lwd = 2,lag.max = 10000)
acf(out$mu[5e5:1e6], lwd = 2,lag.max = 100000)
acf(out$alpha[1,5e5:1e6], lwd = 2,lag.max = 1000)
acf(out$beta[1,5e5:1e6], lwd = 2,lag.max = 1000)
acf(out$beta[1,5e5:1e6], lwd = 2,lag.max = 10000)

acf(outHir$mu, lwd = 2,lag.max = 10)
acf(outHir$alpha[1,], lwd = 2,lag.max = 10)
acf(outHir$beta[1,], lwd = 2,lag.max = 10)
acf(outHir$eta[1,], lwd = 2,lag.max = 10)
acf(outHir$gammaa[1,], lwd = 2,lag.max = 10)

acf(cbind(out$mu,out$alpha[1,]), lwd = 2,lag.max = 1000)
acf(cbind(out$mu,out$beta[1,]), lwd = 2,lag.max = 1000) # smaller
acf(cbind(out$alpha[1,],out$beta[1,]), lwd = 2,lag.max = 1000) # shorter
?plot

### Check if they draw samples from the same distributions? ###

plot(y=outHir$alpha[1,999e3:1e6],x=outHir$x[999e3:1e6], type="l", col="red", main="10000 alpha draws. Base is blue, Hierarchical is red.")
lines(y=out$alpha[1,999e3:1e6],x=out$x[999e3:1e6], type="l", col="blue")

plot(y=outHir$mu[999e3:1e6]+6.8,x=outHir$x[999e3:1e6], type="l", col="red", main="10000 mu draws. Base is blue, Hierarchical is red.")
lines(y=out$mu[999e3:1e6],x=out$x[999e3:1e6], type="l", col="blue")

plot(y=outHir$beta[999e3:1e6]-6.7,x=outHir$x[999e3:1e6], type="l", col="red", main="10000 beta draws. Base is blue, Hierarchical is red.")
lines(y=out$beta[999e3:1e6],x=out$x[999e3:1e6], type="l", col="blue")

summary(out$mu[1,])
summary(outHir$mu[1,]) # very similar
plot(density(out$mu[5e5:1e6]), col="blue", main="Density of 500k mu draws\n(Base is blue, Hierarchical is red)")
lines(density(outHir$mu[5e5:1e6]), col="red")

summary(out$alpha[1,])
summary(outHir$alpha[1,]) # very similar
plot(density(out$alpha[1,5e5:1e6]), col="blue", main="Density of 500k alpha draws\n(Base is blue, Hierarchical is red)")
lines(density(outHir$alpha[1,5e5:1e6]), col="red")

summary(out$beta[1,])
summary(outHir$beta[1,]) # very similar
plot(density(out$beta[1,5e5:1e6]), col="blue", main="Density of 500k beta draws\n(Base is blue, Hierarchical is red)")
lines(density(outHir$beta[1,5e5:1e6]), col="red")

c <- cov(t(rbind(t(as.matrix(drop(out$mu))),out$alpha,out$beta)[,5e5:1e6]))
plot.sociomatrix(c, diaglab=FALSE, main="Covariance matrix of Gibbs samples")

chir <- cov(t(rbind(t(as.matrix(drop(outHir$mu))),outHir$alpha,outHir$beta)[,5e5:1e6]))
plot.sociomatrix(chir, diaglab=FALSE, main="Covariance matrix of Hierarchical Gibbs")



cumStdErr <- function(v){
  std <- numeric(length(v))
  cumsum <- 0
  cum2sum <- 0
  for(i in 1:length(v)){
    cumsum <- cumsum + v[i]
    cum2sum <- cum2sum + (v[i]-cumsum/i)^2
    std[i] <- sqrt(cum2sum/i)
  }
  std
}
addStdErrs <- function(obj,discard=0){
  f <- discard+1
  t <- length(obj$mu)
  list(mu=obj$mu[f:t], alpha=obj$alpha[f:t], beta=obj$alpha[f:t], x=obj$x[f:t],
       mustd=cumStdErr(obj$mu)[f:t], alphastd=cumStdErr(obj$alpha[1,])[f:t], betastd=cumStdErr(obj$beta[1,])[f:t]
  )
}
out1 <- addStdErrs(out, discard=0)
out2 <- addStdErrs(outHir, discard=0)

plot(y=out1$mustd[1e4:1e6],x=out1$x[1e4:1e6], type="l", col="blue")
lines(y=out2$mustd[1e1:1e6],x=out2$x[1e1:1e6], type="l", col="red")

plot(y=out2$alphastd[1e1:1e6],x=seq(1e1:1e6), type="l", col="red")
lines(y=out1$alphastd[1e4:1e6],x=seq(1e4:1e6), type="l", col="blue")
summary(out1$alphastd)
summary(out2$alphastd)

plot(y=out1$betastd[1e4:1e6],x=seq(1e4:1e6), type="l", col="blue")
lines(y=out2$betastd[1e1:1e6],x=seq(1e1:1e6), type="l", col="red")
summary(out1$betastd)
summary(out2$betastd)


cumStd.my1 <- cumStdErr(out$mu)
plot(y=cumStd.my1[1e3:1e6],x=outHir$x[1e3:1e6], type="l")


#  Auto correlation functions

acf(output2$mu, lwd = 2,lag.max = 10)
acf(output2$gammaa[1,], lwd = 2,lag.max = 10)
acf(output2$gammaa[2,], lwd = 2,lag.max = 10)
acf(output2$gammaa[3,], lwd = 2,lag.max = 10)
acf(output2$gammaa[4,], lwd = 2,lag.max = 10)
acf(output2
    $gammaa[5,], lwd = 2,lag.max = 10)
acf(output2$eta[1,], lwd = 2,lag.max = 10)


dim(output2$mu)
dim(output2$gammaa)
dim(output2$eta)
dim(output2$alpha)
acf(output2$alpha[1,], lwd = 2,lag.max = 10)
b <- output2$eta[1:2,] - output2$gammaa[1,]
dim(b)
acf(b[1,], lwd = 2,lag.max = 10)

acf(cbind(output2$mu,output2$eta[1,]), lwd = 2,lag.max = 10)
acf(cbind(output2$mu,output2$gammaa[1,]), lwd = 2,lag.max = 10)
acf(cbind(output2$eta[1,],output2$gammaa[1,]), lwd = 2,lag.max = 10)

cor(output2$mu,output2$gammaa[1,])

#plotting sample paths of eta
plot(seq(1,tmax,1), output2$eta[1,], type='l',col="Blue")
plot(seq(1,tmax,1), output2$eta[2,], type='l',col="Blue")
plot(seq(1,tmax,1), output2$eta[4,], type='l',col="Red")
#this cases there is clearly a burn in period, the samples from the beginning of the chain,
#not representing the well the distribution of the whole !!!

plot(seq(1,tmax,1), output2$gammaa[1,], type='l',col="Blue")
plot(seq(1,tmax,1), output2$gammaa[2,], type='l',col="Blue")

#=====================================================================================================

points(seq(1,tmax,1), output$betaa[5,], type='l',col="Green")

# Correlation matrix
cormat2 <- function(output,corfrom=1,cornum=100){
  corx <- rbind(output$mu,output$gammaa,output$eta)[,corfrom:(corfrom+cornum-1)]
  cor(t(corx))
}
c1000 <- cormat2(output2,corfrom=1000)
c10000 <- cormat2(output2,corfrom=10000)
c100000 <- cormat2(output2,corfrom=100000)
plot(density(abs(c(c100000))))
lines(density(abs(c(c1000))))
lines(density(abs(c(c10000))))

c1 <- cormat2(output2)
c100 <- cormat2(output2,corfrom=100)
c1000 <- cormat2(output2,corfrom=1000)
c10000 <- cormat2(output2,corfrom=10000)
c100000 <- cormat2(output2,corfrom=90000)

# Density

# plot(density(abs(c(c1-c100))))
# lines(density(abs(c(c100-c1000))))
plot(density(abs(c(c1000-c10000))))
lines(density(abs(c(c10000-c100000))))

plot(density(c(c100)))
lines(density(c(c1000)))
lines(density(c(c10000)))
lines(density(c(c100000)))
