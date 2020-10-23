
# From the lectures

gibbsstep <- function(mu, alpha, ybar, J, 
                      sigmaeps = 2, sigmaalph = 4) {
  p <- length(alpha)
  n <- sum(J)
  ymean <- sum(J * ybar) / n
  V <- 1 / (J / sigmaeps^2 + 1 / sigmaalph^2)
  mu <- rnorm(1, ymean - sum(J * alpha) / n, sigmaeps / sqrt(n))
  alpha <- rnorm(p, J * V * (ybar - mu) / sigmaeps^2, sqrt(V))
  list(mu = mu, alpha = alpha)
}

#------------------------------------------------------------------------------------------------------------------------------------------

# I have sokething like this
gibbs <- function(mu0, alpha0, data, t = 1000, ...) {
  mu <- numeric(t)
  alpha <- matrix(0, length(alpha0), t)
  mu[1] <- mu0
  alpha[, 1] <- alpha0
  J <- table(data$lab)
  ybar <- tapply(data$y, data$lab, mean)
  for(i in 2:t) {
    tmp <- gibbsstep(mu[i - 1], alpha[, i - 1], ybar, J, ...)
    mu[i] <- tmp$mu
    alpha[, i] <- tmp$alpha
  }
  list(mu = mu, alpha = alpha) 
}


gibbsfast <- function(mu0, alpha0, data, t = 1000, sigmaeps = 2, 
                      sigmaalph = 4, ...) {
  mu <- numeric(t)
  alpha <- matrix(0, length(alpha0), t)
  mu[1] <- mu0
  alpha[, 1] <- alpha0
  J <- table(data$lab)
  p <- length(alpha0)
  n <- sum(J)
  ybar <- tapply(data$y, data$lab, mean)
  ymean <- sum(J * ybar) / n
  V  <-  1 / (J / sigmaeps^2 + 1 / sigmaalph^2)
  sigmaepsn <-  sigmaeps / sqrt(n)
  sqrtV <- sqrt(V)
  JVsig <- J * V / sigmaeps^2
  for(i in 2:t) {
    mu[i] <- rnorm(1, ymean - sum(J * alpha[, i - 1]) / n, sigmaepsn)
    alpha[, i] <- rnorm(p, JVsig * (ybar - mu[i]), sqrtV)
  }
  list(mu = mu, alpha = alpha) 
}
########################################################################################################################################
#########################################################################################################################################x

gibbsfastest <- function(mu0, alpha0, data, t = 1000, sigmaeps = 2, 
                         sigmaalph = 4, ...) {
  mu <- numeric(t)
  alpha <- matrix(0, length(alpha0), t)
  mu[1] <- mu0
  alpha[, 1] <- alpha0
  J <- table(data$lab)
  
  # it changes here
  p <- length(alpha0)
  q<-length(beta0)
  n <- sum(J)
  ybar <- tapply(data$y, data$lab, mean)
  ymean <- sum(J * ybar) / n
  
  #it calculates the V outside of the step function
  V  <-  1 / (J / sigmaeps^2 + 1 / sigmaalph^2)
  
  sigmaepsn <-  sigmaeps / sqrt(n)
  sqrtV <- sqrt(V)
  JVsig <- J * V / sigmaeps^2
  J <- J / n
  W <- rnorm(t, ymean, sigmaepsn)
  
  U <- matrix(rnorm(t * p, JVsig * ybar, sqrtV), p, t)
  
  for(i in 2:t) {
    mu[i] <- W[i] - sum(J * alpha[, i - 1])
    alpha[, i] <- U[, i] - JVsig * mu[i]
  }
  list(mu = mu, alpha = alpha) 
}
#==================================================================================================================================================

gibbs <- function(mu0, alpha0, data, t = 1000, ...) {
  mu <- numeric(t)
  alpha <- matrix(0, length(alpha0), t)
  mu[1] <- mu0
  alpha[, 1] <- alpha0
  J <- table(data$lab)
  
  
  ybar <- tapply(data$y, data$lab, mean)
  for(i in 2:t) {
    tmp <- gibbsstep(mu[i - 1], alpha[, i - 1], ybar, J, ...)
    mu[i] <- tmp$mu
    alpha[, i] <- tmp$alpha
  }
  list(mu = mu, alpha = alpha) 
}


##################################### My implement - original ###################################################################################xxxxxxx

gibbsstep <- function(mu, alpha,betaa, ybar, J,y,
                      sigmaeps = 1, sigmaalph = sqrt(86),sigmabeta=sqrt(58)) {
  p <- length(alpha) 
  n <- sum(J)
  q<-length(betaa)
  ymean = sum(J * ybar) / n # ybar is the average of the groups
  V1 <- 1 / (J / sigmaeps^2 + 1 / sigmaalph^2)
  V2 <- 1/  (1 / sigmaeps^2+ 1 / sigmabeta^2)
  
  
  group<-rep(1:15, each=2)
 betaa2 <- as.vector(rowsum(betaa, group))
  #betaa2 <- colSums(matrix(beta,nrow=2)) #improved version
  
  mu <- rnorm(1, ymean - sum(J * alpha) / n-sum(betaa)/n, sigmaeps / sqrt(n)) 
  
  alpha <- rnorm(p, J * V1 * ((ybar - mu)-1/J*(betaa2)) / sigmaeps^2, sqrt(V1)) 
  
  betaa<-rnorm(q,V2/sigmaeps^2*(y-mu-rep(alpha,each=2)),sqrt(V2)) 
  
  list(mu = mu, alpha = alpha,betaa=betaa) 
}

### ------------------------ IMPROVED VERSION -----------------------------------------------------------------------------------------------------
profvis({
gibbsimproved <- function(mu0, alpha0,beta0, dataa, t, sigmaeps = 1, sigmaalph = sqrt(86),sigmabeta=sqrt(58)) { 
  mu <- numeric(t)
  
  alpha <- matrix(0, length(alpha0), t)
  betaa<-matrix(0, length(beta0), t)
  mu[1] <- mu0 
  alpha[, 1] <- alpha0 
  betaa[,1]<-beta0
  J <- table(dataa$Batch) 
  y<-dataa$Moisture
  
  p <- length(alpha0)
  q<-length(beta0)
  n <- sum(J)
  ybar <- tapply(dataa$Moisture, dataa$Batch, mean)
  ymean = sum(J * ybar) / n # ybar is the average of the groups
 
  # variances, now we only calculate them once
  V1 <- 1 / (J / sigmaeps^2 + 1 / sigmaalph^2)
  V2 <- 1/  (1 / sigmaeps^2+ 1 / sigmabeta^2)
  
  #taking square roots for std
  sigmaepsn <-  sigmaeps / sqrt(n)
  sqrtV1 <- sqrt(V1)
  sqrtV2 <- sqrt(V2)
  
  JVsig <- J * V1 / sigmaeps^2
  Vsigma<-V2/sigmaeps^2
  Jinv<-1/J # worth mentioning!!
  
  # Generating the random numbers
  W <- rnorm(t, ymean, sigmaepsn)
  
  U <- matrix(rnorm(t * p, JVsig * ybar, sqrtV1), p, t)
  
  Z<-matrix(rnorm(t * q,  Vsigma* y, sqrtV2), q, t)
  
  # actual step process
  for(i in 2:t) {
    
    # Keep the previous lines for PROFILING !!
    #group<-rep(1:15, each=2)
    #betaa2 <- as.vector(rowsum(betaa[, i-1], group)) # it will need improvements, colsum maybe?
    betaa2 <- colSums(matrix(betaa[, i-1],nrow=2))
    
    mu[i] <- W[i] - (1/n)*sum(J * alpha[, i - 1])-sum(betaa[, i-1])/n
    alpha[, i] <- U[,i] - JVsig * mu[i]-JVsig*Jinv*betaa2 # adding the Jinv here also improved the run time
    betaa[, i]<-Z[,i]-Vsigma*mu[i]-Vsigma*alpha[, i]
  }
  
  list(mu = mu, alpha = alpha, betaa=betaa) 
}
#---------------- End of function ------------------------------------------------------------------------------------------------------------------

outputt<-gibbsimproved(0,rep(0, 15),rep(0, 30),pigment,t=10000)
})

# end of profiling

# Sample paths
plot(seq(1,10000,1),outputt$betaa[13,],type = "l")
plot(seq(1,10000,1),outputt$betaa[2,],type = "l")

plot(seq(1,10000,1),outputt$alpha[8,],type = "l")

# Posterior distribution
Y1<-outputt$betaa[1,2:10000]+outputt$mu[2:10000]+outputt$alpha[1,2:10000] #burn in the 1st period, cause the initial values are 0s
plot(density(Y1),lwd=2)
abline(v = mean(Y1),col="Red",lwd=4)

# autocorrelations
acf(outputt$mu, lwd = 2,lag.max = 10000) #lwd is just 


# BYTE compilation #===================================================================================================================
gibbsbyte<-cmpfun(gibbsimproved)

gibssbenchmark<-microbenchmark(
  gibbsimproved(0,rep(0, 15),rep(0, 30),pigment,t=10000),
  gibbsbyte(0,rep(0,15),rep(0, 30),pigment,t=10000)
)
gibssbenchmark # it can be mentioned , very little, almost no improvements
