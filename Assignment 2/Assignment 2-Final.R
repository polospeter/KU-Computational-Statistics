#Assignment 2-----------------------------------------------------------------------------------------------------
  
  # a rejection sampling can be implemented without knowing the normalizing constant
  
  # Rejection sampling
  ##================================= 1st version =====================================================================================================

vMsim <- function(n) {
  y <- numeric(n)
  count <- 0
  for(i in 1:n) {
    ratio <- 0 
    u <- 1
    
    while(u > ratio) {
     
      y0 <- qnorm(runif(1),mean = 1/sqrt(3)-0.03,sd=0.54) 
     
      ratio <-1/2.5*exp(-y0^3+y0)*(y0>=0)/dnorm(y0,mean = 1/sqrt(3)-0.03,sd=0.54)
      u <- runif(1)
      count <- count + 1
    }
    y[i] <- y0
    
  }
  y
}

#----------------------------------------------------------------------------------------------------------------------------------------------------
# general solution- making the envelope an input parameter

vMsim.input <- function(n,func, envelop,env.rand,trace=FALSE) {
  y <- numeric(n)
  count <- 0
  
  for(i in 1:n) {
    ratio <- 0 
    u <- 1
    
    while(u > ratio) {
      y0 <- env.rand(1) 
      ratio<-func(y0)/envelop(y0)
      u <- runif(1)
      count <- count + 1
    }
    y[i] <- y0
    
  }
  y
  
  if(trace)
    cat((count-n)/count) # rejection frequency
}

#==================================================================================================================

#count process
count <- count + 1
c<-(count-n)/count
list(y=y,c=c)

#============================= 2ND VERSION - VECTORIZED FORM =================================================================
#define the function and the envelope outside

func<-function(x) { (x>=0)*exp(-x^3+x) } 
envelop<-function(x){2.5*dnorm(x,mean = 1/sqrt(3)-0.03,sd=0.54)} 
gen.envelop<-function(x){qnorm(runif(x),mean = 1/sqrt(3)-0.03,sd=0.54)}
env.rand<-function(x){qnorm(runif(x),mean = 1/sqrt(3)-0.03,sd=0.54)} #i think this is env.rand

vMsim.vectorized <- function(n,func,envelop,env.rand) {
  m<-n
  fact <- 1
  y <- numeric(n)
  u <- runif(m)
  y0 <-env.rand(m)
  j<-1
  for(i in 1:n) {
    repeat{
      z <- y0[j]
      ratio <-func(z)/envelop(z) 
      accept <- u[j] <=ratio
      j <- j + 1
      if(j > m) {
        if(fact == 1) fact <- n / i
        m <- floor(fact * (n - i + 1))
        y0<-env.rand(m) 
        u <- runif(m)
        j <- 1
      }
      if(accept) break
    }
    y[i] <- z
    
  }
  y
}

#==================================================================================================================
#=============== PROFILING ===============================================================================================

# 1st Version
library("profvis")
profvis({
  vMsim.input <- function(n,func, envelop,env.rand) {
    y <- numeric(n)
    count <- 0
    
    for(i in 1:n) {
      ratio <- 0 
      u <- 1
      
      while(u > ratio) {
        y0 <- env.rand(1) 
        ratio <-func(y0)/envelop(y0)
        u <- runif(1)
        
      }
      y[i] <- y0
      
    }
    y
  }
vMsim.input(100000,func,envelop,gen.envelop)
})

#----------------------------------------------------------------------------------------------------------------
# also have to adjust the two functions to take the same parameters

# 2nd Version-with vectorized form
profvis({
  vMsim.vectorized <- function(n,func,envelop,env.rand) {
    m<-n
    fact <- 1
    y <- numeric(n)
    u <- runif(m)
    y0 <- env.rand(m)
    y0<-y0*(y0>=0)
    j<-1
    for(i in 1:n) {
      repeat{
        z <- y0[j]
        accept <- u[j] <= func(z)/envelop(z)
        j <- j + 1
        if(j > m) {
          if(fact == 1) fact <- n / i
          m <- floor(fact * (n - i + 1))
          y0<- env.rand(m)
          u <- runif(m)
          j <- 1
        }
        if(accept) break
      }
      y[i] <- z
      
    }
    y
  }
  vMsim.vectorized(100000,func,envelop,gen.envelop)
})

# the vectorized form runs about 4 times faster speed now compared to the other one

#================================================================================================================

# ++++ we can add an RCPP implementation

#plots- they look nice now
simde<-vMsim.vectorized(1000000,func,envelop,gen.envelop) # with 1 million iterations pretty good estimate of the density

simde2<-vMsim.input(100000,func,envelop,gen.envelop,trace=TRUE)

vMsim.input(100000,func,envelop,gen.envelop,trace=TRUE)

plot(density(simde))
plot(density(simde2))
#===============================================================================================================================

# Plotting a few outcomes

# target density
targetdens<-function(x) { exp(-x^3+x)}
rangee<-seq(from=0,to=2,by=0.001)
plot(rangee,targetdens(rangee),type="l") 
lines(rangee,targetdens(rangee),type="l")

#add our envelope to the plot
envelope<-dnorm(x=rangee,mean = 1/sqrt(3)-0.03,sd =0.54)
plot(rangee,envelope*2.5,col="red")  

#check its histogram
simdens<-vMsim(100000)
hist(simdens[simdens>0],breaks=100,probability  =TRUE ) #it does not look bad 

#Compare our results witht the target density
plot(rangee,targetdens(rangee)/1.55,type="l") 
lines(density(simdens),col="Red")

# some plots
plot(density(simdens[simdens>0]))
hist(simdens[simdens>0], breaks = seq(0,2, length.out = 20), prob = TRUE)
lines(rangee,targetdens(rangee)/1.55,col="red")

plot(rangee,targetdens(rangee),type="l")

#=======================================================================================================
install.packages("Rfast")
library("Rfast") # for colMins

# the logarithm of our target density
logf <- function(y) log(f(y))
plot(logf,from=0,to=5) # checking if it will be log concave: yes

# slope, derivative of logf 
dlogf <- function(y) { -3*y^2+1 }
plot(dlogf,0,5)

##======================= Adaptive rejection sampling ##================================================================================

prepare_affine <- function(logfunc,dlogfunc,tangents,xmin=0,xmax=Inf){
  f <- function(x) exp(logfunc(x))
  logf <- function(x) logfunc(x)
  dlogf <- function(x) dlogfunc(x)
  pieces <- length(tangents) 
  slope <- dlogfunc(tangents) 
  inter <- logfunc(tangents) - slope*tangents
  aff <- function(x) {
    v <- inter+slope %*% t(as.matrix(x))
    exp( colMins(v,value=TRUE) )
  }
  
  cross <-(inter[2:pieces]-inter[1:(pieces-1)]) / (slope[1:(pieces-1)]-slope[2:pieces])# z(i)
  aucs <- exp(inter)/slope*( exp(slope*c(cross,xmax))-exp(slope*c(xmin,cross)) ) 
  cumauc <- aucs
  for(i in 2:pieces) cumauc[i] <- cumauc[i-1]+aucs[i]
  auc <- cumauc[pieces] 
  
  invcumaff <- function(y){
    m1 <- matrix(rep(y,pieces-1),ncol=pieces-1) 
    m2 <- matrix(rep(cumauc[1:(pieces-1)]/auc,each=length(y)),ncol=pieces-1) 
    fin <- 1+rowsums(m1>m2)  
    log( (y*auc-c(0,cumauc)[fin])*slope[fin]*exp(-inter[fin])+exp(slope[fin]*c(xmin,cross)[fin]) )/slope[fin] 
  }
  
  raff <- function(n){  
    invcumaff(runif(n))
  }
  
  list(aff=aff, cross=cross, cumauc=cumauc, auc=auc, invcumaff=invcumaff, raff=raff)
}


a <- prepare_affine(logf,dlogf,c(0.1,0.6,1,2))

plot(f,0,3)
targetdens<-function(x) { exp(-x^3+x)}
rangee<-seq(from=0,to=2,by=0.001)
# showing the envelope covering the target density
plot(rangee,targetdens(rangee),type="l") 
lines(a$aff(d),x=d,col="Red") 

# Now i have to use these affine function in my original rejections samlping algorithm

a <- prepare_affine(logf,dlogf,c(0.1,0.6,1,2))

a <- prepare_affine(logf,dlogf,optimal$par)

vMa <- vMsim(100000, f, a$aff, a$raff)
vMa<-vMsim.input(100000,func,a$aff, a$raff)
vMa$c # rejection ratio, around 6-7%

#plot the results
plot(rangee,targetdens(rangee)/1.55,type="l") 
lines(density(vMa$y),type="l",col="Red")

# find the best tangent points, by minimizing the area under the curve



auc_of_affine(c(.5,1,2))
auc_of_affine(c(.4,1,2))
auc_of_affine(c(.4,1,1.5))

auc_of_affine <- function(pts) prepare_affine(logf,dlogf,pts)$auc
optimal <- optim(c(.5,1,1.5,2),auc_of_affine)
optimal$par 


# the optimal parameters

#Plot the optimal envelope
a <- prepare_affine(logf,dlogf,optimal$par)
 
plot(a$aff(rangee),x=rangee,col="Red",lwd=3,type="l")
lines(rangee,targetdens(rangee),type="l",lwd=1) 

#compare different number of tangent points, and their performance
env<-function(x){2.5*dnorm(x,mean = 1/sqrt(3)-0.03,sd =0.54)}
env.gen<-function(x){2.5*qnorm(runif(x),mean = 1/sqrt(3)-0.03,sd =0.54)}

env.gen()


a <- prepare_affine(logf,dlogf,optimal$par)
a$


result<-vMsim.input(10000,func,a$aff,a$raff)

#------------------------------------------------------------------------------------------------------------------
# create a version with the envelope as input variable

vMsim.input <- function(n, envelop,enveloperandom) {
  y <- numeric(n)
  count <- 0
  
  for(i in 1:n) {
    ratio <- 0 
    u <- 1
    
    while(u > ratio) {
      y0 <- enveloperandom(1) 
      ratio <-exp(-y0^3+y0)*(y>=0)/envelop(y0)
      u <- runif(1)
      count <- count + 1
    }
    y[i] <- y0
    
  }
  #adding tracing of the rejection ratio
  y
  c<-(count-n)/count
  
  list(y=y,c=c)
}

#=======================================================================================================================

# Benchmarking
install.packages("microbenchmark")
library("microbenchmark")

#compare the performance of the different versions
tmp<-microbenchmark(
vMsim(10000),
vMsim.input(10000,func,envelop,gen.envelop,trace=TRUE),
vMsim.vectorized(10000,func,envelop,gen.envelop),
vMsim_cpp(10000)
)

summary(tmp) # cpp version is way faster than the others

# We should also compare their results in terms of the density function.!!

#====================================================================================================================
install.packages("Rcpp")
library(Rcpp)
# RCPP IMPLEMENTATION
#this loads the function from the rcpp file
sourceCpp("C:/Users/Peter/Desktop/University/KU/Computational statistics/Exam/R-files/Assignment 2/univar.cpp") #it works


trial(1)
trial(2)
trial(3)
trial(4)
trial(5)

plot(density(vMsim(10000)))
plot(density(vMsim_cpp(10000)))


