install.packages("profvis")
library("profvis")

##### ASSIGNMENT 1 ##########

#===================####### 1st VERSION #############=============================================7
#Epanechnikov kernel

kernDensEpan <- function (x, h, m = 512) {
  rg <- range(x)  
  xx <- seq(rg[1] -  3*h, rg[2] +3*h, length.out = m) 
  y <- numeric(m)   
  for (i in seq_along(xx))
    for (j in seq_along(x))
      
      if ((xx[i] - x[j])/h<=1 & (xx[i] - x[j])/h>=-1){   
        y[i] <- y[i] + 0.75-0.75*((xx[i] - x[j])^2 / (h^2))
      }   
  
  y <- y /( h*length(x))
  list(x = xx, y = y)
} 


#define the kernel 
K <- function(x) (x>=-1)*(x<=1)*(1-x*x)*0.75
#=====================================================================================================================

######## 2nd VERSION ######### - a vectorized implementation

kernDensEpan2 <- function (x, h, m = 512) {
  rg <- range(x)  
  xx <- seq(rg[1] - 3 * h, rg[2] + 3 * h, length.out = m) 
  y <- numeric(m)  
  
  for (i in seq_along(xx)){
    
    y[i] <- sum(0.75-0.75*((xx[i] - x[x<=h+xx[i] & x>=xx[i]-h])^2 / (h^2))) 
    
  } 
  y <- y /(h * length(x))
  list(x = xx, y = y)  
}

#===================================================================================================================

#=========####### 3RD VERSION #####---- getting rid of the for cycles/loops completely =======================

kernDensEpan3 <- function (x, h, m = 512) {
  rg <- range(x)
  xx <- seq(rg[1] - 3 * h, rg[2] + 3 * h, length.out = m) 
  y <- numeric(m)   
  # outer: doing all possible product combinations of two sequence
  y <- outer(xx, x,function(xx, x) (0.75-0.75*((xx - x)^2 / (h^2)))) 
  y<-y*(y>0)
  y <- rowMeans(y)/h
  list(x = xx, y = y)  
}

#==============================================================================================================

#----------------------- COMPARISON -----------------------------------------------------------------

# Comparison of the different methods- all of them give the same results

range(kernDensEpan(logF12,0.25*sqrt(5))$y-kernDensEpan2(logF12,0.25*sqrt(5))$y) 

range(kernDensEpan(logF12,0.25*sqrt(5))$y-kernDensEpan3(logF12,0.25*sqrt(5))$y)
# we see that the difference is indeed very small due to only some rounding errors

#check whether the functions give the same results
all.equal(kernDensEpan(logF12,0.25)$y,kernDensEpan2(logF12,0.25)$y)

all.equal(kernDensEpan(logF12,0.25)$y,kernDensEpan3(logF12,0.25)$y)

#--- BENCHMARKING

#simulate dataset for Benchmarking

simul.data<-rnorm(10000,5,5)
plot(density(simul.data))

lines(kernDensEpan(simul.data,0.25*sqrt(5)),col="Red")

simul.data<-rnorm(10000,5,5)

#include later the Rcpp version here as well.
install.packages("microbenchmark")
library("microbenchmark")

Kernel.First.version<-kernDensEpan(simul.data,0.25*sqrt(5))
Kernel.Second.version<-kernDensEpan2(simul.data,0.25*sqrt(5))
Kernel.Third.version<-kernDensEpan3(simul.data,0.25*sqrt(5)) 

# The byte version does not give much improvement here, since all the base R functions are already compiled this way
Byte.Kernel.Second.version <- compiler::cmpfun(kernDensEpan2)

#faster than the first one but slower than the second

microbenchmark(Kernel.First.version,Kernel.Second.version,Kernel.Third.version)

microbenchmark(kernDensEpan(simul.data,0.25),kernDensEpan2(simul.data,0.25),kernDensEpan3(simul.data,0.25 ),Byte.Kernel.Second.version(simul.data,0.25))

#second one is the fastest by far

#### ---------TESTING----------------------------------------------------------------------------------------------

#maybe some examples for different h parameter choices, undersmooth and oversmooth

infrared<- read.table("C:/Users/Peti/Desktop/KU/Computational statistics/CSwR-master/CSwR-master/data/infrared.txt",header = TRUE)
F12 <- infrared$F12
logF12 <- log(F12)
### test all these with some real data now

hist(logF12,type = "l",prob=TRUE,breaks=50) #undersmooth
lines(kernDensEpan(logF12,0.25*sqrt(5)),type = "l",col="Red",lwd=2.5) 

hist(logF12,type = "l",prob=TRUE,breaks=50,main = "h = 1") 
lines(kernDensEpan(logF12,1),type = "l",col="Red",lwd=2.5) 

hist(logF12,type = "l",prob=TRUE,breaks=50,main = "h = 0.5") 
lines(kernDensEpan(logF12,0.6),type = "l",col="Red",lwd=2.5) 

hist(logF12,type = "l",prob=TRUE,breaks=50) 
lines(kernDensEpan(logF12,0.1),type = "l",col="Red",lwd=2.5) 


plot(density(logF12, bw=0.25)) 
lines(kernDensEpan(logF12,0.25*sqrt(5)),type = "l") 
lines(kernDensEpan2(logF12,0.25*sqrt(5)),type = "l",col="Red") 
lines(kernDensEpan3(logF12,0.25*sqrt(5)),type = "l",col="Green") # all cover each other

#Compare the differences between kerndensity and the built-in density function
plot(kernDensEpan(logF12,0.25*sqrt(5))$x,density(logF12,bw = 0.25)$y-kernDensEpan(logF12,0.25*sqrt(5))$y)
# !!!! have to adjust the gridpoint somehow

# The two implementations only compute approximately the same cause density relies on certain approximation for run time efficiency

# getting the exact same as the density

plot(density(logF12,0.249,kernel = "epanechnikov"),lwd=2.5)
lines(kernDensEpan(logF12,0.249*sqrt(5)),type = "l",col="red",lwd=1.5) 

plot(kernDensEpan(logF12,0.249)$y-density(logF12,0.249)$y,type="l")

cbind(kernDensEpan(logF12,0.249*sqrt(5))$y,density(logF12,0.249,kernel = "epanechnikov")$y)
#Finding the ideal value for h

#Silverman rule of Thumb
bw.nrd(x = logF12) #0.2439433 it is a good initial guess for the bandwidth


F12<-as.data.frame(F12) #create format from data

#==============================================================================================================

#### Cross-validation for the bandwith selection #####

ngroups <- 10 
group <- cut(seq(1,length(logF12)),breaks=ngroups,labels=FALSE)
group<-sample(group)
max(group)

kernDensGroups <- function(dat, group, h=1, m = 512) { 
  
  ugroups <- unique(group)
  ngroups <- max(group)
  x <- numeric(ngroups*m)
  y <- numeric(ngroups*m)
  for(i in 1:ngroups){  
    testData <- dat[group == ugroups[i]]
    trainData <- dat[group != ugroups[i]]
    xy <- kernd.crossvalid(testData,trainData, h, m) 
    x[((i-1)*m+1):(i*m)] <- xy$x
    y[((i-1)*m+1):(i*m)] <- xy$y
  }
  list(x=x, y=y)
}

CVlogL <- function(h){
  kdg <- kernDensGroups(logF12, group, h=h)
  -sum(log(pmax(kdg$y,1e-20))) 
}  

hCVmin <- optimize(CVlogL,c(0,1))  #0.2684
hCVmin<-hCVmin$minimum
hCVmin

#------------------------------------------------------------------------------------------------------------------ 

#Plotting and comparison with the optimal bandwidth

plot(density(logF12, hCVmin/sqrt(5), kernel="epanechnikov"))
lines(kernDens(logF12,h=hCVmin),type="l",col="Red")  # why am I not getting this result

# Now it works
hist(logF12,probability = TRUE,breaks = 50)
lines(density(logF12, hCVmin/sqrt(5),kernel="epanechnikov"))
lines(kernDensEpan(logF12,hCVmin),type = "l",col="Red",lwd=2)

plot(kernDensEpan(logF12,hCV),type = "l",col="Red")
lines(kernDens(logF12,hCV),type = "l",col="Red")

######################################################################################################################xx

# Plot differences between the built-in density and my implementation

f_hat_dens <- density(logF12, hCVmin/sqrt(5), kernel="epanechnikov")
f_hat <- kernDens(logF12, hCVmin, minx=min(f_hat_dens$x), maxx=max(f_hat_dens$x))
plot(f_hat, type = "l", lwd = 4, xlab = "x", ylab = "Density")
lines(f_hat_dens, col = "red", lwd = 2)
plot(f_hat$x, f_hat$y - f_hat_dens$y, 
     type = "l", lwd = 2, xlab = "x", ylab = "Difference")


# density relies on certain approximations for run time efficiency.

#------------------------------------------------------------------------------------------------------------------------------
 # RCCP version
install.packages("Rcpp")
library(Rcpp)
# RCPP IMPLEMENTATION
#this loads the function from the rcpp file
sourceCpp("C:/Users/Peter/Desktop/University/KU/Computational statistics/Exam/R-files/Assignment 2/univar.cpp") #it works




## End of Presentation ##

