

library(microbenchmark)
library(ggplot2)
library(compiler)

infrared <- read.table("C:/Users/Jens/Documents/GitHub/CSwR/data/infrared.txt", header = TRUE)
logF12 <- log(infrared$F12)
summary(logF12)
hist(logF12)
hist(logF12, breaks=seq(from=-3,to=8,by=1))
hist(logF12, breaks=seq(from=-3,to=8,by=.5) ,freq=FALSE)
hist(logF12, breaks=seq(from=-3,to=8,by=.1) ,freq=FALSE)
hist(logF12, breaks=seq(from=min(logF12),to=max(logF12),length.out=200),freq=FALSE)



#Bandwidth selection with cross-validation

# Cross validated kernel density

# Epanechnikov kernel
K <- function(x) (x>=-1)*(x<=1)*(1-x*x)*.75

kernDens <- function (x, h=1, m=512, minx=NA, maxx=NA) {
  rg <- range(x)
  if(is.na(minx)) minx <- rg[1]-3*h
  if(is.na(maxx)) maxx <- rg[2]+3*h
  xx <- seq(minx, maxx, length.out = m)
  ## xx is equivalent to grid points in 'density'
  
  y <- numeric(m) ## The evaluations, initialized as a vector of zeroes
  ## The actual computation is done using nested for-loops. The outer loop
  ## is over the grid points, and the inner loop is over the data points.
  for (i in seq_along(x))
    for (j in seq_along(xx))
      y[j] <- y[j] + K( (xx[j]-x[i])/h )
  y <- y / (length(x)*h)
  list(x = xx, y = y)
}
kd <- kernDens(logF12, h=.3)
plot(kd)
length(kd$y)


ngroups <- 10
group <- floor(runif(length(logF12))*ngroups)+1
gr<-1
h<-.3
m<-512
kd <- kernDens(logF12[group==3], h=.3)
plot(kd)


kernDensGroups <- function(dat, group, h=1, m = 512) {
  ugroups <- unique(group)
  ngroups <- length(ugroups)
  x <- numeric(ngroups*m)
  y <- numeric(ngroups*m)
  print(h)
  for(gr in 1:ngroups){
    trainData <- dat[group == ugroups[gr]]
    xy <- kernDens(trainData, h, m)
    x[((gr-1)*m+1):(gr*m)] <- xy$x
    y[((gr-1)*m+1):(gr*m)] <- xy$y
  }
  list(x=x, y=y)
}
kdg <- kernDensGroups(logF12, group, h=.3)
plot(kdg)
summary(kdg$y)
log()


CVlogL <- function(h){
  kdg <- kernDensGroups(logF12, group, h=h)
  -sum(log(pmax(kdg$y,1e-6))) # 1e-3: .21, 1e-6: .25, 1e-10: .27 (larger h when larger penalty)
}
hCVmin <- optimize(CVlogL,c(0,1))
hCV <- hCVmin$minimum
# plot


plot(density(logF12, hCV/sqrt(5), kernel="epanechnikov"))
lines(kernDens(logF12,h=hCV),type="l",col="Red")

?density

f_hat_dens <- density(logF12, 0.25/sqrt(5), kernel="epanechnikov")
f_hat <- kernDens(logF12, 0.25, minx=min(f_hat_dens$x), maxx=max(f_hat_dens$x))
plot(f_hat, type = "l", lwd = 4, xlab = "x", ylab = "Density")
lines(f_hat_dens, col = "red", lwd = 2)
plot(f_hat$x, f_hat$y - f_hat_dens$y, 
     type = "l", lwd = 2, xlab = "x", ylab = "Difference")
# density relies on certain approximations for run time efficiency.


summary(f_hat$x)
summary(f_hat_dens$x)
sum(abs(f_hat$x-f_hat_dens$x))




kernd.crossvalid <- function (testData,trainData, h, m = 512) {
  #returns the max and min of the vector
  ## xx is equivalent to grid points in 'density'
  x<-trainData
  rg <- range(testData) 
  xx <- seq(rg[1] -  3*h, rg[2] +3*h, length.out = m) #giving the starting and ending points of the sequence, not quite sure about the 3h, that is the standard of the x-axis of the density
  y <- numeric(m)   ## The evaluations, initialized as a vector of zeroes
  ## The actual computation is done using nested for-loops. The outer loop
  ## is over the grid points, and the inner loop is over the data points.
  for (j in seq_along(xx)){
    for (i in seq_along(x)){#the length of the dataset, in each loop we add one element to the y vector
      y[j] <- y[j] + K( (xx[j]-x[i])/h )
      # if ( (x[i] - xx[j])/h<=1 & (x[i] - xx[j])/h>=-1) {   #the index condition
      #   y[j] <- y[j] + 0.75-0.75*((xx[j] - x[i])^2 / (h^2))
      # }   
    }
  }
  y<-(y)/(length(x)*h)
  list(x = xx, y = y)
} 
a <- kernd.crossvalid(logF12,logF12, h=.3)
plot(a)

#splitting up the data into equal sized samples

# logF12 <- as.data.frame(logF12) #it needs to be a data.frame for the rest of the code to run
logF12[sample(length(logF12))]

#Randomly shuffle the data
yourdata <- logF12[sample(length(logF12))]

#Create 10 equally sized folds
folds <- cut(seq(1,nrow(yourdata)),breaks=10,labels=FALSE)

#Perform 10 fold cross validation
for(i in 1:10){
  #Segement your data by fold using the which() function 
  testIndexes <- which(folds==i,arr.ind=TRUE)
  testData <- yourdata[testIndexes]
  trainData <- yourdata[-testIndexes ]
  
}


yourdata<-as.data.frame(logF12[sample(nrow(logF12)),])
#Create 10 equally size folds
folds <- cut(seq(1,nrow(yourdata)),breaks=10,labels=FALSE)
#Perform 10 fold cross validation



# trying to find the optimum of the log likelihood
# but not quite sure if I have defined the loglikelihood properly
optimcrossval<-function(h){ 
  for(i in 1:10){
    #Segement your data by fold using the which() function 
    testIndexes <- which(folds==i,arr.ind=TRUE)
    testData <- yourdata[testIndexes, ]
    trainData <- yourdata[-testIndexes, ]
    
    m<-512
    x<-trainData
    rg <- range(testData) 
    xx <- seq(rg[1] -  3*h, rg[2] +3*h, length.out = m) 
    y <- numeric(m)   
    
    for (j in seq_along(xx)){
      for (i in seq_along(x)){
        if ( (x[i] - xx[j])/h<=1 & (x[i] - xx[j])/h>=-1) {   #the index condition
          y[j] <- y[j] + 0.75-0.75*((xx[j] - x[i])^2 / (h^2))
        }   
      }
    }
    y<-(y)/(length(testData)*h)
    
    summa<-summa+sum(log(y[y>0]))
    
    #Use the test and train data partitions however you desire...
  } 
  summa
  
  
  optimise(optimcrossval,interval = c(0.00001,1),maximum = TRUE) # the h value i get does not make sense
  
}