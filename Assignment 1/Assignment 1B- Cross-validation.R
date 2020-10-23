#Bandwidth selection with cross-validation

# Cross validated kernel density

infrared<- read.table("infrared.txt",header=TRUE)
F12 <- infrared$F12
lofF12<-log(F12)
logF12<-as.data.frame(logF12)

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
      if ( (x[i] - xx[j])/h<=1 & (x[i] - xx[j])/h>=-1) {   #the index condition
        y[j] <- y[j] + 0.75-0.75*((xx[j] - x[i])^2 / (h^2))
      }   
    }
  }
  y<-(y)/(length(x)*h)
  list(x = xx, y = y)
} 


#splitting up the data into equal sized samples

logF12<-as.data.frame(logF12) #it needs to be a data.frame for the rest of the code to run
logF12[sample(nrow(logF12))]
#Randomly shuffle the data
yourdata<-logF12[sample(nrow(logF12))]
#Create 10 equally size folds
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