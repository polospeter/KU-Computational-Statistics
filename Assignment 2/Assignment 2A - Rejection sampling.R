# Assignments 2.1: Rejection Sampling

install.packages("Rfast")
library("Rfast") # for colMins
library("profvis")


f <- function(x) { (x>=0)*exp(-x^3+x) } # this is how the function is defined


const <- 1/integrate(f,0,Inf)$value
env <- function(x){ (x>=0)*2*exp(-x^2/2) }

env(-1)

d <- seq(0,3,.01)

plot(f,0,3)
lines(y=env(d),x=d,col="blue")


# env>f because env starts higher and goes to 0 slower than f since x^3 > x^2
ggplot(data.frame(x=d), aes(x = x)) +
  stat_function(fun = eq, geom="line") +
  xlab("x") + ylab("y") +
  stat_function(fun = env, geom="line", colour='red' )


# simplest rejection sampling, many improvements possible
vMsim <- function(n, func, env, randenv) {
  y <- numeric(n)
  for(i in 1:n) {
    ratio <- 0
    u <- 1
    while(u > ratio) {
      y0 <- randenv(1)
      ratio <- func(y0)/env(y0)
      u <- runif(1)
    }
    y[i] <- y0
  }
  y
}

rnormpos <- function(n) qnorm(runif(n)/2+.5)
vMg <- vMsim(100000, f, env, rnormpos)

hist(vMg,breaks=100,prob=TRUE)
plot(density(vMg))
curve(f(x)*const, 0,5, col = "blue", lwd = 2, add = TRUE) # looks proportional :)


# adaptive rejection sampling using a piecewise log-affine envelope

# week 3b lecture, slide 27-

# is f log-concave?
logf <- function(y) log(f(y))
plot(logf,from=0,to=5) # yes

# need the slope of f
dlogf <- function(y) { -3*y^2+1 }
plot(dlogf,0,5)


### set up log-affine envelope

prepare_affine <- function(logfunc,dlogfunc,tangents,xmin=0,xmax=Inf){
  # the tangents inout is basically the x(i) parameter
  f <- function(x) exp(logfunc(x))
  logf <- function(x) logfunc(x)
  dlogf <- function(x) dlogfunc(x)
  pieces <- length(tangents) 
  slope <- dlogfunc(tangents) # slopes for affine tangents, a(i)
  inter <- logfunc(tangents) - slope*tangents # getting b(i), the intercepts with the y-axis
  aff <- function(x) {
    v <- inter+slope %*% t(as.matrix(x))
    exp( colMins(v,value=TRUE) )
  }
  
  # this is the formula of z(i)
  cross <-(inter[2:pieces]-inter[1:(pieces-1)]) / (slope[1:(pieces-1)]-slope[2:pieces])# z(i)
  aucs <- exp(inter)/slope*( exp(slope*c(cross,xmax))-exp(slope*c(xmin,cross)) ) #area under the curve
  cumauc <- aucs
  for(i in 2:pieces) cumauc[i] <- cumauc[i-1]+aucs[i]
  auc <- cumauc[pieces] # getting the total area under the curve, Fi(x)
  
  #getting the inverse cumulative function
  invcumaff <- function(y){
    #find which interval it belongs to
    m1 <- matrix(rep(y,pieces-1),ncol=pieces-1) 
    m2 <- matrix(rep(cumauc[1:(pieces-1)]/auc,each=length(y)),ncol=pieces-1)
    foundin <- 1+rowsums(m1>m2)
    log( (y*auc-c(0,cumauc)[foundin])*slope[foundin]*exp(-inter[foundin]) + exp(slope[foundin]*c(xmin,cross)[foundin]) )/slope[foundin]
  }
  
  raff <- function(n){
    invcumaff(runif(n))
  }
  
  list(f=f, logf=logf, dlogf=dlogf, tangents=tangents, pieces=pieces, slope=slope, inter=inter, aff=aff, cross=cross, cumauc=cumauc, auc=auc, invcumaff=invcumaff, raff=raff)
}

a <- prepare_affine(logf,dlogf,c(0.1,0.6,1,2))

plot(f,0,3)
lines(a$aff(d),x=d,col="Red")
  
a$invcumaff(c(.1,.3,.5,.7,.9,.99,.999))

hist(a$raff(10000),breaks=100,prob=TRUE,xlim=c(0,3))
lines(y=a$aff(d)/a$auc,x=d,col="red")


vMa <- vMsim(100000, f, a$aff, a$raff)

hist(vMa,breaks=100,prob=TRUE,xlim=c(0,3),type="l")
lines(y=const*f(d),x=d,col="red")

# find the best tangent points
auc_of_affine <- function(pts) prepare_affine(logf,dlogf,pts)$auc
auc_of_affine(c(.5,1,2))
auc_of_affine(c(.4,1,2))
auc_of_affine(c(.4,1,1.5))

optimal <- optim(c(.5,1,2),auc_of_affine)
optimal$par

b <- prepare_affine(logf,dlogf,optimal$par)

c(a$auc,b$auc)

plot(f,0,3)
lines(y=a$aff(d) ,x=d, col="blue")
lines(y=b$aff(d) ,x=d, col="red")

vMb <- vMsim(100000, f, b$aff, b$raff) # simulate using optimal log-affine envelope

a$auc/b$auc

hist(vMb,breaks=100,prob=TRUE,xlim=c(0,3))
lines(y=const*f(d),x=d,col="red")



# to do:
# plot-function
# speed comparison:
# - number of affines
# - better vMsim
# - a$raff in Rcpp


