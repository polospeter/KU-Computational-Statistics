

# install.packages("magic")
library(profvis)
library(ggplot2)
library(magic)
library(microbenchmark)

# https://cswr.nrhstat.org/mixed-models.html - example 7.5
# https://cswr.nrhstat.org/em.html
# https://cswr.nrhstat.org/newton-type-algorithms.html - sparsity
# https://cswr.nrhstat.org/fisher-information.html - Fisher Information


#------------------------------- GENERATE DATA -----------------------------------------------------------------------------------

#it is a more generalized solution,
# all because the epsilon ij
simData <- function(n, b0, nu, s2, seed=1234567){
  if(seed != 0) set.seed(seed) #set seed
  m <- length(n) #how many elements does the matrix contains
  nobs <- sum(n) # sum of all the elements of the vector
  lo <- rep(1,m)
  hi <- rep(nobs,m)
  
  for(i in 1:(m-1)){         # the reason why we do this cause the ni is row specific
    lo[i+1] <- lo[i]+n[i]    #low and high point of that interval
    hi[i] <- lo[i+1]-1
  }
  eps <- rnorm(nobs,0,s2) #generate epsilons
  z <- numeric(nobs)
  
  zorg <- rnorm(m) #generate standard normal variables
  
  for(i in 1:m) z[lo[i]:hi[i]] <- zorg[i]
  
  y <- b0 + nu*z + eps # create the Y variable
  
  sumy <- numeric(m) # create a vector 
  for(i in 1:m) sumy[i] <- sum(y[lo[i]:hi[i]])
  
  list(n=n, m=m, nobs=nobs, lo=lo, hi=hi, b0=b0, nu=nu, s2=s2, seed=seed, z=z, eps=eps, y=y, sumy=sumy)
}
a <- simData(rep(c(3,4,5,6,7),100),1,1,1) # so we need to add a vector for n (length of ni in each row)

#-----------------------------------------------------------------------------------------------------------------------------------------------

cbind( a$z[1:10], a$eps[1:10], a$y[1:10] )


# full information ML---- 
x <- cbind(1,a$z) 
bet <- solve(t(x) %*% x) %*% t(x) %*% matrix(a$y,nrow=length(a$y))
b0 <- bet[1]
nu <- bet[2]
yhat <- x %*% bet
s2 <- var(a$y-yhat)
c(bet,s2)


#----------------- Create the M-step and the E-step function ---------------------------------------------------------------------

#profiling

profvis({
# E-step: Obtain conditional expectations of unobserved xi and zeta
Estep <- function(par, obj){
  b0 <- par[1]    # input: we give here our initial parameters 
  nu <- par[2]
  s2 <- par[3]
  mat11 <- diag(obj$m)  #creates an identity matrix
  mat12 <- matrix(0, nrow=obj$m, ncol=obj$nobs)
  for(i in 1:obj$m) mat12[i,obj$lo[i]:obj$hi[i]] <- nu
  mat21 <- t(mat12)  
  mat22 <- matrix(0, nrow=obj$nobs, ncol=obj$nobs)
  for(i in 1:obj$m) mat22[obj$lo[i]:obj$hi[i],obj$lo[i]:obj$hi[i]] <- (diag(drop(s2),obj$n[i]) + nu*nu)
  mat22inv <- solve(mat22)
  xi <- mat12 %*% mat22inv %*% (obj$y-b0) # conditional means
  zeta <- diag(mat11 - mat12 %*% mat22inv %*% mat21) + xi*xi # cond. 2nd moments
  list(xi=drop(xi), zeta=drop(zeta))
}
 e <- Estep(c(1,1,1),a) # the inputs here for beta0, nu and sigma are 1, these are only initial guesses the process will converge to their true values
})


# M-step: Obtain best parameter estimates given conditional expectations
Mstep <- function(xi_zeta, obj) {
  xi <- xi_zeta$xi
  zeta <- xi_zeta$zeta
  mat <- matrix(c(obj$nobs, rep(sum(obj$n*xi),2), sum(obj$n*zeta)) ,2,2)
  betanu <- solve(mat) %*% matrix(c(sum(obj$y),sum(obj$sumy*xi)),nrow=2)
  b0 <- betanu[1]
  nu <- betanu[2]
  s2 <- ( sum(obj$y*obj$y) - sum(obj$n*(b0*b0+nu*nu*zeta+2*b0*nu*xi)) ) / obj$nobs
  c(b0, nu, s2)
}
# par <- Mstep(e, a)


#-----------------------Combine the two steps: generel EM algorithm ------------------------------------------------------------------------

EMest <- function(par, obj, Efct=Estep, Mfct=Mstep, eps=1e-6) {
  repeat{
    par0 <- par
    par <- Mfct(Efct(par, obj), obj)
    if(sum((par - par0)^2) <= eps*(sum(par^2) + eps)) # NB: parameterization
      break
  } 
  par
}


est <- EMest(c(1,1,1),a)



EMrec <- function(par, obj, Efct=Estep, Mfct=Mstep, eps=1e-6) {
  par0 <- par
  par <- Mfct(Efct(par, obj), obj)
  if(sum((par - par0)^2) <= eps * (sum(par^2) + eps))
    return(par)
  EMrec(par, obj, Efct=Efct, Mfct=Mfct, eps=eps)
}
rec <- EMrec(c(1,1,1),a)

# both are equally slow!

# microbenchmark(EMest(c(1,1,1),a), EMrec(c(1,1,1),a), times=1 )
# Unit: seconds
# expr      min       lq     mean   median       uq      max neval
# EMest(c(1, 1, 1), a) 46.10671 46.10671 46.10671 46.10671 46.10671 46.10671     1
# EMrec(c(1, 1, 1), a) 47.41743 47.41743 47.41743 47.41743 47.41743 47.41743     1


profvis({
# Faster E-step
Estep3 <- function(par, obj){
  b0 <- par[1]
  nu <- par[2]
  s2 <- par[3]
  xi3 <- matrix(rep(0,obj$m), nrow=obj$m, ncol=1)
  zeta3 <- matrix(0, nrow=obj$m, ncol=1)
  for(i in 1:obj$m) {
    col_sums <- colSums( solve(diag(drop(s2),obj$n[i]) + nu*nu) )
    xi3[i] <- sum(col_sums * (obj$y[obj$lo[i]:obj$hi[i]]-b0)) * nu
    zeta3[i] <- 1 - sum(col_sums) *nu*nu + xi3[i]*xi3[i]
  }
  list(xi=drop(xi3), zeta=drop(zeta3))
}
e <- Estep3(c(1,1,1),a)
})


# estimation using faster E-step
est3 <- EMest(c(1,1,1),a, Efct=Estep3)
rec3 <- EMrec(c(1,1,1),a, Efct=Estep3)

# recursive and non-recursive are practically the same

# microbenchmark(EMest(c(1,1,1),a,Efct=Estep3), EMrec(c(1,1,1),a,Efct=Estep3), times=100 )
# Unit: milliseconds             expr      min       lq     mean   median       uq      max neval cld
# EMest(c(1, 1, 1), a, Efct = Estep3) 133.5165 139.1904 144.3635 141.9765 143.7739 302.7266   100   a
# EMrec(c(1, 1, 1), a, Efct = Estep3) 133.9082 138.8030 141.0167 140.6150 142.8848 152.0991   100   a








?microbenchmark


mat11 <- diag(a$m)

for(ii in 1:1000){
  print(ii)
  
# E-step: Set up matrixes to obtain conditional expectations of unobserved
  
mat12 <- matrix(0, nrow=a$m, ncol=a$nobs)
  for(i in 1:a$m) mat12[i,a$lo[i]:a$hi[i]] <- nu
  mat21 <- t(mat12)

mat22 <- matrix(0, nrow=a$nobs, ncol=a$nobs)
for(i in 1:a$m) mat22[a$lo[i]:a$hi[i],a$lo[i]:a$hi[i]] <- (diag(drop(s2),a$n[i]) + nu*nu)
mat22inv <- matrix(0, nrow=nobs, ncol=nobs)
for(i in 1:m) mat22inv[lo[i]:hi[i],lo[i]:hi[i]] <- solve(diag(drop(s2),n[i]) + nu*nu)
# sum(abs(solve(mat22)-mat22inv))

mat12_mat22inv <- matrix(0, nrow=m, ncol=nobs)
# for(i in 1:m) mat12_mat22inv[i,lo[i]:hi[i]] <- mat12[i,lo[i]:hi[i]] %*% solve(diag(drop(s2),n[i]) + nu*nu)
for(i in 1:m) mat12_mat22inv[i,lo[i]:hi[i]] <- matrix(nu,nrow=1,ncol=n[i]) %*% solve(diag(drop(s2),n[i]) + nu*nu)

# sum(abs(mat12 %*% mat22inv - mat12_mat22inv))

xi <- (mat12 %*% solve(mat22)) %*% (a$y-b0) # conditional means

xi2 <- matrix(rep(0,m), nrow=m, ncol=1)
y_b <- matrix(a$y-b0, nrow=nobs, ncol=1)
for(i in 1:m) xi2[i] <- matrix(nu,nrow=1,ncol=n[i]) %*% (solve(diag(drop(s2),n[i]) + nu*nu) %*% y_b[lo[i]:hi[i]])
# sum(abs(xi-xi2))

xi3 <- matrix(rep(0,m), nrow=m, ncol=1)
for(i in 1:m) xi3[i] <- colSums(solve(diag(drop(s2),n[i]) + nu*nu)) %*% y_b[lo[i]:hi[i]] * nu
# sum(abs(xi-xi3))

xi4 <- matrix(rep(0,a$m), nrow=a$m, ncol=1)
for(i in 1:a$m) xi4[i] <- sum( colSums(solve(diag(drop(s2),a$n[i]) + nu*nu)) * (a$y[a$lo[i]:a$hi[i]]-b0) )* nu
# sum(abs(xi-xi4))

zeta <- diag(mat11 - mat12 %*% mat22inv %*% mat21) + xi*xi # cond. 2nd moments

zeta2 <- matrix(0, nrow=m, ncol=1)
for(i in 1:m) zeta2[i] <- 1 - sum( solve(diag(drop(s2),n[i])+nu*nu) ) *nu*nu + xi[i]*xi[i]
# sum(abs(zeta-zeta2))

print(c(mean(xi),mean(zeta)))

# M-step: Obtain best parameter estimates given conditional expectations

mat <- matrix(c(nobs, rep(sum(n*xi),2), sum(n*zeta)) ,2,2)
betanu <- solve(mat) %*% matrix(c(sum(a$y),sum(a$sumy*xi)),nrow=2)
b0 <- betanu[1]
nu <- betanu[2]
s2 <- ( sum(a$y*a$y) - sum(n*(b0*b0+nu*nu*zeta+2*b0*nu*xi)) ) / nobs
print(c(b0,nu,s2))
}

mean(zeta)













#diagonal with block-diagonal elements added
covmat <- function(n,nu,s2){
  c <- matrix(0, nrow=sum(n), ncol=sum(n))
  counter <- 1
  for(i in 1:length(n)){
    counter1 <- counter+n[i]-1
    c[counter:counter1,counter:counter1] <- (diag(s2,n[i]) + nu)
    counter <- counter+n[i]
  }
  c
}
covmat(c(1,2,3),.3,.1)
cm <- covmat(a$n,.3,.1)
cm

bet0 <- mean(a$y)
xi0 <- numeric(a$m)
zeta0 <- numeric(a$m)
for(i in 1:a$m){
  xi0[i] <- mean(a$y[a$lo[i]:a$hi[i]]-bet0)
  zeta0[i] <- max(.01,var(a$y[a$lo[i]:a$hi[i]])-s2)
}
bet0
xi0
zeta0

m <- matrix(c(a$m, rep(sum(a$n*xi0),2), sum(a$n*zeta0)),2,2)
m
betanu <- solve(m) %*% matrix(c(sum(a$y),sum(a$y*zeta0)),nrow=2)

s2 <- ( sum(a$y^2) - sum(a$n*(betanu[1]+betanu[2]*zeta0+2*betanu[1]*betanu[2]*xi0)) ) / a$m



mat11 <- diag(m)
