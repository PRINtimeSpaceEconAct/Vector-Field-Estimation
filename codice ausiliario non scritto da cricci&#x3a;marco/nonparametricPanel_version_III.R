#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# The estimate of a nonparametric panel with time and fixed effects
#
# Update: May 3, 2025
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#############
# TODO
# 1) the case of dimension 2

rm(list = ls())

library(sm)
library(MASS)
library(Matrix)

#set.seed(124)

#Setting of the analysis
using.distance.obs <- TRUE
numericalOptimizationProcedure <- FALSE

# parameters #####
N = 100
T = 10

# data Generation ####
#set.seed(1)

# genero i FE ####
FE = mvrnorm(N,mu=c(0,0),Sigma=0.01*diag(2))
FE[,1] = FE[,1] - sum(FE[,1])/N
FE[,2] = FE[,2] - sum(FE[,2])/N

# genero i TE ####
TE = mvrnorm(T,mu=c(0,0),Sigma=0.01*diag(2))
TE[,1] = TE[,1] - sum(TE[,1])/T
TE[,2] = TE[,2] - sum(TE[,2])/T

X0 = mvrnorm(N, mu=c(0,0),Sigma = 0.25*diag(2)) + FE + array(rep(TE[1,],each=N),dim=c(N,2))
X = array(NA,dim = c(N,2,T+1))
X[,,1] = X0

# example 1 - double well ----
VF <- function(X){
  # X = (x,y)
  # U(X) = x^4 - x^2 + y^2
  # VF(X) = -grad U(X) = -(4x^3 - 2x, 2y)
  return( -0.01*c(4*X[1]^3 - 2*X[1], 2*X[2]) )
}

# example 2 -- single well ----
# VF <- function(X){
#     # X = (x,y)
#     # U(X) = x^2 + y^2
#     # VF(X) = -grad U(X) = -(2x, 2y)
#     return( -0.1*c(2*X[1], 2*X[2]) )
# }

# example 3 -- rotation ----
# M = matrix(c(cos(pi/2), -sin(pi/2), sin(pi/2), cos(pi/2)),nrow=2,ncol=2)
# VF <- function(X){
#     # X = (x,y), theta = pi/4
#     return (0.1*(M %*% X - X))
# }


for (t in 1:T){
  X[,,t+1] = X[,,t] + t(apply(X[,,t], 1, VF)) + FE + array(rep(TE[t,],each=N),dim=c(N,2)) + 
    + mvrnorm(N, mu=c(0,0),Sigma = 0.001*diag(2))
}

# calcolo i Delta
Delta = (X[,,2:(T+1)] - X[,,1:T])

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#One dimensional case #####
y.matrix <- Delta[,1,]
dim(y.matrix)

x.matrix <- X[,1,1:T]
dim(x.matrix)

FE <- FE[,1]
TE <- TE[,1]


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Estimate of FE via C-estimator ####
diff.y.array <- array(NA,dim=c(N,T))
diff.x.array  <- array(NA,dim=c(N,T))
diff.FE.array <- array(0,dim=c(N,N,T))

for (j in (1:T)){
  for (i in (1:N)){
    
    #Distance between xs
    diff.x <- abs(x.matrix[i,j]-x.matrix[,j]) 
    
    #Set the distance with itself to the maximum
    diff.x[i] <- max(diff.x) 
    
    #Take the x with the minimum distance
    min.index <- which.min(diff.x) 
    
    #Record the distance between xs
    diff.x.array[i,j] <- diff.x[min.index] 
    
    #Record the distance between ys 
    diff.y.array[i,j] <- y.matrix[i,j]- y.matrix[min.index,j] 
  
    #Set the array for FE estimation
    diff.FE.array[i,i,j] <- 1
    diff.FE.array[i,min.index,j] <- -1
    
    }
  
}

#Stacking of all variables
delta.y.stacked <-  c(diff.y.array)
delta.x.stacked <-  c(diff.x.array)
delta.FE.stacked <- do.call(rbind, lapply(1:T,function(t) diff.FE.array[ , ,t] ))

#Numerical optimization
score <- function(FE.vec,A,delta.Y){
  
  sum( (delta.Y - A%*%FE.vec)^2 )
 
}

## Initial values of FE estimation ####
FE.ini <- matrix(0,nrow=N,ncol=1)
#FE.ini <- FE

if (numericalOptimizationProcedure==TRUE){
t0 <- Sys.time()
res <- optim(par = FE.ini, fn = score,A=delta.FE.stacked, delta.Y=delta.y.stacked, method =  "L-BFGS-B")
print(Sys.time() - t0)
#The outcome of estimation
est.FE <- res$par
est.FE.one.eq <- lm(FE  ~ est.FE)
}else{
  est.FE <- NA
  est.FE.one.eq <- NA
}

#Using the generalized inverse matrix
FE.star <- (ginv(t(delta.FE.stacked)%*%delta.FE.stacked))%*%t(delta.FE.stacked)%*%matrix(delta.y.stacked,nrow=N*T,ncol=1)

#Using the sparse matrix routines
#delta.FE.stacked_sparse <- Matrix(delta.FE.stacked, sparse = TRUE)
#inv_delta.FE.stacked <- solve(t(delta.FE.stacked_sparse)%*%delta.FE.stacked_sparse)
#FE.star.sparse <- inv_delta.FE.stacked%*%t(delta.FE.stacked_sparse)%*%matrix(delta.y.stacked,nrow=N*T,ncol=1)

# #Using a small nudge
# lambda <- 1e-6
# to.be.inverted <- t(delta.FE.stacked)%*%delta.FE.stacked
# to.be.inverted_reg <- to.be.inverted + diag(lambda, nrow(to.be.inverted))
# inv_delta.FE.stacked_reg <- solve(to.be.inverted_reg)
# FE.star_reg <- inv_delta.FE.stacked_reg%*%t(delta.FE.stacked)%*%matrix(delta.y.stacked,nrow=N*T,ncol=1)
FE.star_reg <- NA

## Diagnostics of estimates ####
#Sum of FE
print(c(sum(FE),sum(est.FE),sum(FE.star),sum(FE.star_reg))) #Sum of FE

#MSE
print( c(mean((FE-est.FE)^2), mean((FE-FE.star)^2),mean((FE-FE.star_reg)^2)) ) # MSE

#Explained variance of FE estimation
est.FE.star.one.eq <- lm(FE  ~ FE.star)
#est.FE.star_reg.one.eq <- lm(FE  ~ FE.star_reg )

#print(summary(lm(FE  ~ est.FE)))
#print(summary(lm(FE  ~ FE.star)))
#print(summary(lm(FE  ~ FE.star_reg )))

## Plot of the estimated FE ####
plot(FE.star,FE,pch=19,ylim=range(c(FE,est.FE,FE.star,FE.star_reg),na.rm = T),xlim=range(c(FE,est.FE,FE.star,FE.star_reg),na.rm = T),cex=0.75)
#points(FE,FE,col="green",pch=19)
#points(FE.star_reg,FE,col="red",pch=19)
lines(c(-10,10),c(-10,10))
grid()


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Estimate of TE via C-estimator ####
diff.y.array <- array(NA,dim=c(N,T))
diff.x.array  <- array(NA,dim=c(N,T))
diff.TE.array <- array(0,dim=c(N,T,T))

for (j in (1:T)){
  for (i in (1:N)){
    
    #Distance between xs
    diff.x <- abs(x.matrix[i,j]-x.matrix[i,]) 
    
    #Set the distance with itself to the maximum
    diff.x[j] <- max(diff.x) 
    
    #Take the x with the minimum distance
    min.index <- which.min(diff.x) 
    
    #Record the distance between xs
    diff.x.array[i,j] <- diff.x[min.index] 
    
    #Record the distance between ys 
    diff.y.array[i,j] <- y.matrix[i,j]- y.matrix[i,min.index] 
    
    #Set the array for TE estimation
    diff.TE.array[i,j,j] <- 1
    diff.TE.array[i,min.index,j] <- -1
    
  }
  
}

#Stacking of all variables
delta.y.stacked <-  c(diff.y.array)
delta.x.stacked <-  c(diff.x.array)
delta.TE.stacked <- do.call(rbind, lapply(1:T,function(t) diff.TE.array[ , ,t] ))

#Numerical optimization
score <- function(TE.vec,B,delta.Y){
  
  sum( (delta.Y - B%*%TE.vec)^2 )
  
}

## Initial values of FE estimation ####
TE.ini <- matrix(0,nrow=T,ncol=1)

if (numericalOptimizationProcedure==TRUE){
t0 <- Sys.time()
res <- optim(par = TE.ini, fn = score,B=delta.TE.stacked, delta.Y=delta.y.stacked, method =  "L-BFGS-B")
print(Sys.time() - t0)
#The outcome of estimation
est.TE <- res$par
est.TE.one.eq <- lm(TE  ~ est.TE)
}else{
  est.TE <- NA
  est.TE.one.eq <- NA
  
}

#Using the generalized inverse matrix
TE.star <- (ginv(t(delta.TE.stacked)%*%delta.TE.stacked))%*%t(delta.TE.stacked)%*%matrix(delta.y.stacked,nrow=N*T,ncol=1)

# #Using a small nudge
# lambda <- 1e-6
# to.be.inverted <- t(delta.TE.stacked)%*%delta.TE.stacked
# to.be.inverted_reg <- to.be.inverted + diag(lambda, nrow(to.be.inverted))
# inv_delta.TE.stacked_reg <- solve(to.be.inverted_reg)
# TE.star_reg <- inv_delta.TE.stacked_reg%*%t(delta.TE.stacked)%*%matrix(delta.y.stacked,nrow=N*T,ncol=1)
TE.star_reg <- NA

## Diagnostics of estimates ####

#Sum of FE
print(c(sum(TE),sum(est.TE),sum(TE.star),sum(TE.star_reg))) 

#MSE
print( c(mean((TE-est.TE)^2), mean((TE-TE.star)^2),mean((TE-TE.star_reg)^2)) ) 

#Explained variance of FE estimation
est.TE.star.one.eq <- lm(TE  ~ TE.star)
#est.TE.star_reg.one.eq <- lm(TE  ~ TE.star_reg )
#print(summary(lm(TE  ~ est.TE)))
#print(summary(lm(TE  ~ TE.star)))
#print(summary(lm(TE  ~ TE.star_reg )))

## Plot of the estimated TE ####
plot(TE.star,TE,pch=19,ylim=range(c(TE,est.TE,TE.star,TE.star_reg),na.rm = T),xlim=range(c(TE,est.TE,TE.star,TE.star_reg),na.rm = T),cex=0.75)
#points(TE,TE,col="green",pch=19)
#points(TE.star_reg,TE,col="red",pch=19)
lines(c(-10,10),c(-10,10))
grid()


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Joint estimate of FE and TE via C-estimator ####

diff.y.CS.array <- array(NA,dim=c(N,T))
diff.y.TS.array <- array(NA,dim=c(N,T))

diff.x.CS.array  <- array(NA,dim=c(N,T))
diff.x.TS.array  <- array(NA,dim=c(N,T))

diff.FE.array <- array(0,dim=c(N,N,T))
diff.TE.array <- array(0,dim=c(N,T,T))

diff.total.array <- array(0,dim=c(N,(N+T),T))


for (j in (1:T)){
  for (i in (1:N)){
    
    #Distance between xs cross section
    diff.x.CS <- x.matrix[,j]-x.matrix[i,j] 
    
    #Set the distance with itself to the maximum TS
    diff.x.CS[i] <- max(abs(diff.x.CS)) 
    
    #Distance between xs over time
    diff.x.TS <- x.matrix[i,]-x.matrix[i,j] 
    
    #Set the distance with itself to the maximum TS
    diff.x.TS[j] <- max(abs(diff.x.TS)) 
    
    #Take the x with the minimum distance cross-section
    min.index.CS <- which.min(abs(diff.x.CS)) 
    
    #Take the x with the minimum distance time series
    min.index.TS <- which.min(abs(diff.x.TS)) 
    
    #Record the distance between xs cross section
    diff.x.CS.array[i,j] <- diff.x.CS[min.index.CS] 
    
    #Record the distance between xs time series
    diff.x.TS.array[i,j] <- diff.x.TS[min.index.TS] 
    
    if(using.distance.obs==FALSE){
    
    #Record the distance between ys cross-section
    #diff.y.CS.array[i,j] <- (diff.x.TS.array[i,j])*(y.matrix[i,j]- y.matrix[min.index.CS,j])
    diff.y.CS.array[i,j] <- (y.matrix[i,j]- y.matrix[min.index.CS,j])
    
    #Record the distance between ys time series
    #diff.y.TS.array[i,j] <- (diff.x.CS.array[i,j])*(y.matrix[i,j]- y.matrix[i,min.index.TS]) 
    diff.y.TS.array[i,j] <- (y.matrix[i,j]- y.matrix[i,min.index.TS]) 
    
    #Set the array for FE estimation
    #diff.FE.array[i,i,j] <- (diff.x.TS.array[i,j])
    diff.FE.array[i,i,j] <-  1
    #diff.FE.array[i,min.index.CS,j] <- -(diff.x.TS.array[i,j])
    diff.FE.array[i,min.index.CS,j] <- -1
    
    #Set the array for TE estimation
    #diff.TE.array[i,j,j] <- diff.x.CS.array[i,j]
    diff.TE.array[i,j,j] <- 1
    #diff.TE.array[i,min.index.TS,j] <- - diff.x.CS.array[i,j]
    diff.TE.array[i,min.index.TS,j] <- - 1
    }else{
      #Record the distance between ys cross-section
      diff.y.CS.array[i,j] <- (diff.x.TS.array[i,j])*(y.matrix[i,j]- y.matrix[min.index.CS,j])
      #diff.y.CS.array[i,j] <- (y.matrix[i,j]- y.matrix[min.index.CS,j])
      
      #Record the distance between ys time series
      diff.y.TS.array[i,j] <- (diff.x.CS.array[i,j])*(y.matrix[i,j]- y.matrix[i,min.index.TS]) 
      #diff.y.TS.array[i,j] <- (y.matrix[i,j]- y.matrix[i,min.index.TS]) 
      
      #Set the array for FE estimation
      diff.FE.array[i,i,j] <- (diff.x.TS.array[i,j])
      #diff.FE.array[i,i,j] <-  1
      diff.FE.array[i,min.index.CS,j] <- -(diff.x.TS.array[i,j])
      #diff.FE.array[i,min.index.CS,j] <- -1
      
      #Set the array for TE estimation
      diff.TE.array[i,j,j] <- diff.x.CS.array[i,j]
      #diff.TE.array[i,j,j] <- 1
      diff.TE.array[i,min.index.TS,j] <- - diff.x.CS.array[i,j]
      #diff.TE.array[i,min.index.TS,j] <- - 1
    }
    #Set the array for the estimation of FE and TE jointly
    diff.total.array[i, ,j] <- c(diff.FE.array[i,,j], - diff.TE.array[i,,j])
    
  }
  
}

#Stacking of all variables
delta.y.stacked <-  c(diff.y.CS.array -diff.y.TS.array )
delta.FE.TE.stacked <- do.call(rbind, lapply(1:T,function(t) diff.total.array[ , ,t] ))

#Numerical optimization
score <- function(FE.TE.vec,C,delta.Y){
  
  sum( (delta.Y - C%*%FE.TE.vec)^2 )
  
}

## Initial values of FE estimation ####
FE.TE.ini <- matrix(0,nrow=(N+T),ncol=1)

if (numericalOptimizationProcedure==TRUE){
  t0 <- Sys.time()
  res <- optim(par = FE.TE.ini, fn = score,C=delta.FE.TE.stacked, delta.Y=delta.y.stacked, method =  "L-BFGS-B")
  print(Sys.time() - t0)
  #The outcome of estimation
  est.FE.TE <- res$par
  est.FE <- lm(FE  ~ I(est.FE.TE[1:N]) )
  est.TE <- lm(TE  ~ I(est.FE.TE[(N+1):(T+N)]) )
  
}else{
  
  est.FE.TE <- NA
  est.FE <- NA
  est.TE <- NA
}

#Using the generalized inverse matrix
FE.TE.star <- (ginv(t(delta.FE.TE.stacked)%*%delta.FE.TE.stacked))%*%t(delta.FE.TE.stacked)%*%matrix(delta.y.stacked,nrow=N*T,ncol=1)

# #Using a small nudge
# lambda <- 1e-5
# to.be.inverted <- t(delta.FE.TE.stacked)%*%delta.FE.TE.stacked
# to.be.inverted_reg <- to.be.inverted + diag(lambda, nrow(to.be.inverted))
# inv_delta.FE.TE.stacked_reg <- solve(to.be.inverted_reg)
# FE.TE.star_reg <- inv_delta.FE.TE.stacked_reg%*%t(delta.FE.TE.stacked)%*%matrix(delta.y.stacked,nrow=N*T,ncol=1)
FE.TE.star_reg <- NA


# Diagnostics of estimates

#Sum of FE
print(rbind(
  #c(sum(FE),sum(TE)),
  #c(sum(est.FE.TE[1:N]),sum(est.FE.TE[(N+1):(T+N)])),
  c(sum(FE.TE.star[1:N]),sum(FE.TE.star[(N+1):(T+N)])),
  c(sum(FE.TE.star_reg[1:N]),sum(FE.TE.star_reg[(N+1):(T+N)]))
)
)

#MSE
#print( c(mean((FE-est.FE.TE[1:N])^2), mean((FE-FE.TE.star[1:N])^2),mean((FE-FE.TE.star_reg[1:N])^2)) ) 
#print( c(mean((TE-est.FE.TE[(N+1):(T+N)])^2), mean((TE-FE.TE.star[(N+1):(T+N)])^2),mean((TE-FE.TE.star_reg[(N+1):(T+N)])^2)) )


#Explained variance of FE and TE estimation
est.FE.star <- lm(FE  ~ I(FE.TE.star[1:N]) )
#est.FE.star_reg <- lm(FE  ~ I(FE.TE.star_reg[1:N]) )
est.TE.star <- lm(TE  ~ I(FE.TE.star[(N+1):(T+N)]) )
#est.TE.star_reg <- lm(TE  ~ I(FE.TE.star_reg[(N+1):(T+N)]) )
#print(summary(lm(FE  ~ I(est.FE.TE[1:N]) )))
#print(summary(lm(FE  ~ I(FE.TE.star[1:N]) )))
#print(summary(lm(FE  ~ I(FE.TE.star_reg[1:N]) )))
#print(summary(lm(TE  ~ I(est.FE.TE[(N+1):(T+N)]) )))
#print(summary(lm(TE  ~ I(FE.TE.star[(N+1):(T+N)]) )))
#print(summary(lm(TE  ~ I(FE.TE.star_reg[(N+1):(T+N)]) )))

#Plot of the estimated FE
plot(FE.TE.star[1:N],FE,pch=19,ylim=range(c(FE,est.FE.TE[1:N],FE.TE.star[1:N],FE.TE.star_reg[1:N]),na.rm = T),xlim=range(c(TE,est.FE.TE[1:N],FE.TE.star[1:N],FE.TE.star_reg[1:N]),na.rm = T),cex=0.75)
#points(FE.TE[1:N],FE,col="green",pch=19)
#points(FE.TE.star_reg[1:N],FE,col="red",pch=19)
lines(c(-10,10),c(-10,10))
grid()

plot(FE.TE.star[(N+1):(T+N)],TE,pch=19,ylim=range(c(TE,est.FE.TE[(N+1):(T+N)],FE.TE.star[(N+1):(T+N)],FE.TE.star_reg[(N+1):(T+N)]),na.rm=T),xlim=range(c(TE,est.FE.TE[(N+1):(T+N)],FE.TE.star[(N+1):(T+N)],FE.TE.star_reg[(N+1):(T+N)]),na.rm=T),cex=0.75)
#points(FE.TE[(N+1):(T+N)],TE,col="green",pch=19)
#points(FE.TE.star_reg[(N+1):(T+N)],TE,col="red",pch=19)
lines(c(-10,10),c(-10,10))
grid()

# print("Numerical minimization")
# print(cbind(coef(est.FE.one.eq),coef(est.TE.one.eq),
#   coef(est.FE),coef(est.TE)
#   )
# )

print("Solving the system of equations")
print( 
cbind(coef(est.FE.star.one.eq),coef(est.TE.star.one.eq),
      coef(est.FE.star),coef(est.TE.star)
)
)

print(c(median(diff.x.CS.array),median(diff.x.TS.array)))
print(median(c(diff.x.CS.array/diff.x.TS.array)))
print(c(sd(diff.x.CS.array),sd(diff.x.TS.array)))

library(sm)
summary(c(diff.x.TS.array/diff.x.CS.array))
sm.density(c(diff.x.TS.array/diff.x.CS.array),model="Normal")
summary(c(diff.x.TS.array/diff.x.CS.array))

sm.density(c(diff.x.CS.array),model="Normal",ylim=c(0,1.2))
sm.density(c(diff.x.TS.array),add=T,lwd=1,col="red")

print(cor(c(diff.x.CS.array),c(diff.x.TS.array)))
print(c(var(c(diff.x.CS.array)),var(c(diff.x.TS.array))))
