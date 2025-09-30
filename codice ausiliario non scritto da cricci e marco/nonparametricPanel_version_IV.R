#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# The estimate of a nonparametric panel with time and fixed effects
#
# Update: June 4, 2025
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#############
# TODO
# 1) the case of dimension 2

rm(list = ls())

library(sm)
library(MASS)
library(Matrix)
library(KernSmooth)

#set.seed(124)

#Setting of the analysis
using.distance.obs <- TRUE
numericalOptimizationProcedure <- FALSE

# parameters #####
N = 1000
T = 20

#Nonparametric function
m_x <- function(x) {
  1/(1+exp(-x))
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Data generation from Lee, Mukkherjee and Ullah (2019)####

## errors ####
u <- matrix(rnorm(N*T,0,1),nrow=N,ncol=T)

## X #####
X <- matrix(rnorm(N*T,2,4^2),nrow=N,ncol=T)


## FE ####
v <- runif(N,min=-0.5,max=0.5)
FE <- 0.5*rowMeans(X)+v
FE <- FE - mean(FE) #our normalization

## TE ####
TE <- rnorm(T,0,0.3)
TE <- TE - mean(TE) #our normalization

## y 
y <- m_x(X) + matrix(FE,nrow=N,ncol=T) + matrix(TE,nrow=N,ncol=T,byrow=T) + u

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#One dimensional case #####
y.matrix <- y
dim(y.matrix)

x.matrix <- X
dim(x.matrix)

mean(apply(y,FUN=sd,MARGIN=1))
mean(apply(y,FUN=sd,MARGIN=2))

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
print("@@@@@@@@@@@@@@@@@@@@@@@@@@@",quote=FALSE)
print("",quote=FALSE)
print("Estimation of FE alone via C-estimator")

print("Sum")
print(c(sum(FE),sum(est.FE),sum(FE.star),sum(FE.star_reg)),quote=FALSE) #Sum of FE

#MSE
print("MSE")
print( c(mean((FE-est.FE)^2), mean((FE-FE.star)^2),mean((FE-FE.star_reg)^2)) ) # MSE

#Explained variance of FE estimation
est.FE.star.one.eq <- lm(FE  ~ FE.star)
#est.FE.star_reg.one.eq <- lm(FE  ~ FE.star_reg )

#print(summary(lm(FE  ~ est.FE)))
#print(summary(lm(FE  ~ FE.star)))
#print(summary(lm(FE  ~ FE.star_reg )))

## Plot of the estimated FE ####
plot(FE.star,FE,pch=19,ylim=range(c(FE,est.FE,FE.star,FE.star_reg),na.rm = T),xlim=range(c(FE,est.FE,FE.star,FE.star_reg),na.rm = T),cex=0.75,main="FE alone")
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

#Sum of TE
print("@@@@@@@@@@@@@@@@@@@@@@@@@@@",quote=FALSE)
print("",quote=FALSE)
print("Estimation of TE alone via C-estimator")
print("Sum")
print(c(sum(TE),sum(est.TE),sum(TE.star),sum(TE.star_reg))) 

#MSE
print("MSE")
print( c(mean((TE-est.TE)^2), mean((TE-TE.star)^2),mean((TE-TE.star_reg)^2)) ) 

#Explained variance of FE estimation
est.TE.star.one.eq <- lm(TE  ~ TE.star)
#est.TE.star_reg.one.eq <- lm(TE  ~ TE.star_reg )
#print(summary(lm(TE  ~ est.TE)))
#print(summary(lm(TE  ~ TE.star)))
#print(summary(lm(TE  ~ TE.star_reg )))

## Plot of the estimated TE ####
plot(TE.star,TE,pch=19,ylim=range(c(TE,est.TE,TE.star,TE.star_reg),na.rm = T),xlim=range(c(TE,est.TE,TE.star,TE.star_reg),na.rm = T),cex=0.75,main="TE alone")
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
print("@@@@@@@@@@@@@@@@@@@@@@@@@@@",quote=FALSE)
print("",quote=FALSE)
print("Joint estimation of FE and TE via C-estimator")

print("Sum FE")
print(rbind(
  #c(sum(FE),sum(TE)),
  #c(sum(est.FE.TE[1:N]),sum(est.FE.TE[(N+1):(T+N)])),
  c(sum(FE.TE.star[1:N]),sum(FE.TE.star[(N+1):(T+N)])),
  c(sum(FE.TE.star_reg[1:N]),sum(FE.TE.star_reg[(N+1):(T+N)]))
)
)

#MSE
print("MSE of FE")
print( c(mean((FE-est.FE.TE[1:N])^2), mean((FE-FE.TE.star[1:N])^2),mean((FE-FE.TE.star_reg[1:N])^2)) ) 
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
plot(FE.TE.star[1:N],FE,pch=19,ylim=range(c(FE,est.FE.TE[1:N],FE.TE.star[1:N],FE.TE.star_reg[1:N]),na.rm = T),xlim=range(c(TE,est.FE.TE[1:N],FE.TE.star[1:N],FE.TE.star_reg[1:N]),na.rm = T),cex=0.75,main="FE (joint estimated with TE)")
#points(FE.TE[1:N],FE,col="green",pch=19)
#points(FE.TE.star_reg[1:N],FE,col="red",pch=19)
lines(c(-10,10),c(-10,10))
grid()

plot(FE.TE.star[(N+1):(T+N)],TE,pch=19,ylim=range(c(TE,est.FE.TE[(N+1):(T+N)],FE.TE.star[(N+1):(T+N)],FE.TE.star_reg[(N+1):(T+N)]),na.rm=T),xlim=range(c(TE,est.FE.TE[(N+1):(T+N)],FE.TE.star[(N+1):(T+N)],FE.TE.star_reg[(N+1):(T+N)]),na.rm=T),cex=0.75,main="TE (joint estimated with FE)")
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

#print(c(median(diff.x.CS.array),median(diff.x.TS.array)))
#print(median(c(diff.x.CS.array/diff.x.TS.array)))
#print(c(sd(diff.x.CS.array),sd(diff.x.TS.array)))

#library(sm)
#summary(c(diff.x.TS.array/diff.x.CS.array))
#sm.density(c(diff.x.TS.array/diff.x.CS.array),model="Normal")
#summary(c(diff.x.TS.array/diff.x.CS.array))

#sm.density(c(diff.x.CS.array),model="Normal",ylim=c(0,1.2))
#sm.density(c(diff.x.TS.array),add=T,lwd=1,col="red")

#print(cor(c(diff.x.CS.array),c(diff.x.TS.array)))
#print(c(var(c(diff.x.CS.array)),var(c(diff.x.TS.array))))


y <- m_x(X) + matrix(FE,nrow=N,ncol=T) + matrix(TE,nrow=N,ncol=T,byrow=T) + u


y_filtered <- y - matrix(FE.TE.star[1:N],nrow=N,ncol=T) - matrix(FE.TE.star[(N+1):(T+N)],nrow=N,ncol=T,byrow=T)
plot(c(X),c(y_filtered),pch=19,cex=0.25)
fit <- locpoly(c(y_filtered), c(X), bandwidth = 0.5)
lines(fit)
curve(m_x,from=-40,to=40,add=TRUE,col="red")
