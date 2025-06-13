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

#Dimensions of the panel
N <- 200 #Number of units
T <- 20 #Number of periods

#Generation of exogenous
x1 <- rnorm(n=T*N, m=1,sd=1) 
x2 <- rnorm(n=T*N, m=1,sd=1) 
x <- cbind(x1,x2)

#Generation of Fixed effects
FE <- x[1:N,]+ cbind(exp(rnorm(n=N,m=0,sd=0.5)),exp(rnorm(n=N,m=0,sd=0.5))) #Generation of fixed effects
#FE <- exp(rnorm(n=N,m=0,sd=0.5)) #Generation of fixed effects
FE <- FE - mean(FE) #ATTENTION: normalization of FE in order that the sum is equal to 0, i.e. they are distributed around 0. This is also compatible with a common intercept to be estimated in some way.

print(paste("Correlation between exogenous and FE = ", round(cor(x[1:N,1],FE[,1]),digits=3)))
print(paste("Correlation between exogenous and FE = ", round(cor(x[1:N,2],FE[,2]),digits=3)))

#Time effetcs
TE <- seq(1:T)/T #Generation of time effects
TE <-  TE - mean(TE)

#Generation of random components
epsilon <- rnorm(n=T*N, m=0,sd=0.1) 

#Parameters of nonlinear relationship between exogenous and endogenous
beta.1.1 <- 1
beta.2.1 <- 0.2
beta.3.1 <- -0.02
beta.4.1 <- 0.001

beta.1.2 <- 0.5
beta.2.2 <- -0.5
beta.3.2 <- 0.02
beta.4.2 <- -0.001

FE.matrix <-  t(matrix(FE,ncol=N,nrow=T,byrow=T))
TE.matrix <-  t(matrix(TE,ncol=N,nrow=T))
x1.matrix <- t(matrix(x1,ncol=N,nrow=T))
x2.matrix <- t(matrix(x2,ncol=N,nrow=T))
epsilon.matrix <- t(matrix(epsilon,ncol=N,nrow=T))
  
y.matrix <- FE.matrix + TE.matrix + beta.1*x + beta.2*x^2 +beta.3*x^3 + beta.4*x^4  + epsilon.matrix

matplot(t(y.matrix),type="l")


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
    diff.x.CS <- abs(x.matrix[i,j]-x.matrix[,j]) 
    
    #Set the distance with itself to the maximum TS
    diff.x.CS[i] <- max(diff.x.CS) 
    
    #Distance between xs over time
    diff.x.TS <- abs(x.matrix[i,j]-x.matrix[i,]) 

    #Set the distance with itself to the maximum TS
    diff.x.TS[j] <- max(diff.x.TS) 
    
    #Take the x with the minimum distance cross-section
    min.index.CS <- which.min(diff.x.CS) 
    
    #Take the x with the minimum distance time series
    min.index.TS <- which.min(diff.x.TS) 
    
    #Record the distance between xs cross section
    diff.x.CS.array[i,j] <- diff.x.CS[min.index.CS] 
    
    #Record the distance between xs time series
    diff.x.TS.array[i,j] <- diff.x.TS[min.index.TS] 
    
    #Record the distance between ys cross-section
    #diff.y.CS.array[i,j] <- (diff.x.TS.array[i,j]/diff.x.CS.array[i,j])*(y.matrix[i,j]- y.matrix[min.index.CS,j])
    diff.y.CS.array[i,j] <- (y.matrix[i,j]- y.matrix[min.index.CS,j])
    
    #Record the distance between ys time series
    diff.y.TS.array[i,j] <- y.matrix[i,j]- y.matrix[i,min.index.TS] 
    
    #Set the array for FE estimation
    #diff.FE.array[i,i,j] <- - (diff.x.TS.array[i,j]/diff.x.CS.array[i,j])
    diff.FE.array[i,i,j] <- - 1
    
    #diff.FE.array[i,min.index.CS,j] <- (diff.x.TS.array[i,j]/diff.x.CS.array[i,j])
    diff.FE.array[i,min.index.CS,j] <- 1
    
    #Set the array for TE estimation
    diff.TE.array[i,j,j] <- 1
    
    diff.TE.array[i,min.index.TS,j] <- - 1
    
    #Set the array for the estimation of FE and TE jointly
    diff.total.array[i, ,j] <- c(diff.FE.array[i,,j], diff.TE.array[i,,j])
    
  }
  
}

#Stacking of all variables
delta.y.stacked <-  c(diff.y.TS.array - diff.y.CS.array)
delta.FE.TE.stacked <- do.call(rbind, lapply(1:T,function(t) diff.total.array[ , ,t] ))

#Numerical optimization
score <- function(FE.TE.vec,C,delta.Y){
  
  sum( (delta.Y - C%*%FE.TE.vec)^2 )
  
}

## Initial values of FE estimation ####
FE.TE.ini <- matrix(0,nrow=(N+T),ncol=1)

t0 <- Sys.time()
res <- optim(par = FE.TE.ini, fn = score,C=delta.FE.TE.stacked, delta.Y=delta.y.stacked, method =  "L-BFGS-B")
print(Sys.time() - t0)

#The outcome of estimation
est.FE.TE <- res$par

#Using the generalized inverse matrix
FE.TE.star <- (ginv(t(delta.FE.TE.stacked)%*%delta.FE.TE.stacked))%*%t(delta.FE.TE.stacked)%*%matrix(delta.y.stacked,nrow=N*T,ncol=1)

#Using a small nudge
lambda <- 1e-6
to.be.inverted <- t(delta.FE.TE.stacked)%*%delta.FE.TE.stacked
to.be.inverted_reg <- to.be.inverted + diag(lambda, nrow(to.be.inverted))
inv_delta.FE.TE.stacked_reg <- solve(to.be.inverted_reg)
FE.TE.star_reg <- inv_delta.FE.TE.stacked_reg%*%t(delta.FE.TE.stacked)%*%matrix(delta.y.stacked,nrow=N*T,ncol=1)


# Diagnostics of estimates

#Sum of FE
print(rbind(
  c(sum(FE),sum(TE)),
  c(sum(est.FE.TE[1:N]),sum(est.FE.TE[(N+1):(T+N)])),
  c(sum(FE.TE.star[1:N]),sum(FE.TE.star[(N+1):(T+N)])),
  c(sum(FE.TE.star_reg[1:N]),sum(FE.TE.star_reg[(N+1):(T+N)]))
)
)

#MSE
print( c(mean((FE-est.FE.TE[1:N])^2), mean((FE-FE.TE.star[1:N])^2),mean((FE-FE.TE.star_reg[1:N])^2)) ) 
print( c(mean((TE-est.FE.TE[(N+1):(T+N)])^2), mean((TE-FE.TE.star[(N+1):(T+N)])^2),mean((TE-FE.TE.star_reg[(N+1):(T+N)])^2)) )


#Explained variance of FE and TE estimation
print(summary(lm(FE  ~ I(est.FE.TE[1:N]) )))
print(summary(lm(FE  ~ I(FE.TE.star[1:N]) )))
print(summary(lm(FE  ~ I(FE.TE.star_reg[1:N]) )))
print(summary(lm(TE  ~ I(est.FE.TE[(N+1):(T+N)]) )))
print(summary(lm(TE  ~ I(FE.TE.star[(N+1):(T+N)]) )))
print(summary(lm(TE  ~ I(FE.TE.star_reg[(N+1):(T+N)]) )))

#Plot of the estimated FE
plot(est.FE.TE[1:N],FE,pch=19,ylim=range(c(FE,est.FE.TE[1:N],FE.TE.star[1:N],FE.TE.star_reg[1:N])),xlim=range(c(TE,est.FE.TE[1:N],FE.TE.star[1:N],FE.TE.star_reg[1:N])),cex=0.75)
points(FE.TE.star[1:N],FE,col="green",pch=19)
points(FE.TE.star_reg[1:N],FE,col="red",pch=19)
lines(c(-10,10),c(-10,10))
grid()

plot(est.FE.TE[(N+1):(T+N)],TE,pch=19,ylim=range(c(TE,est.FE.TE[(N+1):(T+N)],FE.TE.star[(N+1):(T+N)],FE.TE.star_reg[(N+1):(T+N)])),xlim=range(c(TE,est.FE.TE[(N+1):(T+N)],FE.TE.star[(N+1):(T+N)],FE.TE.star_reg[(N+1):(T+N)])),cex=0.75)
points(FE.TE.star[(N+1):(T+N)],TE,col="green",pch=19)
points(FE.TE.star_reg[(N+1):(T+N)],TE,col="red",pch=19)
lines(c(-10,10),c(-10,10))
grid()
