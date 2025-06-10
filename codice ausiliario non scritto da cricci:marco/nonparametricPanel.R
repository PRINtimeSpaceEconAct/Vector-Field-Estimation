#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# The estimate of a nonparametric panel with time and fixed effects
#
# Update: April 24, 2025
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#############
# TODO
# 1) the case of dimension 2
# 2) estimate of time effects after washing out the FE
# 3) check if it more convenient washing out the TE before FE 
# 4) check why weights are not working in the optmization
# 

rm(list = ls())

library(sm)

#set.seed(124)

#Dimensions of the panel
N <- 50 #Number of units
T <- 10 #Number of periods

#Generation of exogenous
x <- rnorm(n=T*N, m=1,sd=1) 

#Generation of Fixed effects
FE <- x[1:N]+exp(rnorm(n=N,m=0,sd=0.5)) #Generation of fixed effects
#FE <- exp(rnorm(n=N,m=0,sd=0.5)) #Generation of fixed effects
FE <- FE - mean(FE) + 1 #ATTENTION: normalization of FE in order that the sum is equal to N, i.e. they are distributed around 1. This is an outcome of maximization process to be investigated

print(paste("Correlation between exogenous and FE = ", round(cor(x[1:N],FE),digits=3)))

#Time effetcs
TE <- seq(1:T)/T*2 #Generation of time effects

#Generation of random components
epsilon <- rnorm(n=T*N, m=0,sd=0.1) 

#Parameters of nonlinear relationship between exogenous and endogenous
beta.1 <- 1
beta.2 <- 0.2
beta.3 <- -0.02
beta.4 <- 0.001

FE.matrix <-  t(matrix(FE,ncol=N,nrow=T,byrow=T))
TE.matrix <-  t(matrix(TE,ncol=N,nrow=T))
x.matrix <- t(matrix(x,ncol=N,nrow=T))
epsilon.matrix <- t(matrix(epsilon,ncol=N,nrow=T))
  
y.matrix <- FE.matrix + TE.matrix + beta.1*x + beta.2*x^2 +beta.3*x^3 + beta.4*x^4  + epsilon.matrix

matplot(t(y.matrix),type="l")

# Estimate of FE via C-estimator ####
diff.y.array <- array(NA,dim=c(N,N,T))
diff.x.array  <- array(NA,dim=c(N,N,T))
for (j in (1:T)){
  for (i in (1:N)){
    diff.x <- abs(x.matrix[i,j]-x.matrix[,j]) #Distance between xs
    diff.x[i] <- max(diff.x) #Set the distance with itself to the maximum
    min.index <- which.min(diff.x) #Take the x with the minimum distance
    diff.x.array[i,min.index,j] <- diff.x[min.index] #Calculate the distance
    diff.y.array[i,min.index,j] <- y.matrix[i,j]- y.matrix[min.index,j] #Distance between ys 
  }
  diag(diff.y.array[,,j]) <- 0 #FE minus the same
  diag(diff.x.array[,,j]) <- NA #For avoiding issues in the weighting in the score function
}

#average.diff.alpha <- apply(diff.alpha.array,FUN=mean,MARGIN = c(1,2),na.rm=TRUE)
#average.diff.alpha[is.nan(average.diff.alpha)] <- NA
#diag(average.diff.alpha) <- 0

#Initial values of FE
ones <- matrix(1,nrow=N,ncol=1)

# score <- function(alpha.vec,A){
#   B <- alpha.vec%*%t(ones)
#   log(sum((B - t(B) - A)^2,na.rm = TRUE))
# }

score <- function(FE.vec,delta.Y,W){
  
  #alpha.vec <- alpha.ini
  #A <- diff.alpha.array
  #W <- diff.x.array
  
  #Vector of FE \times transposed vector of ones for a matrix with row with the same FE
  B <- FE.vec%*%t(ones)
  diff.B.array <- array(NA,dim=dim(delta.Y))
  for (j in 1:T){
    diff.B.array[,,j] <- B - t(B)
  }
  
  #Weights to be used in the calculation of score
  #wei <- ((1/W^2)/sum((1/W^2),na.rm = TRUE))*sum(!is.na(A))
  #wei <-  exp(( log(1/W)/prod(log(1/W), na.rm = TRUE))*sum(!is.na(delta.Y)))
  wei <- array(1,dim=dim(delta.Y))
  
  #Initialization of score matrix
  #C <- 0
  #
  #Calculation of score year by year and its cumulation over time. Use of weights inverse to the distance
  #for (j in 1:T){
  #  C <- C + sum( ((B - t(B) - delta.Y[,,j])^2)*(wei[, , j]), na.rm = TRUE)
  #}
  #log(C)
  
  log(sum( ((diff.B.array - delta.Y)^2)*(wei), na.rm = TRUE))
 
  }

## Initial values of FE estimation ####
FE.ini <- matrix(1,nrow=N,ncol=1)
#FE.ini <- FE

t0 <- Sys.time()
res <- optim(par = FE.ini, fn = score,delta.Y=diff.y.array,W=diff.x.array, method =  "L-BFGS-B")
print(Sys.time() - t0)

#The outcome of estimation
est.FE <- res$par

# Diagnostics of estimates
#print(cbind(est.FE,FE))
print(c(sum(FE),sum(est.FE))) #Sum of FE
print( mean((FE-est.FE)^2)) # MSE

#Explained variance of FE estimation
print(summary(lm(FE  ~ est.FE)))

#Plot of the esitmated FE
plot(est.FE,FE,pch=19,ylim=range(c(FE,est.FE)),xlim=range(c(FE,est.FE)),cex=0.75)
lines(c(-10,10),c(-10,10))
grid()

