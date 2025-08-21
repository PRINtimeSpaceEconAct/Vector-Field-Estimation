# Clear workspace and load dependencies
rm(list = ls())
DEBUG = FALSE
source("src/libs/loadLib.R")


# parameters ----
nObs = 10000
nEval = 1
x = matrix(c(0.5,0.5,-0.1,-0.1),nrow=2,ncol=2,byrow = T)
# data Generation ----
X0 = mvrnorm(nObs, mu=c(0,0),Sigma = 1*diag(2))
fx = dmvnorm(x, mu=c(0,0),Sigma = 1*diag(2))

A = matrix(runif(4),nrow=2,ncol=2)
SIGMA = t(A) %*% A
# SIGMA = 0.01^2 * diag(2)


# example 1 - double well ----
VF <- function(X){
    # X = (x,y)
    # U(X) = x^4 - x^2 + y^2
    # VF(X) = -grad U(X) = -(4x^3 - 2x, 2y)
    return( -0.01*c(4*X[1]^3 - 2*X[1], 2*X[2]) )
}

# example 2 -- single well ----
VF <- function(X){
    # X = (x,y)
    # U(X) = x^2 + y^2
    # VF(X) = -grad U(X) = -(2x, 2y)
    return( -0.01*c(2*X[1], 2*X[2]) )
}

# example 3 -- rotation ----
M = matrix(c(cos(pi/2), -sin(pi/2), sin(pi/2), cos(pi/2)),nrow=2,ncol=2)
VF <- function(X){
    # X = (x,y), theta = pi/4
    return (M %*% X)
}

# stima 1 ----
h = 0.14724848498449
X1 = X0 + t(apply(X0, 1, VF)) + mvrnorm(nObs, mu=c(0,0),Sigma = 1*diag(2)) %*% sqrtm(SIGMA)
t0 = Sys.time()
est_field = NWfield(X0, X1, x=x, kernel.type="gauss",h = 0.14724848498449,chunk_size=1000,gc=TRUE)
t = Sys.time() - t0

# stima MC
nMC = 1000
VFEstEval1 = matrix(nrow=nMC,ncol=2)
VFEstEval2 = matrix(nrow=nMC,ncol=2)

for(i in 1:nMC){
    print(i)
    X0 = mvrnorm(nObs, mu=c(0,0),Sigma = 1*diag(2))
    X1 = X0 + t(apply(X0, 1, VF)) + mvrnorm(nObs, mu=c(0,0),Sigma = 1*diag(2)) %*% sqrtm(SIGMA)    
    
    est_field = NWfield(X0, X1, x=x, kernel.type="gauss",h = 0.14724848498449,chunk_size=1000,gc=TRUE)    
    VFEstEval1[i,] = est_field$estimator[1,]
    VFEstEval2[i,] = est_field$estimator[2,]
}


# check variance covariance ----
k = 1/(2*sqrt(pi))

VFEval1 = VF(x[1,])
VFEval2 = VF(x[2,])

EstCovEv1 = cov(sqrt(nObs*h^2)*VFEstEval1)
EstCovEv2 = cov(sqrt(nObs*h^2)*VFEstEval2)

nu11 = sum(sqrtm(SIGMA)[1,1]^2 + sqrtm(SIGMA)[1,2]^2)
nu22 = sum(sqrtm(SIGMA)[2,1]^2 + sqrtm(SIGMA)[2,2]^2)
nu12 = sum(sqrtm(SIGMA)[1,1]*sqrtm(SIGMA)[2,1] + sqrtm(SIGMA)[1,2]*sqrtm(SIGMA)[2,2])
nu21 = nu12
NU = matrix(c(nu11,nu12,nu21,nu22),nrow=2,ncol=2)

TeoCovEv1 = NU/fx[1]*k^2
TeoCovEv2 = NU/fx[2]*k^2





