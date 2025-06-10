rm(list = ls())

# parameters ----
nObs = 100
nT = 10

# data Generation ----
set.seed(1)

# genero i FE
alpha_i = mvrnorm(nObs,mu=c(0,0),Sigma=0.01*diag(2))
alpha_i[,1] = alpha_i[,1] - sum(alpha_i[,1])/nObs
alpha_i[,2] = alpha_i[,2] - sum(alpha_i[,2])/nObs

# genero i TE
gamma_t = mvrnorm(nT,mu=c(0,0),Sigma=0.01*diag(2))
gamma_t[,1] = gamma_t[,1] - sum(gamma_t[,1])/nT
gamma_t[,2] = gamma_t[,2] - sum(gamma_t[,2])/nT

X0 = mvrnorm(nObs, mu=c(0,0),Sigma = 0.25*diag(2)) + alpha_i + array(rep(gamma_t[1,],each=nObs),dim=c(nObs,2))
X = array(NA,dim = c(nObs,2,nT+1))
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


for (t in 1:nT){
    X[,,t+1] = X[,,t] + t(apply(X[,,t], 1, VF)) + alpha_i + array(rep(gamma_t[t,],each=nObs),dim=c(nObs,2)) + 
        + mvrnorm(nObs, mu=c(0,0),Sigma = 0.001*diag(2))
}

# calcolo i Delta
Delta = (X[,,2:(nT+1)] - X[,,1:nT])

