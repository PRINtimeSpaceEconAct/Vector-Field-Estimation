rm(list = ls())
setwd("~/Library/CloudStorage/OneDrive-UniversityofPisa/timeSpaceEvolutionEcAct/RVF/R code/Vector Field Estimation/")
source("src/libs/loadLib.R")
library(scales)
library(fields)
library(sm)
DEBUG = TRUE

# parameters ----
nObs = 200
nT = 10
nEval = 1024

# data Generation ----
set.seed(3)

# genero i FE
alpha_i = mvrnorm(nObs,mu=c(0,0),Sigma=0.01*diag(2))
alpha_i[,1] = alpha_i[,1] - sum(alpha_i[,1])/nObs
alpha_i[,2] = alpha_i[,2] - sum(alpha_i[,2])/nObs

# alpha_i[,1] = alpha_i[,1] - alpha_i[1,1]
# alpha_i[,2] = alpha_i[,2] - alpha_i[1,2]
# alpha_i[-1,1] = alpha_i[-1,1] - sum(alpha_i[,1])/(nObs-1)
# alpha_i[-1,2] = alpha_i[-1,2] - sum(alpha_i[,2])/(nObs-1)

# genero i TE
gamma_t = mvrnorm(nT,mu=c(0,0),Sigma=0.01*diag(2))
gamma_t[,1] = gamma_t[,1] - sum(gamma_t[,1])/nT
gamma_t[,2] = gamma_t[,2] - sum(gamma_t[,2])/nT
# gamma_t[,1] = gamma_t[,1] - gamma_t[1,1]
# gamma_t[,2] = gamma_t[,2] - gamma_t[1,2]
# gamma_t[-1,1] = gamma_t[-1,1] - sum(gamma_t[,1])/(nT-1)
# gamma_t[-1,2] = gamma_t[-1,2] - sum(gamma_t[,2])/(nT-1)


X0 = mvrnorm(nObs, mu=c(0,0),Sigma = 0.25*diag(2)) + alpha_i + array(rep(gamma_t[1,],each=nObs),dim=c(nObs,2))
X = array(NA,dim = c(nObs,2,nT+1))
X[,,1] = X0

# example 1 - double well ----
# VF <- function(X){
#     # X = (x,y)
#     # U(X) = x^4 - x^2 + y^2
#     # VF(X) = -grad U(X) = -(4x^3 - 2x, 2y)
#     return( -0.01*c(4*X[1]^3 - 2*X[1], 2*X[2]) )
# }
# 
# 
# JVF1 <- function(X){
#     # X = (x,y)
#     # U(X) = x^4 - x^2 + y^2
#     # VF(X) = -grad U(X) = -(4x^3 - 2x, 2y)
#     return( -0.01*c(12*X[1]^2 - 2, 0) )
# }
# 
# JVF2 <- function(X){
#     # X = (x,y)
#     # U(X) = x^4 - x^2 + y^2
#     # VF(X) = -grad U(X) = -(4x^3 - 2x, 2y)
#     return( c(0, 2) )
# }



# example 3 -- rotation ----
M = matrix(c(cos(pi/40), -sin(pi/40), sin(pi/40), cos(pi/40)),nrow=2,ncol=2)
VF <- function(X){
    # X = (x,y), theta = pi/4
    return ( (diag(2)-M) %*% (c(1,-1)-X) )
}

JVF1 <- function(X){
    return( M[1,]-diag(2)[1,] )
}

JVF2 <- function(X){
    return( M[2,]-diag(2)[2,] )
}

noise = array(rnorm(nObs*nT*2,sd = 0.01),dim=c(nObs,2,nT))
for (t in 1:nT){
    X[,,t+1] = X[,,t] + t(apply(X[,,t], 1, VF)) + alpha_i + array(rep(gamma_t[t,],each=nObs),dim=c(nObs,2)) + 
        + noise[,,t]
}

# plot a scatter plot of X at all time points
plot(X[,1,1], X[,2,1], type="p", col="red", pch=16,
     xlab="x", ylab="y", main="Scatter Plot of X")
points(X[,1,nT], X[,2,nT], type="p", col="blue", pch=16,
     xlab="x", ylab="y", main="Scatter Plot of X")

# eval points
xGrid = seq(from=min(X[,1,]), to=max(X[,1,]), length.out=round(sqrt(nEval)))
yGrid = seq(from=min(X[,2,]), to=max(X[,2,]), length.out=round(sqrt(nEval)))
x = as.matrix(expand.grid(xGrid, yGrid))
# x = rbind(x,(X[1,,1]))
# nEval = nEval + 1

Filtered = within_transform(X, FE = TRUE, TE = TRUE,
                            uniform_weights = TRUE, nEval_chunk = nEval,
                            x = x, kernel.type = "gauss",
                            method.h = "silverman", chunk_size = 512)


X0_raw = Filtered$X0_raw_unrolled
X0_star = Filtered$X0_star_unrolled
Y1_star = Filtered$Y1_star_unrolled
Y2_star = Filtered$Y2_star_unrolled
Y1 = Filtered$Y1_unrolled
Y2 = Filtered$Y2_unrolled


derivative_estimator_1 = compute_derivative_term(X0_raw, X0_star, x=x,
                                              kernel.type="gauss", D=NULL, 
                                              method.h="silverman", h=NULL, lambda=NULL, 
                                              sparse=FALSE, gc=FALSE, chunk_size=512, Y=Y1_star)

derivative_estimator_2 = compute_derivative_term(X0_raw, X0_star, x=x,
                                                 kernel.type="gauss", D=NULL,
                                                 method.h="silverman", h=NULL, lambda=NULL,
                                                 sparse=FALSE, gc=FALSE, chunk_size=512, Y=Y2_star)

derivative_obs_1 = compute_derivative_term(X0_raw, X_star=X0_star[,,rep(1,nrow(X0_raw))], x=X0_raw,
                                           kernel.type="gauss", D=NULL, 
                                           method.h="silverman", h=NULL, lambda=NULL, 
                                           sparse=FALSE, gc=FALSE, chunk_size=512, Y=Y1_star[,rep(1,nrow(Y1_star))])

derivative_obs_2 = compute_derivative_term(X0_raw, X_star=X0_star[,,rep(1,nrow(X0_raw))], x=X0_raw,
                                           kernel.type="gauss", D=NULL, 
                                           method.h="silverman", h=NULL, lambda=NULL, 
                                           sparse=FALSE, gc=FALSE, chunk_size=512, Y=Y2_star[,rep(1,nrow(Y2_star))])


meanPoint = apply(X0_raw, MARGIN = c(2), FUN = sum)/((nT-1)*nObs)
iBest = which.min(sqrt((x[,1]-meanPoint[1])^2 + (x[,2]-meanPoint[2])^2))


m10 = compute_m0(X_unrolled=X0_raw, Y_unrolled=Y1, beta=derivative_obs_1$estimator, x0=x[iBest,], beta_0=derivative_estimator_1$estimator[iBest,])
m20 = compute_m0(X_unrolled=X0_raw, Y_unrolled=Y2, beta=derivative_obs_2$estimator, x0=x[iBest,], beta_0=derivative_estimator_2$estimator[iBest,])

# VF_hat1 = compute_m(X0_raw, x, beta=derivative_estimator_1$estimator, m_0=VF(x[iBest,])[1], x0=x[iBest,], beta_0=derivative_estimator_1$estimator[iBest,])
# VF_hat2 = compute_m(X0_raw, x, beta=derivative_estimator_2$estimator, m_0=VF(x[iBest,])[2], x0=x[iBest,], beta_0=derivative_estimator_2$estimator[iBest,])
VF_hat1 = compute_m(X0_raw, x, beta=derivative_estimator_1$estimator, m_0=m10, x0=x[iBest,], beta_0=derivative_estimator_1$estimator[iBest,])
VF_hat2 = compute_m(X0_raw, x, beta=derivative_estimator_2$estimator, m_0=m20, x0=x[iBest,], beta_0=derivative_estimator_2$estimator[iBest,])

# Stitch together the two components of the vector field
VF_hat = cbind(VF_hat1, VF_hat2)
lengthArrows=1e-1
plot(x, type = "n", xlab = "X1", ylab="X2", main = " ")
arrows(x[,1],x[,2],x[,1]+lengthArrows*VF_hat[,1],x[,2]+lengthArrows*VF_hat[,2],angle=15,col="purple",length=0.05)
abline(h=0)
abline(v=0)
points(X[,1,],X[,2,])


VFx = t(apply(x, 1, VF))
ErrorVF = VFx - VF_hat
plot(x, type = "n", xlab = "X1", ylab="X2", main = "Error")
arrows(derivative_estimator_1$x[,1], derivative_estimator_1$x[,2],
       derivative_estimator_1$x[,1] + lengthArrows*ErrorVF[,1],
       derivative_estimator_1$x[,2] + lengthArrows*ErrorVF[,2],
       length = 0.05, angle = 15, col = "red")

est.dens = sm.density(X0_raw,display="none")


errorNorm = sqrt((VF_hat[,1] - VFx[,1])^2 + (VF_hat[,2] - VFx[,2])^2)/sqrt(VFx[,1]^2+VFx[,2]^2)
errorNorm[which.min(errorNorm)] = min(errorNorm[-c(which.min(errorNorm))])
image.plot(x = sort(unique(x[,1])), y = sort(unique(x[,2])), z = matrix(log10(errorNorm), nrow=sqrt(nEval), ncol=sqrt(nEval)),xlab="x",ylab="y",main="error norm rel (log10)")
contour(est.dens$eval.points[,1], est.dens$eval.points[,2], est.dens$estimate,add=T)

errorAbs = sqrt((VF_hat[,1] - VFx[,1])^2 + (VF_hat[,2] - VFx[,2])^2)
errorAbs[which.min(errorAbs)] = min(errorAbs[-c(which.min(errorAbs))])
image.plot(x = sort(unique(x[,1])), y = sort(unique(x[,2])), z = matrix(log10(errorAbs), nrow=sqrt(nEval), ncol=sqrt(nEval)),xlab="x",ylab="y",main="error norm abs (log10)")
contour(est.dens$eval.points[,1], est.dens$eval.points[,2], est.dens$estimate,add=T)

# reconstruct FE
VF_hat1 = compute_m(X0_raw, X0_raw, beta=derivative_obs_1$estimator, m_0=m10, x0=x[iBest,], beta_0=derivative_estimator_1$estimator[iBest,])
VF_hat2 = compute_m(X0_raw, X0_raw, beta=derivative_obs_2$estimator, m_0=m20, x0=x[iBest,], beta_0=derivative_estimator_2$estimator[iBest,])

 
YObs = aperm(array(cbind(Y1,Y2),dim = c(nT,nObs,2)),c(2,3,1))
VFObs = aperm(array(cbind(VF_hat1,VF_hat2),dim = c(nT,nObs,2)),c(2,3,1))

alpha_i_hat = apply(YObs - VFObs, MARGIN = c(1, 2), FUN = sum)/nT
gamma_t_hat = t(apply(YObs - VFObs, MARGIN = c(2, 3), FUN = sum)/nObs)

# test con regressione
lm_alpha = lm(c(alpha_i) ~ c(alpha_i_hat) + 0); summary(lm_alpha)
lm_gamma = lm(c(gamma_t) ~ c(gamma_t_hat) + 0); summary(lm_gamma)





