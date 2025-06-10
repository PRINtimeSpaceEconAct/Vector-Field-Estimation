# Clear workspace and load dependencies
setwd("~/Library/CloudStorage/OneDrive-UniversityofPisa/timeSpaceEvolutionEcAct/RVF/R code/Vector Field Estimation/")
rm(list = ls())
graphics.off()
dev.new()
DEBUG = FALSE
source("src/libs/loadLib.R")
library(fields)
library(latex2exp)


# parameters ----
nObs = 100
nT = 10
nEval = 2500

# eval points
xGrid = seq(from=-1, to=1, length.out=round(sqrt(nEval)))
yGrid = seq(from=-1, to=1, length.out=round(sqrt(nEval)))
x = as.matrix(expand.grid(xGrid, yGrid))

# data Generation ----
set.seed(1)

# genero i FE
alpha_i = mvrnorm(nObs,mu=c(0,0),Sigma=0.01*diag(2))
alpha_i[,1] = alpha_i[,1] - sum(alpha_i[,1])/nObs
alpha_i[,2] = alpha_i[,2] - sum(alpha_i[,2])/nObs
# alpha_i = matrix(0,nrow=nObs,ncol=2)

# genero i TE
gamma_t = mvrnorm(nT,mu=c(0,0),Sigma=0.01*diag(2))
gamma_t[,1] = gamma_t[,1] - sum(gamma_t[,1])/nT
gamma_t[,2] = gamma_t[,2] - sum(gamma_t[,2])/nT
# gamma_t = matrix(0,nrow = nT, ncol = 2)

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

# example 2 -- single well ----
# VF <- function(X){
#     # X = (x,y)
#     # U(X) = x^2 + y^2
#     # VF(X) = -grad U(X) = -(2x, 2y)
#     return( -0.05*c(2*X[1], 2*X[2]) )
# }

# example 3 -- rotation ----
M = matrix(c(cos(pi/2), -sin(pi/2), sin(pi/2), cos(pi/2)),nrow=2,ncol=2)
VF <- function(X){
    # X = (x,y), theta = pi/4
    return (0.1*(M %*% X - X))
}


for (t in 1:nT){
    X[,,t+1] = X[,,t] + t(apply(X[,,t], 1, VF)) + alpha_i + array(rep(gamma_t[t,],each=nObs),dim=c(nObs,2)) + 
        + mvrnorm(nObs, mu=c(0,0),Sigma = 0.001*diag(2))
}

# calcolo i Delta
Delta = (X[,,2:(nT+1)] - X[,,1:nT])

# demeaning per i FE
MeanDelta_i = apply(Delta, c(1, 2), mean)
MeanDelta_t = apply(Delta, c(2, 3), mean)
DeltaFE = Delta - array(rep(MeanDelta_i, nT), dim = c(nObs, 2, nT))  +
    - array(rep(MeanDelta_t, each = nObs), dim = c(nObs, 2, nT))  +
     + aperm(array(rep(apply(Delta, 2, mean),nObs*nT),dim = c(2,nObs,nT)),c(2,1,3))

# DeltaFE = Delta

# metto in formato input per il VF
DeltaFEStack = do.call(rbind, lapply(1:nT, function(t) DeltaFE[,,t]))
X0Stack = do.call(rbind, lapply(1:nT, function(t) X[,,t]))
X1Stack = X0Stack + DeltaFEStack

# diagnostic on means
MX = array(NA, dim = c(nObs, 2, nT))
for (t in 1:nT){
    MX[,,t] = t(apply(X[,,t], 1, VF))
}
meanMX_i = apply(MX, c(1, 2), mean)
meanMX_t = apply(MX, c(2, 3), mean)
summary(meanMX_i)
summary(t(meanMX_t))

# stima ----
# est_field_adaptive = NWfieldAdaptive(X0Stack, X1Stack, x=x, kernel.type="gauss", hOpt = FALSE, alphaOpt = FALSE, chunk_size = 1000)

est_field_adaptive = NWfieldAdaptive(X0Stack, X1Stack, x=x, kernel.type="epa",
                    hOpt = TRUE, alphaOpt = TRUE)

# bootstrap
# est_field_adaptive_bootstrap = bootstrapKernelFieldErrors(est_field_adaptive, B = 10)
# signifBoot = significanceBootstrap(est_field_adaptive,est_field_adaptive_bootstrap)

# est = est_field_adaptive

# plot ----
# plot campo vero ----
# dev.new()
# op <- par(family = "mono") #Possible families: "mono", "Helvetica","Palatino" or "Times"
# VFx = t(apply(x, 1, VF))
# lengthArrows=1.0
# plot(x, type = "n", xlab = TeX(r'($X_1$)'), ylab=TeX(r'($X_2$)'), main = "")
# arrows(x[,1],x[,2],x[,1]+lengthArrows*VFx[,1],x[,2]+lengthArrows*VFx[,2],angle=15,col="black",length=0.05)
# abline(h=0)
# abline(v=0)
# # dev.copy2pdf(file="testPics/doubleWellcampoVero.pdf")
# par(op)

## plot campo stimato ----
dev.new()
op <- par(family = "mono") #Possible families: "mono", "Helvetica","Palatino" or "Times"
signifVFest = significanceVF(est_field_adaptive)
lengthArrows=1.0
plot(x, type = "n", xlab = TeX(r'($X_1$)'), ylab=TeX(r'($X_2$)'), main = "Tutte frecce")
arrows(est_field_adaptive$x[,1], est_field_adaptive$x[,2],
       est_field_adaptive$x[,1] + lengthArrows*est_field_adaptive$estimator[,1],
       est_field_adaptive$x[,2] + lengthArrows*est_field_adaptive$estimator[,2],
       length = 0.05, angle = 15, col = "blue")
points(X0Stack)
abline(h=0)
abline(v=0)
# dev.copy2pdf(file="testPics/doubleWellcampoStimato.pdf")
par(op)

## plot errore ----
dev.new()
VFx = t(apply(x, 1, VF))
plot(est_field_adaptive$x, type = "n", xlab = "X", ylab = "Y", main = "Error Vector Field")
arrows(est_field_adaptive$x[,1], est_field_adaptive$x[,2],
       est_field_adaptive$x[,1] + est_field_adaptive$estimator[,1] - VFx[,1],
       est_field_adaptive$x[,2] + est_field_adaptive$estimator[,2] - VFx[,2],
       length = 0.05, angle = 15, col = "red")

## image errore realtivo ----
dev.new()
errorNorm = sqrt((est_field_adaptive$estimator[,1] - VFx[,1])^2 + (est_field_adaptive$estimator[,2] - VFx[,2])^2)/sqrt(VFx[,1]^2+VFx[,2]^2)
image.plot(x = unique(x[,1]), y = unique(x[,2]), z = matrix(log10(errorNorm), nrow=sqrt(nEval), ncol=sqrt(nEval)),xlab="x",ylab="y",main="error norm rel (log10)")


# stima dei FE
est_field_adaptiveOBS = NWfieldAdaptive(X0Stack, X1Stack, x=X0Stack, kernel.type="epa",
                                      hOpt = TRUE, alphaOpt = TRUE)

# est_field_adaptiveOBS = NWfieldAdaptive(X0Stack, X1Stack, x=X0Stack, kernel.type="gauss", hOpt = FALSE, alphaOpt = FALSE, chunk_size = 1000)
                                     
DeltaEstStack = est_field_adaptiveOBS$estimator
DeltaEst1 = array(DeltaEstStack[,1], dim = c(nObs,nT))
DeltaEst2 = array(DeltaEstStack[,2], dim = c(nObs,nT))
DeltaEst = aperm(array(c(DeltaEst1,DeltaEst2),dim=c(nObs,nT,2)),c(1,3,2))

epsi1 = rowMeans((Delta[,1,]-DeltaEst[,1,]),na.rm=T) 
epsi2 = rowMeans((Delta[,2,]-DeltaEst[,2,]),na.rm=T)
epst1 = colMeans((Delta[,1,]-DeltaEst[,1,]),na.rm=T) 
epst2 = colMeans((Delta[,2,]-DeltaEst[,2,]),na.rm=T)
eps1 = mean((Delta[,1,]-DeltaEst[,1,]),na.rm=T)
eps2 = mean((Delta[,2,]-DeltaEst[,2,]),na.rm=T)

FE1 = epsi1-eps1
FE2 = epsi2-eps2
TE1 = epst1-eps1
TE2 = epst2-eps2

dev.new()
par(mfrow=c(2,2))
plot(alpha_i[,1],FE1,xlab="FE1",ylab="alpha1",asp=1,main="FE1")
abline(lm(FE1~alpha_i[,1]),col="red")
plot(alpha_i[,2],FE2,xlab="FE2",ylab="alpha2",asp=1,main="FE2")
abline(lm(FE2~alpha_i[,2]),col="red")
plot(gamma_t[,1],TE1,xlab="TE1",ylab="gamma1",asp=1,main = "TE1")
abline(lm(TE1~gamma_t[,1]),col="red")
plot(gamma_t[,2],TE2,xlab="TE2",ylab="gamma2",asp=1,main = "TE2")
abline(lm(TE2~gamma_t[,2]),col="red")


c(coef(lm(FE1~alpha_i[,1]))[2], summary(lm(FE1~alpha_i[,1]))$r.squared)
c(coef(lm(FE2~alpha_i[,2]))[2], summary(lm(FE2~alpha_i[,2]))$r.squared)
c(coef(lm(TE1~gamma_t[,1]))[2], summary(lm(TE1~gamma_t[,1]))$r.squared)
c(coef(lm(TE2~gamma_t[,2]))[2], summary(lm(TE2~gamma_t[,2]))$r.squared)

# provo modello stimato con i FE contro modello vero
X1StackFEHat = X0Stack + est_field_adaptiveOBS$estimator + 
    + array(rep(cbind(FE1,FE2), nT), dim = c(nObs*nT, 2))  +
    + array(rep(cbind(TE1,TE2), each = nObs), dim = c(nObs*nT, 2)) 

summary((X1StackFEHat - X1Stack))
