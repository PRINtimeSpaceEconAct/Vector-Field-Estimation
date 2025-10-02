# Clear workspace and load dependencies
setwd("~/Library/CloudStorage/OneDrive-UniversityofPisa/timeSpaceEvolutionEcAct/RVF/R code/Vector Field Estimation/")
rm(list = ls())
DEBUG = FALSE
source("src/libs/loadLib.R")
library(fields)
library(latex2exp)


# parameters ----
nObs = 1000
nEval = 2500

# data Generation ----
set.seed(1)
X0 = mvrnorm(nObs, mu=c(0,0),Sigma = 0.5*diag(2))

# example 1 - double well ----
VF <- function(X){
    # X = (x,y)
    # U(X) = x^4 - x^2 + y^2
    # VF(X) = -grad U(X) = -(4x^3 - 2x, 2y)
    return( -0.1*c(4*X[1]^3 - 2*X[1], 2*X[2]) )
}

# example 2 -- single well ----
# VF <- function(X){
#     # X = (x,y)
#     # U(X) = x^2 + y^2
#     # VF(X) = -grad U(X) = -(2x, 2y)
#     return( -0.01*c(2*X[1], 2*X[2]) )
# }

# example 3 -- rotation ----
# M = matrix(c(cos(pi/2), -sin(pi/2), sin(pi/2), cos(pi/2)),nrow=2,ncol=2)
# VF <- function(X){
#     # X = (x,y), theta = pi/4
#     return (M %*% X)
# }

# apply VF
X1 = X0 + t(apply(X0, 1, VF)) +  matrix(rnorm(2*nObs),nrow=nObs) %*% matrix(c(0.01,0.005,0.005,0.02),nrow=2)

# eval points
X_unrolled <- rbind(X0, X1)
x <- defineEvalPoints(X_unrolled, nEval)

# stima ----

est_field_adaptiveLL = LLfieldAdaptive(X0, X1, x=x, kernel.type="epa",method.h = "silverman",
                                     chunk_size=3000,
                                     gc=TRUE, hOpt = TRUE, alphaOpt = TRUE)

# est_field_adaptiveNW = NWfieldAdaptive(X0, X1, x=x, kernel.type="epa",method.h = "silverman",
#                                      chunk_size=1000,
#                                      gc=TRUE,
#                                      hOpt = TRUE, alphaOpt = TRUE)
# est_field = NWfield(X0, X1, x=x, kernel.type="gauss",h = 0.14724848498449,chunk_size=1000,gc=TRUE)
# save(est_field_adaptiveLL, file = "est_field_adaptiveLL.RData")

# bootstrap
# est_field_adaptive_bootstrap = bootstrapKernelFieldErrors(est_field_adaptive, B = 10)
# signifBoot = significanceBootstrap(est_field_adaptive,est_field_adaptive_bootstrap)

# est = est_field_adaptive


# plot ----
# plot campo vero ----
dev.new()
op <- par(family = "mono") #Possible families: "mono", "Helvetica","Palatino" or "Times"
VFx = t(apply(x, 1, VF))
lengthArrows=0.1
plot(x, type = "n", xlab = TeX(r'($X_1$)'), ylab=TeX(r'($X_2$)'), main = "")
arrows(x[,1],x[,2],x[,1]+lengthArrows*VFx[,1],x[,2]+lengthArrows*VFx[,2],angle=15,col="black",length=0.05)
abline(h=0)
abline(v=0)
dev.copy2pdf(file="testPics/doubleWellcampoVero.pdf")
# par(op)

## plot campo stimato ----
dev.new()
op <- par(family = "mono") #Possible families: "mono", "Helvetica","Palatino" or "Times" 
signifVFest = significanceVF(est_field_adaptiveLL)
lengthArrows=0.1
plot(x, type = "n", xlab = TeX(r'($X_1$)'), ylab=TeX(r'($X_2$)'), main = "Tutte frecce")
arrows(est_field_adaptiveLL$x[,1], est_field_adaptiveLL$x[,2],
       est_field_adaptiveLL$x[,1] + lengthArrows*est_field_adaptiveLL$estimator[,1], 
       est_field_adaptiveLL$x[,2] + lengthArrows*est_field_adaptiveLL$estimator[,2],
       length = 0.05, angle = 15, col = "blue")
abline(h=0)
abline(v=0)
dev.copy2pdf(file="testPics/doubleWellcampoStimato.pdf")
par(op)

## plot campo stimato significativo ----
signifInd = which(signifBoot)
dev.new()
op <- par(family = "mono") #Possible families: "mono", "Helvetica","Palatino" or "Times" 
lengthArrows=0.1
plot(x, type = "n", xlab = TeX(r'($X_1$)'), ylab=TeX(r'($X_2$)'), main = "Solo significative")
arrows(est_field_adaptiveLL$x[signifInd,1], est_field_adaptiveLL$x[signifInd,2],
       est_field_adaptiveLL$x[signifInd,1] + lengthArrows*est_field_adaptiveLL$estimator[signifInd,1], 
       est_field_adaptiveLL$x[signifInd,2] + lengthArrows*est_field_adaptiveLL$estimator[signifInd,2],
       length = 0.05, angle = 15, col = "blue")
abline(h=0)
abline(v=0)
# dev.copy2pdf(file="testPics/doubleWellcampoStimato.pdf")
par(op)


## plot errore ----
VFx = t(apply(x, 1, VF))
plot(est_field_adaptiveLL$x, type = "n", xlab = "X", ylab = "Y", main = "Error Vector Field")
arrows(est_field_adaptiveLL$x[,1], est_field_adaptiveLL$x[,2],
       est_field_adaptiveLL$x[,1] + est_field_adaptiveLL$estimator[,1] - VFx[,1],
       est_field_adaptiveLL$x[,2] + est_field_adaptiveLL$estimator[,2] - VFx[,2],
       length = 0.05, angle = 15, col = "red")

## image errore assoluto ----
errorNorm = sqrt((est_field_adaptiveLL$estimator[,1] - VFx[,1])^2 + (est_field_adaptiveLL$estimator[,2] - VFx[,2])^2)
image.plot(x = unique(x[,1]), y = unique(x[,2]), z = matrix(log10(errorNorm), nrow=sqrt(nEval), ncol=sqrt(nEval)),xlab="x",ylab="y",main="error norm (log10)")

## image errore relativo ----
VFx = t(apply(x, 1, VF))
dev.new()
op <- par(family = "mono") #Possible families: "mono", "Helvetica","Palatino" or "Times" 
errorNormRel = sqrt((est_field_adaptiveLL$estimator[,1] - VFx[,1])^2 + (est_field_adaptiveLL$estimator[,2] - VFx[,2])^2)/sqrt(VFx[,1]^2 + VFx[,2]^2)
errorNormRel[is.na(errorNormRel)] = 1; errorNormRel[!is.finite(errorNormRel)] = 1; errorNormRel[errorNormRel==0] = 1
image.plot(x = unique(x[,1]), y = unique(x[,2]), z = matrix(log10(errorNormRel), nrow=sqrt(nEval), ncol=sqrt(nEval)),xlab=TeX(r'($X_1$)'),ylab=TeX(r'($X_2$)'),main="")
abline(h=0)
abline(v=0)
dev.copy2pdf(file="testPics/LogdoubleWellNormRelAdaptive.pdf")
par(op)


# stima VAR ----
data = data.frame(X0 = X0[,1], Y0 = X0[,2], X1 = X1[,1], Y1 = X1[,2])
lm1 = lm(X1 ~ X0 + Y0 + 0, data = data)
lm2 = lm(Y1 ~ X0 + Y0 + 0, data = data)
A = matrix(c(coef(lm1)[1:2],coef(lm2)[1:2]),nrow=2,byrow=TRUE)

VFVAR = function(x){
    return((A-diag(2)) %*% x)
}
VFVARx = t(apply(x, 1, VFVAR))

dev.new()
op <- par(family = "mono") #Possible families: "mono", "Helvetica","Palatino" or "Times" 
plot(x, type = "n", xlab = "X1", ylab ="X2", main = "")
arrows(x[,1], x[,2],
       x[,1] +  VFVARx[,1], 
       x[,2] +  VFVARx[,2],
       length = 0.05, angle = 15, col = "blue")
abline(h=0)
abline(v=0)
dev.copy2pdf(file="testPics/VARdoubleWellcampoStimato.pdf")
par(op)

## image errore relativo ----
dev.new()
op <- par(family = "mono") #Possible families: "mono", "Helvetica","Palatino" or "Times" 
VFx = t(apply(x, 1, VF))
errorNormRel = sqrt((VFVARx[,1] - VFx[,1])^2 + (VFVARx[,2] - VFx[,2])^2)/sqrt(VFx[,1]^2 + VFx[,2]^2)
image.plot(x = unique(x[,1]), y = unique(x[,2]), z = matrix(log10(errorNormRel), nrow=sqrt(nEval), ncol=sqrt(nEval)),xlab=TeX(r'($X_1$)'),ylab=TeX(r'($X_2$)'),main="")
abline(h=0)
abline(v=0)
dev.copy2pdf(file="testPics/LogVARdoubleWellNormRel.pdf")
par(op)


