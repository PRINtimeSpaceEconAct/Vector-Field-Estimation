# Clear workspace and load dependencies
setwd("~/Library/CloudStorage/OneDrive-UniversityofPisa/timeSpaceEvolutionEcAct/RVF/R code/Vector Field Estimation/")
rm(list = ls())
DEBUG = TRUE
source("src/libs/loadLib.R")
library(fields)


# parameters ----
nObs = 10000
nEval = 2500

# data Generation ----
X0 = mvrnorm(nObs, mu=c(0,0),Sigma = 0.5*diag(2))

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
#     return( -0.01*c(2*X[1], 2*X[2]) )
# }

# example 3 -- rotation ----
# M = matrix(c(cos(pi/2), -sin(pi/2), sin(pi/2), cos(pi/2)),nrow=2,ncol=2)
# VF <- function(X){
#     # X = (x,y), theta = pi/4
#     return (M %*% X)
# }

# apply VF
X1 = X0 + t(apply(X0, 1, VF))

# eval points
xGrid = seq(from=min(c(X0[,1],X1[,1])), to=max(c(X0[,1],X1[,1])), length.out=round(sqrt(nEval)))
yGrid = seq(from=min(c(X0[,2],X1[,2])), to=max(c(X0[,2],X1[,2])), length.out=round(sqrt(nEval)))
x = as.matrix(expand.grid(xGrid, yGrid))

# stima ----
t0 = Sys.time()
# est_field_adaptive = NWfield(X0, X1, x=x, kernel.type="epa",method.h = "sj",
#                                      chunk_size=1000,
#                                      sparse=FALSE, gc=TRUE)
est_field_adaptive = NWfieldAdaptive(X0, X1, x=x, kernel.type="epa",method.h = "sj",
                                     chunk_size=1000,
                                     sparse=FALSE, gc=TRUE, alpha=0.5)
est = est_field_adaptive
t = Sys.time() - t0

# plot ----
## plot campo vero ----
# dev.new()
VFx = t(apply(x, 1, VF))
plot(x, type = "n", xlab = "X1", ylab = "X2", main = "")
arrows(x[,1],x[,2],x[,1]+VFx[,1],x[,2]+VFx[,2],angle=15,col="black",length=0.05)
points(X0, col="red", pch=19,cex=0.01)
# dev.copy2pdf(file="testPics/doubleWellcampoVero.pdf")

## plot campo stimato ----
# dev.new()
signifVFest = significanceVF(est_field_adaptive,X0,X1,0.05)
plot(est_field_adaptive$x, type = "n", xlab = "X1", ylab = "X2", main = "")
arrows(est_field_adaptive$x[,1], est_field_adaptive$x[,2],
       est_field_adaptive$x[,1] + signifVFest$signif*est_field_adaptive$estimator[,1], 
       est_field_adaptive$x[,2] + signifVFest$signif*est_field_adaptive$estimator[,2],
       length = 0.05, angle = 15, col = "blue")
points(X0, col="red", pch=19,cex=0.01)
# dev.copy2pdf(file="testPics/doubleWellcampoStimato.pdf")


## plot errore ----
VFx = t(apply(x, 1, VF))
plot(est_field_adaptive$x, type = "n", xlab = "X", ylab = "Y", main = "Error Vector Field")
arrows(est_field_adaptive$x[,1], est_field_adaptive$x[,2],
       est_field_adaptive$x[,1] + est_field_adaptive$estimator[,1] - VFx[,1], 
       est_field_adaptive$x[,2] + est_field_adaptive$estimator[,2] - VFx[,2],
       length = 0.05, angle = 15, col = "red")

## image errore assoluto ----
errorNorm = sqrt((est_field_adaptive$estimator[,1] - VFx[,1])^2 + (est_field_adaptive$estimator[,2] - VFx[,2])^2)
image.plot(x = unique(x[,1]), y = unique(x[,2]), z = matrix(log10(errorNorm), nrow=sqrt(nEval), ncol=sqrt(nEval)),xlab="x",ylab="y",main="error norm (log10)")

## image errore relativo ----
dev.new()
errorNormRel = sqrt((est_field_adaptive$estimator[,1] - VFx[,1])^2 + (est_field_adaptive$estimator[,2] - VFx[,2])^2)/sqrt(VFx[,1]^2 + VFx[,2]^2)
image.plot(x = unique(x[,1]), y = unique(x[,2]), z = matrix(log10(errorNormRel), nrow=sqrt(nEval), ncol=sqrt(nEval)),xlab="X1",ylab="X2",main="")
points(X0, col="black", pch=19,cex=0.01)
dev.copy2pdf(file="testPics/LogdoubleWellNormRel.pdf")


# stima VAR ----
data = data.frame(X0 = X0[,1], Y0 = X0[,2], X1 = X1[,1], Y1 = X1[,2])
lm1 = lm(X1 ~ X0 + Y0, data = data)
lm2 = lm(Y1 ~ X0 + Y0, data = data)
A = matrix(c(coef(lm1)[2:3],coef(lm2)[2:3]),nrow=2,byrow=TRUE)
C = c(coef(lm1)[1],coef(lm2)[1])

VFVAR = function(x){
    return((A-diag(2)) %*% x + C)
}
VFVARx = t(apply(x, 1, VFVAR))

dev.new()
plot(x, type = "n", xlab = "X1", ylab ="X2", main = "")
arrows(x[,1], x[,2],
       x[,1] +  VFVARx[,1], 
       x[,2] +  VFVARx[,2],
       length = 0.05, angle = 15, col = "blue")
points(X0, col="red", pch=19,cex=0.01)
dev.copy2pdf(file="testPics/VARdoubleWellcampoStimato.pdf")


## image errore relativo ----
dev.new()
VFx = t(apply(x, 1, VF))
errorNormRel = sqrt((VFVARx[,1] - VFx[,1])^2 + (VFVARx[,2] - VFx[,2])^2)/sqrt(VFx[,1]^2 + VFx[,2]^2)
image.plot(x = unique(x[,1]), y = unique(x[,2]), z = matrix(log10(errorNormRel), nrow=sqrt(nEval), ncol=sqrt(nEval)),xlab="X1",ylab="X2",main="")
points(X0, col="black", pch=19,cex=0.01)
dev.copy2pdf(file="testPics/LogVARdoubleWellNormRel.pdf")

