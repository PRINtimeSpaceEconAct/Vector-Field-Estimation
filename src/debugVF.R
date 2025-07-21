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

# apply VF
X1 = X0 + t(apply(X0, 1, VF)) +  matrix(rnorm(2*nObs),nrow=nObs) %*% matrix(c(0.01,0.005,0.005,0.02),nrow=2)

# eval points
xGrid = seq(from=-1, to=1, length.out=round(sqrt(nEval)))
yGrid = seq(from=-1, to=1, length.out=round(sqrt(nEval)))
# xGrid = seq(from=min(c(X0[,1],X1[,1])), to=max(c(X0[,1],X1[,1])), length.out=round(sqrt(nEval)))
# yGrid = seq(from=min(c(X0[,2],X1[,2])), to=max(c(X0[,2],X1[,2])), length.out=round(sqrt(nEval)))
x = as.matrix(expand.grid(xGrid, yGrid))

# stima ----

# est_field_adaptive = LLfieldAdaptive(X0, X1, x=x, kernel.type="epa",
#                                      chunk_size=3000,
#                                      sparse=FALSE, gc=TRUE, h = 0.1519322, alpha = 0.225)

# est_field_adaptive = LLfieldAdaptive(X0, X1, x=x, kernel.type="epa",
#                                      chunk_size=3000,
#                                      sparse=FALSE, gc=TRUE, h = 0.1, alpha = 0.225)

# est_field_adaptiveLL = LLfield(X0, X1, x=x, kernel.type="epa",
#                                      chunk_size=3000,
#                                      sparse=FALSE, gc=TRUE, h = 0.1)
# 
# est_field_adaptiveNW = NWfield(X0, X1, x=x, kernel.type="epa",
#                              chunk_size=3000,
#                              sparse=FALSE, gc=TRUE, h = 0.1)

est_field_adaptiveLL = LLfieldAdaptive(X0, X1, x=x, kernel.type="epa",method.h = "silverman",
                                     chunk_size=3000,
                                     sparse=FALSE, gc=TRUE, hOpt = TRUE, alphaOpt = TRUE)

est_field_adaptiveNW = NWfieldAdaptive(X0, X1, x=x, kernel.type="epa",method.h = "silverman",
                                     chunk_size=1000,
                                     sparse=FALSE, gc=TRUE,
                                     hOpt = TRUE, alphaOpt = TRUE)



# dev.new()
op <- par(family = "mono") 
lengthArrows=1.0
VFx = t(apply(x, 1, VF))
plot(x, type = "n", xlab = TeX(r'($X_1$)'), ylab=TeX(r'($X_2$)'), main = "Tutte frecce")
arrows(x[,1],x[,2],x[,1]+lengthArrows*VFx[,1],x[,2]+lengthArrows*VFx[,2],angle=15,col="black",length=0.05)
arrows(est_field_adaptiveLL$x[,1], est_field_adaptiveLL$x[,2],
       est_field_adaptiveLL$x[,1] + lengthArrows*est_field_adaptiveLL$estimator[,1], 
       est_field_adaptiveLL$x[,2] + lengthArrows*est_field_adaptiveLL$estimator[,2],
       length = 0.05, angle = 15, col = "blue")
arrows(est_field_adaptiveNW$x[,1], est_field_adaptiveNW$x[,2],
       est_field_adaptiveNW$x[,1] + lengthArrows*est_field_adaptiveNW$estimator[,1], 
       est_field_adaptiveNW$x[,2] + lengthArrows*est_field_adaptiveNW$estimator[,2],
       length = 0.05, angle = 15, col = "red")
abline(h=0)
abline(v=0)
par(op)

errorNormNW = sqrt((est_field_adaptiveNW$estimator[,1] - VFx[,1])^2 + (est_field_adaptiveNW$estimator[,2] - VFx[,2])^2)
image.plot(x = unique(x[,1]), y = unique(x[,2]), z = matrix(log10(errorNormNW), nrow=sqrt(nEval), ncol=sqrt(nEval)),xlab="x",ylab="y",main="error norm NW (log10)")

errorNormLL = sqrt((est_field_adaptiveLL$estimator[,1] - VFx[,1])^2 + (est_field_adaptiveLL$estimator[,2] - VFx[,2])^2)
image.plot(x = unique(x[,1]), y = unique(x[,2]), z = matrix(log10(errorNormLL), nrow=sqrt(nEval), ncol=sqrt(nEval)),xlab="x",ylab="y",main="error norm LL (log10)")

