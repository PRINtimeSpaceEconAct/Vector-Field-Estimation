# Clear workspace and load dependencies
rm(list = ls())
DEBUG = TRUE
source("src/libs/loadLib.R")


# parameters ----
nObs = 10000
nEval = 2500

# data Generation ----
X0 = mvrnorm(nObs, mu=c(0,0),Sigma = 1*diag(2))

# eval points
xGrid = seq(from=-1, to=1, length.out=round(sqrt(nEval)))
yGrid = seq(from=-1, to=1, length.out=round(sqrt(nEval)))
x = as.matrix(expand.grid(xGrid, yGrid))

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

# stima ----
X1 = X0 + t(apply(X0, 1, VF))
t0 = Sys.time()
est_field_adaptive = NWfieldAdaptive(X0, X1, x=x, kernel.type="epa",method.h = "sj",
                                     chunk_size=1000,
                                     sparse=FALSE, gc=TRUE, alpha=0.5)
# est_field_adaptive = LLfieldAdaptive(X0, X1, x=x, kernel.type="epa",method.h = "sj",
                                     # chunk_size=1000,
                                     # sparse=FALSE, gc=TRUE, alpha=0.5)

t = Sys.time() - t0

# plot ----
## plot campo vero ----
VFx = t(apply(x, 1, VF))
plot(x, type = "n", xlab = "X", ylab = "Y", main = "True Vector Field")
arrows(x[,1],x[,2],x[,1]+VFx[,1],x[,2]+VFx[,2],angle=15,col="black",length=0.05)


## plot campo stimato ----
plot(est_field_adaptive$x, type = "n", xlab = "X", ylab = "Y", main = "Estimated Vector Field")
arrows(est_field_adaptive$x[,1], est_field_adaptive$x[,2],
       est_field_adaptive$x[,1] + est_field_adaptive$estimator[,1], 
       est_field_adaptive$x[,2] + est_field_adaptive$estimator[,2],
       length = 0.05, angle = 15, col = "blue")

## plot errore ----
VFx = t(apply(x, 1, VF))
plot(est_field_adaptive$x, type = "n", xlab = "X", ylab = "Y", main = "Error Vector Field")
arrows(est_field_adaptive$x[,1], est_field_adaptive$x[,2],
       est_field_adaptive$x[,1] + est_field_adaptive$estimator[,1] - VFx[,1], 
       est_field_adaptive$x[,2] + est_field_adaptive$estimator[,2] - VFx[,2],
       length = 0.05, angle = 15, col = "red")

## image errore ----
library(fields)
errorNorm = sqrt((est_field_adaptive$estimator[,1] - VFx[,1])^2 + (est_field_adaptive$estimator[,2] - VFx[,2])^2)
image.plot(x = unique(x[,1]), y = unique(x[,2]), z = matrix(log10(errorNorm), nrow=sqrt(nEval), ncol=sqrt(nEval)),xlab="x",ylab="y",main="error norm (log10)")


