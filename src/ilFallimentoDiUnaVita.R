rm(list = ls())
setwd("~/Library/CloudStorage/OneDrive-UniversityofPisa/timeSpaceEvolutionEcAct/RVF/R code/Vector Field Estimation/")
source("src/libs/loadLib.R")
DEBUG = TRUE

# parameters ----
nObs = 1000
nT = 10
nEval = 1024

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

X0 = mvrnorm(nObs, mu=c(0,0),Sigma = 0.25*diag(2)) #+ alpha_i + array(rep(gamma_t[1,],each=nObs),dim=c(nObs,2))
X = array(NA,dim = c(nObs,2,nT+1))
X[,,1] = X0

# example 1 - double well ----
VF <- function(X){
    # X = (x,y)
    # U(X) = x^4 - x^2 + y^2
    # VF(X) = -grad U(X) = -(4x^3 - 2x, 2y)
    return( -0.01*c(4*X[1]^3 - 2*X[1], 2*X[2]) )
}


JVF1 <- function(X){
    # X = (x,y)
    # U(X) = x^4 - x^2 + y^2
    # VF(X) = -grad U(X) = -(4x^3 - 2x, 2y)
    return( -0.01*c(12*X[1]^2 - 2, 0) )
}

JVF2 <- function(X){
    # X = (x,y)
    # U(X) = x^4 - x^2 + y^2
    # VF(X) = -grad U(X) = -(4x^3 - 2x, 2y)
    return( c(0, 2) )
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
    X[,,t+1] = X[,,t] + t(apply(X[,,t], 1, VF)) #+ alpha_i + array(rep(gamma_t[t,],each=nObs),dim=c(nObs,2)) + 
        + mvrnorm(nObs, mu=c(0,0),Sigma = 0.001*diag(2))
}
# plot a scatter plot of X at all time points
plot(X[,1,1], X[,2,1], type="p", col="red", pch=16,
     xlab="x", ylab="y", main="Scatter Plot of X")
points(X[,1,10], X[,2,10], type="p", col="blue", pch=16,
     xlab="x", ylab="y", main="Scatter Plot of X")

# eval points
xGrid = seq(from=min(X), to=max(X), length.out=round(sqrt(nEval)))
yGrid = seq(from=min(X), to=max(X), length.out=round(sqrt(nEval)))
x = as.matrix(expand.grid(xGrid, yGrid))


Filtered = within_transform(X, FE = TRUE, TE = TRUE,
                            uniform_weights = FALSE, nEval_chunk = nEval,
                            x = x, kernel.type = "gauss",
                            method.h = "silverman", chunk_size = 512)


X0_raw = Filtered$X0_raw_unrolled
X0_star = Filtered$X0_star_unrolled
Y1 = Filtered$Y1_unrolled
Y2 = Filtered$Y2_unrolled

derivative_estimator = compute_derivative_term(X0_raw, X0_star, x=x,
                                              kernel.type="gauss", D=NULL, 
                                              method.h="silverman", h=NULL, lambda=NULL, 
                                              sparse=FALSE, gc=FALSE, chunk_size=512, Y=Y1)

# plot true VF
VFx = t(apply(x, 1, VF))
lengthArrows=0.1
plot(x, type = "n", xlab = "X1", ylab="X2", main = " ")
arrows(x[,1],x[,2],x[,1]+lengthArrows*VFx[,1],x[,2]+lengthArrows*VFx[,2],angle=15,col="black",length=0.05)
abline(h=0)
abline(v=0)

# plot the derivative estimator
lengthArrows=1.0
JVF1x = t(apply(derivative_estimator$x, 1, JVF1))
plot(x, type = "n", xlab = "X1", ylab="X2", main = "Tutte frecce")
arrows(derivative_estimator$x[,1], derivative_estimator$x[,2],
       derivative_estimator$x[,1] + lengthArrows*JVF1x[,1],
       derivative_estimator$x[,2] + lengthArrows*JVF1x[,2],
       length = 0.05, angle = 15, col = "black")

arrows(derivative_estimator$x[,1], derivative_estimator$x[,2],
       derivative_estimator$x[,1] + lengthArrows*derivative_estimator$estimator[,1],
       derivative_estimator$x[,2] + lengthArrows*derivative_estimator$estimator[,2],
       length = 0.05, angle = 15, col = "blue")


ErrorVF1 = JVF1x - derivative_estimator$estimator
plot(x, type = "n", xlab = "X1", ylab="X2", main = "Error")
arrows(derivative_estimator$x[,1], derivative_estimator$x[,2],
       derivative_estimator$x[,1] + lengthArrows*ErrorVF1[,1],
       derivative_estimator$x[,2] + lengthArrows*ErrorVF1[,2],
       length = 0.05, angle = 15, col = "red")

errorNorm = sqrt((derivative_estimator$estimator[,1] - JVF1x[,1])^2 + (derivative_estimator$estimator[,2] - JVF1x[,2])^2)/sqrt(JVF1x[,1]^2+JVF1x[,2]^2)
image.plot(x = unique(x[,1]), y = unique(x[,2]), z = matrix(log10(errorNorm), nrow=sqrt(nEval), ncol=sqrt(nEval)),xlab="x",ylab="y",main="error norm rel (log10)")


plot(x[,1], derivative_estimator$estimator[,1], type="p", col="red", pch=16,
     xlab="x", ylab="derivative", main="Derivative Estimator", ylim=c(-0.5,0.5))
lines(x[,1], -0.01*(12*x[,1]^2-2), col="black", lty=2)


derivative_estimator = compute_derivative_term(X0_raw, X0_star, x=x,
                                              kernel.type="gauss", D=NULL, 
                                              method.h="silverman", h=NULL, lambda=NULL, 
                                              sparse=FALSE, gc=FALSE, chunk_size=512, Y=Y2)

# plot the derivative estimator
plot(x[,2], derivative_estimator$estimator[,2], type="p", col="red", pch=16,
     xlab="x", ylab="derivative", main="Derivative Estimator", ylim=c(-0.1,0.1))
abline(h=-0.01*2, col="black", lty=2)

# # --- Testing within_transform ---

# # Generate clean test data
# nObs_test <- 50
# nT_test <- 5
# set.seed(123)

# # Fixed effects (sum to zero)
# FE_test <- matrix(rnorm(nObs_test * 2), nObs_test, 2)
# FE_test[,1] <- FE_test[,1] - mean(FE_test[,1])
# FE_test[,2] <- FE_test[,2] - mean(FE_test[,2])

# # Time effects (sum to zero)
# TE_test <- matrix(rnorm(nT_test * 2), nT_test, 2)
# TE_test[,1] <- TE_test[,1] - mean(TE_test[,1])
# TE_test[,2] <- TE_test[,2] - mean(TE_test[,2])

# # Random noise
# noise <- array(rnorm(nObs_test * 2 * nT_test, 0, 0.1), dim = c(nObs_test, 2, nT_test))

# # Construct test data: X_it = FE_i + TE_t + noise_it
# X_test <- array(0, dim = c(nObs_test, 2, nT_test))
# for (t in 1:nT_test) {
#     X_test[,,t] <- FE_test + matrix(rep(TE_test[t,], each=nObs_test), nrow=nObs_test)
# }
# X_test <- X_test + noise


# cat("\n--- Running Tests for within_transform ---\n")

# # Test 1: Fixed Effects (within-individual)
# cat("Test 1: FE=TRUE, TE=FALSE (within-individual transformation)\n")
# X_within_i <- within_transform(X_test, FE = TRUE, TE = FALSE, uniform_weights = TRUE, nEvals = 1)
# individual_means <- apply(X_within_i[,,,1], c(1, 2), mean)
# cat("Max absolute individual mean:", max(abs(individual_means)), "\n")


# # Test 2: Time Effects (within-time)
# cat("Test 2: FE=FALSE, TE=TRUE (within-time transformation)\n")
# X_within_t <- within_transform(X_test, FE = FALSE, TE = TRUE, uniform_weights = TRUE, nEvals = 1)
# time_means <- apply(X_within_t[,,,1], c(3, 2), mean)
# cat("Max absolute time mean:", max(abs(time_means)), "\n")

# # Test 3: Two-ways Effects
# cat("Test 3: FE=TRUE, TE=TRUE (two-ways transformation)\n")
# X_within_it <- within_transform(X_test, FE = TRUE, TE = TRUE, uniform_weights = TRUE, nEvals = 1)
# individual_means_2way <- apply(X_within_it[,,,1], c(1, 2), mean)
# time_means_2way <- apply(X_within_it[,,,1], c(3, 2), mean)
# cat("Max absolute individual mean (two-way):", max(abs(individual_means_2way)), "\n")
# cat("Max absolute time mean (two-way):", max(abs(time_means_2way)), "\n")

