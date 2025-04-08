# Clear workspace and load dependencies
#setwd("~/Library/CloudStorage/OneDrive-UniversityofPisa/timeSpaceEvolutionEcAct/RVF/R code/Vector Field Estimation/")
rm(list = ls())
DEBUG = TRUE
source("src/libs/loadLib.R")
library(fields)
library(latex2exp)


# parameters ----
nObs = 1000
nEval = 2500

# data Generation ----
set.seed(1)
X0 = mvrnorm(nObs, mu=c(0,0),Sigma = 0.05*diag(2))
# compute min and max distance between points
# get distances component wise
#distances = computeDcomponents(X0, X0)
#distance_matrix = sqrt(mahalanobis(distances$z1, distances$z2, A = diag(2), 1))
# print the non zero minimum and maximum distance
#print(min(distance_matrix[distance_matrix > 0]))
#print(max(distance_matrix[distance_matrix > 0]))

# example 1 - double well ----
# VF <- function(X){
#     # X = (x,y)
#     # U(X) = x^4 - x^2 + y^2
#     # VF(X) = -grad U(X) = -(4x^3 - 2x, 2y)
#     return( -0.1*c(4*X[1]^3 - 2*X[1], 2*X[2]) )
# }

# example 2 -- single well ----
# VF <- function(X){
#     # X = (x,y)
#     # U(X) = x^2 + y^2
#     # VF(X) = -grad U(X) = -(2x, 2y)
#     return( -0.01*c(2*X[1], 2*X[2]) )
# }

# example 3 -- rotation ----
M = matrix(c(cos(pi/2), -sin(pi/2), sin(pi/2), cos(pi/2)),nrow=2,ncol=2)
VF <- function(X){
    # X = (x,y), theta = pi/4
    return (M %*% X)
}

# apply VF
X1 = X0 + t(apply(X0, 1, VF)) +  matrix(rnorm(2*nObs),nrow=nObs) %*%  matrix(c(0.01,0.005,0.005,0.02),nrow=2) * 10

# eval points
xGrid = seq(from=-1, to=1, length.out=round(sqrt(nEval)))
yGrid = seq(from=-1, to=1, length.out=round(sqrt(nEval)))
x = as.matrix(expand.grid(xGrid, yGrid))

VFx = t(apply(x, 1, VF))
# stima ----
# est_field_LL = LLfield(X0, X1, x=x, kernel.type="epa",method.h = "sj",
#                         chunk_size=1000,
#                         sparse=FALSE, gc=TRUE)
est_field_NW_opt = NWfield(X0, X1, x=x, kernel.type="gauss",method.h = "sj",
                       chunk_size=1000,
                       sparse=FALSE, gc=TRUE, hOpt = TRUE)

est_field_NW_adaptive = NWfieldAdaptive(X0, X1, x=x, kernel.type="gauss",method.h = "silverman",
                       chunk_size=1000,
                       sparse=FALSE, gc=TRUE, hOpt = TRUE, h = NULL, alpha=NULL, alphaOpt = TRUE, nGridAlpha=10)

contour(x=est_field_NW_adaptive$hGrid,  y=est_field_NW_adaptive$alphaGrid, z=est_field_NW_adaptive$AICc)
#indMin = which(est_field_NW_adaptive$AICc==min(est_field_NW_adaptive$AICc),arr.ind = TRUE)
points(est_field_NW_adaptive$h,est_field_NW_adaptive$alpha, pch=19)

kernel.type = "gauss"
list.h = define_h_method.h(X0, NULL ,"silverman", kernel.type)
hStart = list.h$h/10
hEnd = list.h$h*2
nGridh = 10
hGrid = exp(log(hStart) + (log(hEnd) - log(hStart)) * (0:(nGridh-1))/(nGridh-1))

est_dim = cbind(2500,2)

# Initialize array with correct dimensions
estimators = array(NA, dim = c(nGridh, est_dim[1], est_dim[2]))

for (i in 1:nGridh){
    est_field_NW = NWfield(X0, X1, x=x, kernel.type="gauss",
                       chunk_size=1000,
                       sparse=FALSE, gc=TRUE, hOpt = FALSE, h = hGrid[i])
    estimators[i,,] = est_field_NW$estimator
}
# Compute the mean squared error between the estimated and true VF
mse = array(NA, dim = nGridh)
error_cricci = array(NA, dim = nGridh)
for (i in 1:nGridh){
    mse[i] = sum(rowSums((estimators[i,,] - VFx)^2,na.rm=T)) / sum(is.finite(rowSums((estimators[i,,] - VFx))))
    # Substitute the NAs with 0
    estimators[i,,][is.na(estimators[i,,])] = 0
    error_cricci[i] = sum(rowSums((estimators[i,,] - VFx)^2)/rowSums(VFx^2)) 
}
print(paste("MSE: ",paste(mse,collapse=" ")))
print(paste("error_cricci: ",paste(error_cricci,collapse=" ")))
# Print all the optimal bandwidths
print(paste("Optimal h for AICc", hGrid[which.min(est_field_NW_opt$AICc)], "at index", which.min(est_field_NW_opt$AICc)))
print(paste("Optimal h for mse", hGrid[which.min(mse)], "at index", which.min(mse)))
print(paste("Optimal h for error_cricci", hGrid[which.min(error_cricci)], "at index", which.min(error_cricci)))
# Normalize the errors
mse = mse / max(mse)
error_cricci = error_cricci / max(error_cricci)
# Handle AICc normalization considering it could have both positive and negative values
AICc_values = est_field_NW_opt$AICc
if(min(AICc_values) < 0) {
  # If there are negative values, shift to make all values positive before normalizing
  AICc_shifted = AICc_values - min(AICc_values) + 1e-10  # Add small constant to avoid zeros
  AICc = AICc_shifted / max(AICc_shifted)
} else {
  # If all values are positive, normalize directly
  AICc = AICc_values / max(AICc_values)
}
# Plot the errors
plot(hGrid, mse, type = "l", col = "red", lwd = 2, xlab = "Bandwidth", ylab = "Normalized Error", main = "Normalized Errors",ylim=range(c(mse,error_cricci,AICc)))
lines(hGrid, error_cricci, col = "blue", lwd = 2)
lines(hGrid, AICc, col = "green", lwd = 2)
legend("topright", legend = c("MSE", "Error Cricci", "AICc"), col = c("red", "blue", "green"), lwd = 2)

# plot ----
## plot campo vero ----
# dev.new()
# op <- par(family = "mono") #Possible families: "mono", "Helvetica","Palatino" or "Times" 
# VFx = t(apply(x, 1, VF))
# plot(x, type = "n", xlab = TeX(r'($X_1$)'), ylab=TeX(r'($X_2$)'), main = "")
# arrows(x[,1],x[,2],x[,1]+VFx[,1],x[,2]+VFx[,2],angle=15,col="blue",length=0.05)
# abline(h=0)
# abline(v=0)
# dev.copy2pdf(file="testPics/doubleWellcampoVero.pdf")
# par(op)

## plot campo stimato ----
# est_field = est_field_NW_SJ
# lengthArrows=0.1

# # dev.new()
# plot(x, type = "n", xlab = TeX(r'($X_1$)'), ylab=TeX(r'($X_2$)'), main = "")
# arrows(est_field$x[,1], est_field$x[,2],
#        est_field$x[,1] + lengthArrows*est_field$estimator[,1], 
#        est_field$x[,2] + lengthArrows*est_field$estimator[,2],
#        length = 0.05, angle = 15, col = "blue")
# abline(h=0)
# abline(v=0)
# points(X0)


# est_field_NW_SJ$estimator[is.na(est_field_NW_SJ$estimator)] = 0
# est_field_NW$estimator[is.na(est_field_NW$estimator)] = 0

#sum(rowSums((est_field_NW_SJ$estimator - VFx)^2,na.rm=T)) / sum(is.finite(rowSums((est_field_NW_SJ$estimator - VFx))))
#sum(rowSums((est_field_NW$estimator - VFx)^2,na.rm=T)) / sum(is.finite(rowSums((est_field_NW$estimator - VFx))))

   
