# Clear workspace
rm(list = ls())
# Load custom libraries and suppress plotly startup messages
source("src/libs/loadLib.R")
suppressPackageStartupMessages(library(plotly))

# Generate random normal data
nObs = 10000
set.seed(123)
X_0 = matrix(nrow=nObs,rnorm(2*nObs))

# Create evaluation grid
nEval = 2500
xGrid = seq(from=min(X_0[,1]),to=max(X_0[,1]),length.out=round(sqrt(nEval)))
yGrid = seq(from=min(X_0[,2]),to=max(X_0[,2]),length.out=round(sqrt(nEval)))
x = as.matrix(expand.grid(xGrid,yGrid))

# Adaptive bandwidth parameter
alpha = 0.5

# Test the adaptive bandwidth density estimation using Gaussian kernel
est <- densityEst2dAdaptive(X_0,x=x, kernel = "gauss", sparse=FALSE,gc=TRUE, alpha = alpha)

# Print memory usage
print(Sys.procmem())

# Create a 3D surface plot of the density estimate
xCoord = est$x[,1]
yCoord = est$x[,2]
z = est$densityEst
print(plot_ly(x = xCoord, y = yCoord, z = z, intensity = z, type = "mesh3d") %>%
    layout(title = "Adaptive Bandwidth Density Estimation"))

# Compare with sm package implementation
library(sm)
# Time the sm.density calculation
system.time(est.sm <- sm.density(X_0,eval.points = x,eval.grid=FALSE,nbins=0))
print(paste("sm bandwidth:",est.sm$h))
# Plot difference between sm and our implementation
z = est.sm$estimate - est$densityEst
print(plot_ly(x = xCoord, y = yCoord, z = z, intensity = z, type = "mesh3d") %>%
    layout(title = paste("Difference between sm and Custom Adaptive Implementation with alpha =", alpha)))

# Create synthetic target data with affine transformations plus noise
X_1 = matrix(nrow=nObs, ncol=2)
X_1[,1] = 2*X_0[,1] + 1  # First component: 2x + 1 
X_1[,2] = -0.5*X_0[,2] + 3  # Second component: -0.5x + 3


# Predict each component using NW regression
est_comp1 = NWregression(X_0, X_1[,1], x=x, h=0.5, kernel="gauss", sparse=FALSE, gc=TRUE)
est_comp2 = NWregression(X_0, X_1[,2], x=x, h=0.5, kernel="gauss", sparse=FALSE, gc=TRUE)

#est_comp1 = NWregressionAdaptive(X_0, X_1[,1], x=x, alpha=0.5, kernel="gauss", sparse=FALSE, gc=TRUE)
#est_comp2 = NWregressionAdaptive(X_0, X_1[,2], x=x, alpha=0.5, kernel="gauss", sparse=FALSE, gc=TRUE)

print(plot_ly() %>%
    add_trace(x=est_comp1$x[,1], y=est_comp1$NWest, name="Predicted", mode="markers") %>%
    add_trace(x=xGrid, y=2*xGrid + 1, name="True", mode="lines") %>%
    layout(title = "NW Regression: First Component (2x + 1)",
           xaxis = list(title = "X"),
           yaxis = list(title = "Y")))

print(plot_ly() %>%
    add_trace(x=est_comp2$x[,2], y=est_comp2$NWest, name="Predicted", mode="markers") %>%
    add_trace(x=xGrid, y=-0.5*xGrid + 3, name="True", mode="lines") %>%
    layout(title = "NW Regression: Second Component (-0.5x + 3)",
           xaxis = list(title = "X"),
           yaxis = list(title = "Y")))

# Calculate true values on the evaluation grid
true_comp1 = 2*x[,1] + 1
true_comp2 = -0.5*x[,2] + 3

# Calculate total squared error at each point
error = sqrt((est_comp1$NWest - true_comp1)^2 + (est_comp2$NWest - true_comp2)^2)

# Create 3D surface plot of the error
print(plot_ly(x = xCoord, y = yCoord, z = error, intensity = error, type = "mesh3d") %>%
    layout(title = "Total Regression Error (Euclidean Distance)",
           scene = list(
               zaxis = list(title = "Error"),
               xaxis = list(title = "X"),
               yaxis = list(title = "Y")
           )))