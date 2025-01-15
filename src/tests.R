# Clear workspace and load dependencies
rm(list = ls())
source("src/libs/loadLib.R")
suppressPackageStartupMessages(library(plotly))

# ---- Data Generation ----
# Generate random normal data for source distribution
nObs = 10000
set.seed(123)
X_0 = matrix(nrow=nObs, rnorm(2*nObs))

# Create evaluation grid for density estimation and regression
nEval = 2500
xGrid = seq(from=min(X_0[,1]), to=max(X_0[,1]), length.out=round(sqrt(nEval)))
yGrid = seq(from=min(X_0[,2]), to=max(X_0[,2]), length.out=round(sqrt(nEval)))
x = as.matrix(expand.grid(xGrid, yGrid))

# ---- Density Estimation Tests ----
# Test adaptive bandwidth density estimation
alpha = 0.5
est <- densityEst2dAdaptive(X_0, x=x, kernel="gauss", sparse=FALSE, gc=TRUE, alpha=alpha)

# Monitor memory usage
print(Sys.procmem())

# Visualize adaptive density estimate
xCoord = est$x[,1]
yCoord = est$x[,2]
z = est$densityEst
print(plot_ly(x=xCoord, y=yCoord, z=z, intensity=z, type="mesh3d") %>%
    layout(title="Adaptive Bandwidth Density Estimation (Gaussian Kernel)",
           scene=list(
               zaxis=list(title="Density"),
               xaxis=list(title="X"),
               yaxis=list(title="Y")
           )))

# Compare with sm package implementation
library(sm)
system.time(est.sm <- sm.density(X_0, eval.points=x, eval.grid=FALSE, nbins=0))
print(paste("sm package bandwidth:", est.sm$h))

# Visualize difference between implementations
z = est.sm$estimate - est$densityEst
print(plot_ly(x=xCoord, y=yCoord, z=z, intensity=z, type="mesh3d") %>%
    layout(title=paste("Density Estimation: sm package vs Custom Adaptive (alpha =", alpha, ")"),
           scene=list(
               zaxis=list(title="Difference"),
               xaxis=list(title="X"),
               yaxis=list(title="Y")
           )))

# ---- Bandwidth Comparison Tests ----
# Compare different bandwidth approaches
est_h025 <- densityEst2d(X_0, x=x, h=0.25, kernel="gauss", sparse=FALSE, gc=TRUE)
est_h05_lambda <- densityEst2d(X_0, x=x, h=0.5, kernel="gauss", 
                              lambda=rep(0.5, nrow(X_0)), sparse=FALSE, gc=TRUE)

# Calculate and visualize differences
z_diff = est_h025$densityEst - est_h05_lambda$densityEst
print(plot_ly(x=xCoord, y=yCoord, z=z_diff, intensity=z_diff, type="mesh3d") %>%
    layout(title="Density Estimation: h=0.25 vs h=0.5 with λ=0.5",
           scene=list(
               zaxis=list(title="Difference"),
               xaxis=list(title="X"),
               yaxis=list(title="Y")
           )))

# Print maximum difference
max_abs_diff = max(abs(z_diff))
print(paste("Maximum absolute difference between bandwidth approaches:", max_abs_diff))

# ---- Regression Tests ----
# Generate synthetic target data with known transformations
X_1 = matrix(nrow=nObs, ncol=2)
X_1[,1] = 2*X_0[,1] + 1      # First component: linear transformation y = 2x + 1 
X_1[,2] = -0.5*X_0[,2] + 3   # Second component: linear transformation y = -0.5x + 3

# Nadaraya-Watson regression for each component
est_comp1 = NWregression(X_0, X_1[,1], x=x, h=0.5, kernel="gauss", sparse=FALSE, gc=TRUE)
est_comp2 = NWregression(X_0, X_1[,2], x=x, h=0.5, kernel="gauss", sparse=FALSE, gc=TRUE)

# Visualize first component regression
print(plot_ly() %>%
    add_trace(x=est_comp1$x[,1], y=est_comp1$NWest, name="Predicted", mode="markers", 
              marker=list(size=2, opacity=0.6)) %>%
    add_trace(x=xGrid, y=2*xGrid + 1, name="True", mode="lines", line=list(width=2)) %>%
    layout(title="NW Regression: First Component (y = 2x + 1)",
           xaxis=list(title="X"),
           yaxis=list(title="Y")))

# Visualize second component regression
print(plot_ly() %>%
    add_trace(x=est_comp2$x[,2], y=est_comp2$NWest, name="Predicted", mode="markers",
              marker=list(size=2, opacity=0.6)) %>%
    add_trace(x=xGrid, y=-0.5*xGrid + 3, name="True", mode="lines", line=list(width=2)) %>%
    layout(title="NW Regression: Second Component (y = -0.5x + 3)",
           xaxis=list(title="X"),
           yaxis=list(title="Y")))

# Calculate and visualize regression error
true_comp1 = 2*x[,1] + 1
true_comp2 = -0.5*x[,2] + 3
error = sqrt((est_comp1$NWest - true_comp1)^2 + (est_comp2$NWest - true_comp2)^2)

print(plot_ly(x=xCoord, y=yCoord, z=error, intensity=error, type="mesh3d") %>%
    layout(title="Total Regression Error (Euclidean Distance)",
           scene=list(
               zaxis=list(title="Error"),
               xaxis=list(title="X"),
               yaxis=list(title="Y")
           )))