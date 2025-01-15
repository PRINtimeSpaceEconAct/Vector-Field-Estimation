# Clear workspace and load dependencies
rm(list = ls())
DEBUG = TRUE
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
estAdaptive <- densityEst2dAdaptive(X_0, x=x, kernel.type="gauss", sparse=FALSE, gc=TRUE, chunk_size=1024, alpha=alpha)

# Visualize adaptive density estimate
xCoord = estAdaptive$x[,1]
yCoord = estAdaptive$x[,2]
z = estAdaptive$estimator
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
z = est.sm$estimate - est$estimator
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
X_1 = matrix(nrow=nObs, ncol=1)
X_1[,1] = 2*X_0[,1] + 1      # First component: linear transformation y = 2x + 1 

# Nadaraya-Watson regression for each component
est_comp = NWregression(X_0, X_1[,1], x=x, h=0.5, kernel.type="gauss", sparse=FALSE, gc=TRUE)

# Visualize first component regression
print(plot_ly() %>%
    add_trace(x=est_comp$x[,1], y=est_comp$estimator, name="Predicted", mode="markers", 
              marker=list(size=2, opacity=0.6)) %>%
    add_trace(x=xGrid, y=2*xGrid + 1, name="True", mode="lines", line=list(width=2)) %>%
    layout(title="NW Regression: First Component (y = 2x + 1)",
           xaxis=list(title="X"),
           yaxis=list(title="Y")))

# Compare with np package implementation using fixed bandwidth
library(np)
np_bw <- npregbw(xdat=X_0, ydat=X_1[,1], bws=c(0.5, 0.5), 
                 bandwidth.compute=FALSE, 
                 ckertype="gaussian",  # Specify Gaussian kernel
                 ckerorder=2)  # Second-order kernel
np_est <- npreg(bws=np_bw, 
                exdat=x)

# Visualize comparison between custom NW and np package
print(plot_ly() %>%
    add_trace(x=x[,1], y=np_est$mean, name="np package", mode="markers",
              marker=list(size=2, opacity=0.6)) %>%
    add_trace(x=xGrid, y=2*xGrid + 1, name="True", mode="lines", line=list(width=2)) %>%
    layout(title="Regression Comparison: Custom NW vs np package (h=0.5)",
           xaxis=list(title="X"),
           yaxis=list(title="Y")))

# Make a 3D plot with the difference between the two estimators
z_diff = est_comp$estimator - np_est$mean
print(plot_ly(x=x[,1], y=x[,2], z=z_diff, intensity=z_diff, type="mesh3d") %>%
    layout(title="Difference between Custom NW and np package (h=0.5)",
           scene=list(
               zaxis=list(title="Difference"),
               xaxis=list(title="X"),
               yaxis=list(title="Y")
           )))

