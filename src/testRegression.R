# Clear workspace and load dependencies
rm(list = ls())
DEBUG = TRUE
source("src/libs/loadLib.R")
source("src/libs/codeMunicipalities.R")
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(sm))
suppressPackageStartupMessages(library(mvtnorm))

# ---- Data Generation ----
# Generate random normal data for source distribution
nObs = 10000
set.seed(123)
X_0_Gauss = matrix(nrow=nObs, rnorm(2*nObs))

# Create evaluation grid for density estimation and regression
nEval = 2500
xGrid = seq(from=min(X_0_Gauss[,1]), to=max(X_0_Gauss[,1]), length.out=round(sqrt(nEval)))
yGrid = seq(from=min(X_0_Gauss[,2]), to=max(X_0_Gauss[,2]), length.out=round(sqrt(nEval)))
x = as.matrix(expand.grid(xGrid, yGrid))


# ---- Regression Tests ----
# Generate synthetic target data with known transformations
X_1 = matrix(nrow=nObs, ncol=2)
X_1[,1] = 0.5*X_0_Gauss[,1]^2     # First component: linear transformation y = 2x + 1 
X_1[,2] = 0.5*X_0_Gauss[,2]     # Second component: linear transformation y = 3x + 2
target = 0.5*x[,1]^2

# Nadaraya-Watson regression for each component
est_comp = NWregression(X_0_Gauss, X_1[,1], x=x, chunk_size=1024, kernel.type="gauss", sparse=FALSE, gc=TRUE)
est_comp_adaptive = NWregressionAdaptive(X_0_Gauss, X_1[,1], x=x, chunk_size=1024, kernel.type="gauss", sparse=FALSE, gc=TRUE)
est_comp_LL = LLregression(X_0_Gauss, X_1[,1], x=x, chunk_size=1024, kernel.type="gauss", sparse=FALSE, gc=TRUE)
est_comp_LL_adaptive = LLregressionAdaptive(X_0_Gauss, X_1[,1], x=x, chunk_size=1024, kernel.type="gauss", sparse=FALSE, gc=TRUE)

# Do a scatter plot of the forecasted values
plot_ly(x=est_comp$x[,1], y=est_comp$estimator, type="scatter", mode="markers") %>%
    add_lines(x=x[,1], y=target, type="scatter", mode="lines")

plot_ly(x=est_comp_adaptive$x[,1], y=est_comp_adaptive$estimator, type="scatter", mode="markers") %>%
    add_lines(x=x[,1], y=target, type="scatter", mode="lines")

# Do a scatter plot of the forecasted values
plot_ly(x=est_comp_LL$x[,1], y=est_comp_LL$estimator, type="scatter", mode="markers") %>%
    add_lines(x=x[,1], y=target, type="scatter", mode="lines")

plot_ly(x=est_comp_LL_adaptive$x[,1], y=est_comp_LL_adaptive$estimator, type="scatter", mode="markers") %>%
    add_lines(x=x[,1], y=target, type="scatter", mode="lines")

# Compare with np package implementation using fixed bandwidth
suppressPackageStartupMessages(library(np))
np_bw <- npregbw(xdat=X_0_Gauss, ydat=X_1[,1], bws=c(0.5, 0.5), 
                 bandwidth.compute=FALSE, 
                 ckertype="gaussian",  # Specify Gaussian kernel
                 ckerorder=2)  # Second-order kernel
np_est <- npreg(bws=np_bw, 
                exdat=x)

# Make a 3D plot with the difference between the two estimators
z_diff = np_est$mean - est_comp$estimator
print(plot_ly(x=x[,1], y=x[,2], z=z_diff, intensity=z_diff, type="mesh3d") %>%
    layout(title="Difference: np Package Implementation - Custom NW ",
           scene=list(
               zaxis=list(title="Difference"),
               xaxis=list(title="X"),
               yaxis=list(title="Y")
           )))

# Make a 3D plot with the difference between the true target and the estimated target
z_diff = target - np_est$mean
print(plot_ly(x=x[,1], y=x[,2], z=z_diff, intensity=z_diff, type="mesh3d") %>%
    layout(title="Difference: True Linear Transform - np Package Estimate",
           scene=list(
               zaxis=list(title="Difference"),
               xaxis=list(title="X"),
               yaxis=list(title="Y")
           )))

est_comp_adaptive = NWregressionAdaptive(X_0_Gauss, X_1[,1], x=x, kernel.type="gauss", chunk_size=1024, sparse=FALSE, gc=TRUE, alpha=0.5)

z_diff = target - est_comp_adaptive$estimator
print(plot_ly(x=x[,1], y=x[,2], z=z_diff, intensity=z_diff, type="mesh3d") %>%
    layout(title="Difference: True Linear Transform - Adaptive NW Estimate",
           scene=list(
               zaxis=list(title="Difference"),
               xaxis=list(title="X"),
               yaxis=list(title="Y")
           )))

