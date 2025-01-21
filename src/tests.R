# Clear workspace and load dependencies
rm(list = ls())
DEBUG = TRUE
source("src/libs/loadLib.R")
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

mean = c(0, 0)
sigma = matrix(c(1, 0, 0, 1), nrow=2)
trueGaussian = dmvnorm(x, mean, sigma)

# ---- Density Estimation Tests ----

print("#---- Density Estimation Tests ----")
print("We start by testing the fixed bandwidth density estimation against the sm package")
print("================================================")
print(paste("Number of observations:", nObs))
print(paste("Number of evaluation points:", nEval))
print("================================================")
print("Estimating a standard Gaussian")
bandwidth = 0.5
print(paste(" Gaussian kernel", "bandwidth =", bandwidth))
estGaussian <- densityEst2d(X_0_Gauss, x=x, h=bandwidth, kernel="gauss", sparse=FALSE, gc=TRUE)
estGaussian.sm <- sm.density(X_0_Gauss, h=c(bandwidth, bandwidth), eval.points=x, eval.grid=FALSE, nbins=0)
print(paste("Maximum absolute difference between sm and custom:", max(abs(estGaussian$estimator - estGaussian.sm$estimate))))
# Test adaptive bandwidth density estimation
print("================================================")
print("Estimating a standard Gaussian with adaptive bandwidth")
alpha = 0.5
estAdaptive <- densityEst2dAdaptive(X_0_Gauss, x=x, kernel.type="gauss", sparse=FALSE, gc=TRUE, chunk_size=1024, alpha=alpha)
print(paste("Maximum absolute difference between true and adaptive estimate:", max(abs(trueGaussian - estAdaptive$estimator))))

z_diff = trueGaussian - estAdaptive$estimator
print(plot_ly(x=x[,1], y=x[,2], z=z_diff, intensity=z_diff, type="mesh3d") %>%
    layout(title="Difference between true and adaptive estimate",
           scene=list(
               zaxis=list(title="Difference"),
               xaxis=list(title="X"),
               yaxis=list(title="Y")
           )))
print("================================================")
print("We now test our estimation against a bimodal distribution")
# Generate bimodal data from mixture of two normals
n1 = round(nObs/2)
n2 = nObs - n1
mean1 = c(-1, -1)
mean2 = c(1, 1)
sigma1_1 = 0.6
sigma1_2 = 0.6
rho1 = 0.6
sigma1 = matrix(c(sigma1_1^2, rho1*sigma1_1*sigma1_2, rho1*sigma1_1*sigma1_2, sigma1_2^2), nrow=2)
print(sigma1)
sigma2_1 = 0.4
sigma2_2 = 0.5
rho2 = -0.2
sigma2 = matrix(c(sigma2_1^2, rho2*sigma2_1*sigma2_2, rho2*sigma2_1*sigma2_2, sigma2_2^2), nrow=2)
print(sigma2)

X_0_Bimodal = rbind(
    rmvnorm(n1, mean1, sigma1),
    rmvnorm(n2, mean2, sigma2)
)

# True density for comparison
trueBimodal = 0.5 * dmvnorm(x, mean1, sigma1) + 0.5 * dmvnorm(x, mean2, sigma2)
# Plot the true density
print(plot_ly(x=x[,1], y=x[,2], z=trueBimodal, intensity=trueBimodal, type="mesh3d") %>%
    layout(title="True Bimodal Density",
           scene=list(
               zaxis=list(title="Density"),
               xaxis=list(title="X"),
               yaxis=list(title="Y")
           )))
           
# Scatter plot of the bimodal data
print(plot_ly(x=X_0_Bimodal[,1], y=X_0_Bimodal[,2], type="scatter", mode="markers", marker=list(size=2, opacity=0.6)) %>%
    layout(title="Bimodal Data Scatter Plot",
           xaxis=list(title="X"),
           yaxis=list(title="Y")))

print("Estimating the bimodal density")
estBimodal <- densityEst2d(X_0_Bimodal, x=x, kernel="gauss", sparse=FALSE, gc=TRUE)
# Plot the estimated density
print(plot_ly(x=x[,1], y=x[,2], z=estBimodal$estimator, intensity=estBimodal$estimator, type="mesh3d") %>%
    layout(title="Estimated Bimodal Density",
           scene=list(
               zaxis=list(title="Density"),
               xaxis=list(title="X"),
               yaxis=list(title="Y")
           )))

# Plot the difference between the true and estimated density
z_diff = trueBimodal - estBimodal$estimator
print(plot_ly(x=x[,1], y=x[,2], z=z_diff, intensity=z_diff, type="mesh3d") %>%
    layout(title="Difference between true and estimated density",
           scene=list(
               zaxis=list(title="Difference"),
               xaxis=list(title="X"),
               yaxis=list(title="Y")
           )))
# Estimate the density with an adaptive bandwidth
print("Estimating the bimodal density with adaptive bandwidth")
alpha = 0.5
estBimodalAdaptive <- densityEst2dAdaptive(X_0_Bimodal, x=x, kernel.type="gauss", sparse=FALSE, gc=TRUE, chunk_size=1024, alpha=alpha)
# Plot the estimated density
print(plot_ly(x=x[,1], y=x[,2], z=estBimodalAdaptive$estimator, intensity=estBimodalAdaptive$estimator, type="mesh3d") %>%
    layout(title="Estimated Bimodal Density with Adaptive Bandwidth",
           scene=list(
               zaxis=list(title="Density"),
               xaxis=list(title="X"),
               yaxis=list(title="Y")
           )))

# Plot the difference between the true and estimated density
z_diff = trueBimodal - estBimodalAdaptive$estimator
print(plot_ly(x=x[,1], y=x[,2], z=z_diff, intensity=z_diff, type="mesh3d") %>%
    layout(title="Difference between true and estimated density with adaptive bandwidth",
           scene=list(
               zaxis=list(title="Difference"),
               xaxis=list(title="X"),
               yaxis=list(title="Y")
           )))

print("================================================")
print("We now test our estimation against a trimodal distribution")
# Generate trimodal data from mixture of three normals
n1 = round(nObs/3)
n2 = round(nObs/3)
n3 = nObs - n1 - n2
mean1 = c(-1.5, -1.5)
mean2 = c(1.5, 1.5)
mean3 = c(0, 0)
sigma1_1 = 0.4
sigma1_2 = 0.4
rho1 = 0.6
sigma1 = matrix(c(sigma1_1^2, rho1*sigma1_1*sigma1_2, rho1*sigma1_1*sigma1_2, sigma1_2^2), nrow=2)
print(sigma1)
sigma2_1 = 0.3
sigma2_2 = 0.3
rho2 = -0.2
sigma2 = matrix(c(sigma2_1^2, rho2*sigma2_1*sigma2_2, rho2*sigma2_1*sigma2_2, sigma2_2^2), nrow=2)
print(sigma2)
sigma3_1 = 0.5
sigma3_2 = 0.5
rho3 = 0
sigma3 = matrix(c(sigma3_1^2, rho3*sigma3_1*sigma3_2, rho3*sigma3_1*sigma3_2, sigma3_2^2), nrow=2)
print(sigma3)

X_0_Trimodal = rbind(
    rmvnorm(n1, mean1, sigma1),
    rmvnorm(n2, mean2, sigma2),
    rmvnorm(n3, mean3, sigma3)
)

# True density for comparison
trueTrimodal = (1/3) * dmvnorm(x, mean1, sigma1) + (1/3) * dmvnorm(x, mean2, sigma2) + (1/3) * dmvnorm(x, mean3, sigma3)
# Plot the true density
print(plot_ly(x=x[,1], y=x[,2], z=trueTrimodal, intensity=trueTrimodal, type="mesh3d") %>%
    layout(title="True Trimodal Density",
           scene=list(
               zaxis=list(title="Density"),
               xaxis=list(title="X"),
               yaxis=list(title="Y")
           )))
           
# Scatter plot of the trimodal data
print(plot_ly(x=X_0_Trimodal[,1], y=X_0_Trimodal[,2], type="scatter", mode="markers", marker=list(size=2, opacity=0.6)) %>%
    layout(title="Trimodal Data Scatter Plot",
           xaxis=list(title="X"),
           yaxis=list(title="Y")))

print("Estimating the trimodal density")
estTrimodal <- densityEst2d(X_0_Trimodal, x=x, kernel="gauss", sparse=FALSE, gc=TRUE)
# Plot the estimated density
print(plot_ly(x=x[,1], y=x[,2], z=estTrimodal$estimator, intensity=estTrimodal$estimator, type="mesh3d") %>%
    layout(title="Estimated Trimodal Density",
           scene=list(
               zaxis=list(title="Density"),
               xaxis=list(title="X"),
               yaxis=list(title="Y")
           )))

# Plot the difference between the true and estimated density
z_diff = trueTrimodal - estTrimodal$estimator
print(plot_ly(x=x[,1], y=x[,2], z=z_diff, intensity=z_diff, type="mesh3d") %>%
    layout(title="Difference between true and estimated density",
           scene=list(
               zaxis=list(title="Difference"),
               xaxis=list(title="X"),
               yaxis=list(title="Y")
           )))
# Estimate the density with an adaptive bandwidth
print("Estimating the trimodal density with adaptive bandwidth")
alpha = 0.5
estTrimodalAdaptive <- densityEst2dAdaptive(X_0_Trimodal, x=x, kernel.type="gauss", sparse=FALSE, gc=TRUE, chunk_size=1024, alpha=alpha)
# Plot the estimated density
print(plot_ly(x=x[,1], y=x[,2], z=estTrimodalAdaptive$estimator, intensity=estTrimodalAdaptive$estimator, type="mesh3d") %>%
    layout(title="Estimated Trimodal Density with Adaptive Bandwidth",
           scene=list(
               zaxis=list(title="Density"),
               xaxis=list(title="X"),
               yaxis=list(title="Y")
           )))

# Plot the difference between the true and estimated density
z_diff = trueTrimodal - estTrimodalAdaptive$estimator
print(plot_ly(x=x[,1], y=x[,2], z=z_diff, intensity=z_diff, type="mesh3d") %>%
    layout(title="Difference between true and estimated density with adaptive bandwidth",
           scene=list(
               zaxis=list(title="Difference"),
               xaxis=list(title="X"),
               yaxis=list(title="Y")
           )))
# # Compare with sm package implementation
# library(sm)
# system.time(est.sm <- sm.density(X_0, eval.points=x, eval.grid=FALSE, nbins=0))
# print(paste("sm package bandwidth:", est.sm$h))

# # Visualize difference between implementations
# z = est.sm$estimate - est$estimator
# print(plot_ly(x=xCoord, y=yCoord, z=z, intensity=z, type="mesh3d") %>%
#     layout(title=paste("Density Estimation: sm package vs Custom Adaptive (alpha =", alpha, ")"),
#            scene=list(
#                zaxis=list(title="Difference"),
#                xaxis=list(title="X"),
#                yaxis=list(title="Y")
#            )))

# # ---- Bandwidth Comparison Tests ----
# # Compare different bandwidth approaches
# est_h025 <- densityEst2d(X_0, x=x, h=0.25, kernel="gauss", sparse=FALSE, gc=TRUE)
# est_h05_lambda <- densityEst2d(X_0, x=x, h=0.5, kernel="gauss", 
#                               lambda=rep(0.5, nrow(X_0)), sparse=FALSE, gc=TRUE)

# # Calculate and visualize differences
# z_diff = est_h025$densityEst - est_h05_lambda$densityEst
# print(plot_ly(x=xCoord, y=yCoord, z=z_diff, intensity=z_diff, type="mesh3d") %>%
#     layout(title="Density Estimation: h=0.25 vs h=0.5 with λ=0.5",
#            scene=list(
#                zaxis=list(title="Difference"),
#                xaxis=list(title="X"),
#                yaxis=list(title="Y")
#            )))

# # Print maximum difference
# max_abs_diff = max(abs(z_diff))
# print(paste("Maximum absolute difference between bandwidth approaches:", max_abs_diff))

# # ---- Regression Tests ----
# # Generate synthetic target data with known transformations
# X_1 = matrix(nrow=nObs, ncol=1)
# X_1[,1] = 2*X_0[,1] + 1      # First component: linear transformation y = 2x + 1 

# # Nadaraya-Watson regression for each component
# est_comp = NWregression(X_0, X_1[,1], x=x, h=0.5, kernel.type="gauss", sparse=FALSE, gc=TRUE)

# # Visualize first component regression
# print(plot_ly() %>%
#     add_trace(x=est_comp$x[,1], y=est_comp$estimator, name="Predicted", mode="markers", 
#               marker=list(size=2, opacity=0.6)) %>%
#     add_trace(x=xGrid, y=2*xGrid + 1, name="True", mode="lines", line=list(width=2)) %>%
#     layout(title="NW Regression: First Component (y = 2x + 1)",
#            xaxis=list(title="X"),
#            yaxis=list(title="Y")))

# # Compare with np package implementation using fixed bandwidth
# library(np)
# np_bw <- npregbw(xdat=X_0, ydat=X_1[,1], bws=c(0.5, 0.5), 
#                  bandwidth.compute=FALSE, 
#                  ckertype="gaussian",  # Specify Gaussian kernel
#                  ckerorder=2)  # Second-order kernel
# np_est <- npreg(bws=np_bw, 
#                 exdat=x)

# # Visualize comparison between custom NW and np package
# print(plot_ly() %>%
#     add_trace(x=x[,1], y=np_est$mean, name="np package", mode="markers",
#               marker=list(size=2, opacity=0.6)) %>%
#     add_trace(x=xGrid, y=2*xGrid + 1, name="True", mode="lines", line=list(width=2)) %>%
#     layout(title="Regression Comparison: Custom NW vs np package (h=0.5)",
#            xaxis=list(title="X"),
#            yaxis=list(title="Y")))

# # Make a 3D plot with the difference between the two estimators
# z_diff = est_comp$estimator - np_est$mean
# print(plot_ly(x=x[,1], y=x[,2], z=z_diff, intensity=z_diff, type="mesh3d") %>%
#     layout(title="Difference between Custom NW and np package (h=0.5)",
#            scene=list(
#                zaxis=list(title="Difference"),
#                xaxis=list(title="X"),
#                yaxis=list(title="Y")
#            )))

