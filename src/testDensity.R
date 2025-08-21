# Clear workspace and load dependencies
rm(list = ls())
DEBUG = TRUE
source("libs/loadLib.R")
source("libs/codeMunicipalities.R")
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(sm))
suppressPackageStartupMessages(library(mvtnorm))
# ---- Data Generation ----
# Generate random normal data for source distribution
nObs = 10000
set.seed(123)
# Sigma = matrix(c(1, 0, 0, 1), nrow=2)
Sigma = matrix(c(0.75, 0.5, 0.5, 0.75), nrow=2)

X_0_Gauss = mvrnorm(nObs, mu=c(0,0),Sigma = Sigma)



# Create evaluation grid for density estimation and regression
nEval = 2500
xGrid = seq(from=min(X_0_Gauss[,1]), to=max(X_0_Gauss[,1]), length.out=round(sqrt(nEval)))
yGrid = seq(from=min(X_0_Gauss[,2]), to=max(X_0_Gauss[,2]), length.out=round(sqrt(nEval)))
x = as.matrix(expand.grid(xGrid, yGrid))

mean = c(0, 0)
trueGaussian = dmvnorm(x, mean, Sigma)

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
estGaussian <- densityEst2d(X_0_Gauss, x=x, h=bandwidth, kernel.type = "gauss", gc=TRUE)
estGaussian.sm <- sm.density(X_0_Gauss, h=c(bandwidth, bandwidth), eval.points=x, eval.grid=FALSE, nbins=0)
print(paste("Maximum absolute difference between sm and custom:", max(abs(estGaussian$estimator - estGaussian.sm$estimate))))
# Test adaptive bandwidth density estimation
print("================================================")
print("Estimating a standard Gaussian with adaptive bandwidth")
alpha = 0.5
estAdaptive <- densityEst2dAdaptive(X_0_Gauss, x=x, kernel.type="gauss", gc=TRUE, chunk_size=1024, alpha=alpha)
print(paste("Maximum absolute difference between true and adaptive estimate:", max(abs(trueGaussian - estAdaptive$estimator))))

z_diff = trueGaussian - estAdaptive$estimator
print(plot_ly(x=x[,1], y=x[,2], z=z_diff, intensity=z_diff, type="mesh3d") %>%
          layout(title="Difference: True Standard Gaussian - Adaptive Bandwidth Estimate",
                 scene=list(
                     zaxis=list(title="Difference"),
                     xaxis=list(title="X"),
                     yaxis=list(title="Y")
                 )))

# Estimate the density with the municipalities code
estMunicipalities <- normKernelBivAdaptive(x, X_0_Gauss, alpha = alpha)
# Plot the difference between the adaptive and the municipalities code
z_diff = estMunicipalities - estAdaptive$estimator
print(plot_ly(x=x[,1], y=x[,2], z=z_diff, intensity=z_diff, type="mesh3d") %>%
          layout(title="Difference: Municipalities Implementation (Standard Gaussian) - Custom Adaptive",
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
# # Plot the true density
# print(plot_ly(x=x[,1], y=x[,2], z=trueBimodal, intensity=trueBimodal, type="mesh3d") %>%
#     layout(title="True Bimodal Density",
#            scene=list(
#                zaxis=list(title="Density"),
#                xaxis=list(title="X"),
#                yaxis=list(title="Y")
#            )))

# Scatter plot of the bimodal data
print(plot_ly(x=X_0_Bimodal[,1], y=X_0_Bimodal[,2], type="scatter", mode="markers", marker=list(size=2, opacity=0.6)) %>%
          layout(title="Bimodal Data Scatter Plot",
                 xaxis=list(title="X"),
                 yaxis=list(title="Y")))

print("Estimating the bimodal density")
estBimodal <- densityEst2d(X_0_Bimodal, x=x, kernel="gauss", gc=TRUE)
# # Plot the estimated density
# print(plot_ly(x=x[,1], y=x[,2], z=estBimodal$estimator, intensity=estBimodal$estimator, type="mesh3d") %>%
#     layout(title="Estimate: Mixture of Two Gaussians - Fixed Bandwidth Estimate",
#            scene=list(
#                zaxis=list(title="Difference"),
#                xaxis=list(title="X"),
#                yaxis=list(title="Y")
#            )))

# Plot the difference between the true and estimated density
z_diff = trueBimodal - estBimodal$estimator
print(plot_ly(x=x[,1], y=x[,2], z=z_diff, intensity=z_diff, type="mesh3d") %>%
          layout(title="Difference: Mixture of Two Gaussians - Fixed Bandwidth Estimate (Silverman)",
                 scene=list(
                     zaxis=list(title="Difference"),
                     xaxis=list(title="X"),
                     yaxis=list(title="Y")
                 )))
# Estimate the density with an adaptive bandwidth
print("Estimating the bimodal density with adaptive bandwidth")
alpha = 0.5
estBimodalAdaptive <- densityEst2dAdaptive(X_0_Bimodal, x=x, kernel.type="gauss", gc=TRUE, chunk_size=1024, alpha=alpha)
# # Plot the estimated density
# print(plot_ly(x=x[,1], y=x[,2], z=estBimodalAdaptive$estimator, intensity=estBimodalAdaptive$estimator, type="mesh3d") %>%
#     layout(title="Estimated Bimodal Density with Adaptive Bandwidth",
#            scene=list(
#                zaxis=list(title="Density"),
#                xaxis=list(title="X"),
#                yaxis=list(title="Y")
#            )))

# Plot the difference between the true and estimated density
z_diff = trueBimodal - estBimodalAdaptive$estimator
print(plot_ly(x=x[,1], y=x[,2], z=z_diff, intensity=z_diff, type="mesh3d") %>%
          layout(title="Difference: Mixture of Two Gaussians - Adaptive Bandwidth Estimate",
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
# # Plot the true density
# print(plot_ly(x=x[,1], y=x[,2], z=trueTrimodal, intensity=trueTrimodal, type="mesh3d") %>%
#     layout(title="True Trimodal Density",
#            scene=list(
#                zaxis=list(title="Density"),
#                xaxis=list(title="X"),
#                yaxis=list(title="Y")
#            )))

# Scatter plot of the trimodal data
print(plot_ly(x=X_0_Trimodal[,1], y=X_0_Trimodal[,2], type="scatter", mode="markers", marker=list(size=2, opacity=0.6)) %>%
          layout(title="Trimodal Data Scatter Plot",
                 xaxis=list(title="X"),
                 yaxis=list(title="Y")))

print("Estimating the trimodal density")
estTrimodal <- densityEst2d(X_0_Trimodal, x=x, kernel="gauss", gc=TRUE)
# Plot the estimated density
# print(plot_ly(x=x[,1], y=x[,2], z=estTrimodal$estimator, intensity=estTrimodal$estimator, type="mesh3d") %>%
#     layout(title="Difference: Mixture of Three Gaussians - Fixed Bandwidth Estimate",
#            scene=list(
#                zaxis=list(title="Difference"),
#                xaxis=list(title="X"),
#                yaxis=list(title="Y")
#            )))

# Plot the difference between the true and estimated density
z_diff = trueTrimodal - estTrimodal$estimator
print(plot_ly(x=x[,1], y=x[,2], z=z_diff, intensity=z_diff, type="mesh3d") %>%
          layout(title="Difference: Mixture of Three Gaussians - Fixed Bandwidth Estimate (Silverman)",
                 scene=list(
                     zaxis=list(title="Difference"),
                     xaxis=list(title="X"),
                     yaxis=list(title="Y")
                 )))
# Estimate the density with an adaptive bandwidth
print("Estimating the trimodal density with adaptive bandwidth")
alpha = 0.5
estTrimodalAdaptive <- densityEst2dAdaptive(X_0_Trimodal, x=x, kernel.type="gauss", gc=TRUE, chunk_size=1024, alpha=alpha)
# # Plot the estimated density
# print(plot_ly(x=x[,1], y=x[,2], z=estTrimodalAdaptive$estimator, intensity=estTrimodalAdaptive$estimator, type="mesh3d") %>%
#     layout(title="Estimated Trimodal Density with Adaptive Bandwidth",
#            scene=list(
#                zaxis=list(title="Density"),
#                xaxis=list(title="X"),
#                yaxis=list(title="Y")
#            )))

# Plot the difference between the true and estimated density
z_diff = trueTrimodal - estTrimodalAdaptive$estimator
print(plot_ly(x=x[,1], y=x[,2], z=z_diff, intensity=z_diff, type="mesh3d") %>%
          layout(title="Difference: Mixture of Three Gaussians - Adaptive Bandwidth Estimate",
                 scene=list(
                     zaxis=list(title="Difference"),
                     xaxis=list(title="X"),
                     yaxis=list(title="Y")
                 )))

