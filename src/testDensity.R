# Clear workspace and load dependencies
rm(list = ls())
DEBUG = FALSE
source("src/libs/loadLib.R")
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(mvtnorm))

# ---- Data Generation Parameters ----
nObs = 10000
nEval = 2500
set.seed(123)
alpha = 0.5

# ---- Unimodal (Gaussian) Data Generation ----
mean_unimodal = c(0,0)
sigma_unimodal = matrix(c(0.75, 0.5, 0.5, 0.75), nrow=2)
X_0_Gauss = mvrnorm(nObs, mu=mean_unimodal, Sigma = sigma_unimodal)
df_gauss <- as.data.frame(X_0_Gauss)
colnames(df_gauss) <- c("V1", "V2")

# ---- Bimodal Data Generation ----
n1_bimodal = round(nObs/2)
n2_bimodal = nObs - n1_bimodal
mean1_bimodal = c(-1, -1); mean2_bimodal = c(1, 1)
sigma1_bimodal = matrix(c(0.6^2, 0.6*0.6*0.6, 0.6*0.6*0.6, 0.6^2), nrow=2)
sigma2_bimodal = matrix(c(0.4^2, -0.2*0.4*0.5, -0.2*0.4*0.5, 0.5^2), nrow=2)
X_0_Bimodal = rbind(rmvnorm(n1_bimodal, mean1_bimodal, sigma1_bimodal), rmvnorm(n2_bimodal, mean2_bimodal, sigma2_bimodal))
df_bimodal <- as.data.frame(X_0_Bimodal)
colnames(df_bimodal) <- c("V1", "V2")

# ---- Trimodal Data Generation ----
n1_trimodal = round(nObs/3); n2_trimodal = round(nObs/3); n3_trimodal = nObs - n1_trimodal - n2_trimodal
mean1_trimodal = c(-1.5, -1.5); mean2_trimodal = c(1.5, 1.5); mean3_trimodal = c(0, 0)
sigma1_trimodal = matrix(c(0.4^2, 0.6*0.4*0.4, 0.6*0.4*0.4, 0.4^2), nrow=2)
sigma2_trimodal = matrix(c(0.3^2, -0.2*0.3*0.3, -0.2*0.3*0.3, 0.3^2), nrow=2)
sigma3_trimodal = matrix(c(0.5^2, 0, 0, 0.5^2), nrow=2)
X_0_Trimodal = rbind(rmvnorm(n1_trimodal, mean1_trimodal, sigma1_trimodal), rmvnorm(n2_trimodal, mean2_trimodal, sigma2_trimodal), rmvnorm(n3_trimodal, mean3_trimodal, sigma3_trimodal))
df_trimodal <- as.data.frame(X_0_Trimodal)
colnames(df_trimodal) <- c("V1", "V2")


# ---- Reusable Test Function ----
perform_density_test <- function(df, true_density_func, dist_name) {
    
    cat(paste("\n\n================================================\n"))
    cat(paste("--- Testing", dist_name, "Distribution ---\n"))
    cat(paste("================================================\n"))

    # --- Fixed Bandwidth Test ---
    cat("\n--- Fixed Bandwidth Analysis ---\n")
    analysis_fixed <- runDensityAnalysis(df = df, var_cols = colnames(df), nEval = nEval,
                                         adaptive = FALSE, show_plots = TRUE, 
                                         component_names = c("V1", "V2"))
    
    true_density_fixed <- true_density_func(analysis_fixed$x)
    z_diff_fixed <- true_density_fixed - analysis_fixed$density_results$estimator
    
    cat("\n--- Error Summary (Fixed Bandwidth) ---\n")
    cat(paste("Max Absolute Error:", max(abs(z_diff_fixed)), "\n"))
    cat(paste("Mean Absolute Error:", mean(abs(z_diff_fixed)), "\n"))
    
    print(plot_ly(x=analysis_fixed$x[,1], y=analysis_fixed$x[,2], z=z_diff_fixed, intensity=z_diff_fixed, type="mesh3d") %>%
              layout(title=paste("Difference:", dist_name, "(Fixed) vs True Density"),
                     scene=list(zaxis=list(title="Difference"), xaxis=list(title="V1"), yaxis=list(title="V2"))))

    # --- Adaptive Bandwidth Test ---
    cat("\n--- Adaptive Bandwidth Analysis ---\n")
    analysis_adaptive <- runDensityAnalysis(df = df, var_cols = colnames(df), nEval = nEval,
                                            adaptive = TRUE, alpha = alpha, show_plots = TRUE,
                                            component_names = c("V1", "V2"))

    true_density_adaptive <- true_density_func(analysis_adaptive$x)
    z_diff_adaptive <- true_density_adaptive - analysis_adaptive$density_results$estimator

    cat("\n--- Error Summary (Adaptive Bandwidth) ---\n")
    cat(paste("Max Absolute Error:", max(abs(z_diff_adaptive)), "\n"))
    cat(paste("Mean Absolute Error:", mean(abs(z_diff_adaptive)), "\n"))

    print(plot_ly(x=analysis_adaptive$x[,1], y=analysis_adaptive$x[,2], z=z_diff_adaptive, intensity=z_diff_adaptive, type="mesh3d") %>%
          layout(title=paste("Difference:", dist_name, "(Adaptive) vs True Density"),
                 scene=list(zaxis=list(title="Difference"), xaxis=list(title="V1"), yaxis=list(title="V2"))))
}

# ---- Execute Tests ----

# 1. Unimodal Test
true_unimodal <- function(x) dmvnorm(x, mean_unimodal, sigma_unimodal)
perform_density_test(df_gauss, true_unimodal, "Unimodal Gaussian")

# 2. Bimodal Test
true_bimodal <- function(x) 0.5 * dmvnorm(x, mean1_bimodal, sigma1_bimodal) + 0.5 * dmvnorm(x, mean2_bimodal, sigma2_bimodal)
perform_density_test(df_bimodal, true_bimodal, "Bimodal Mixture")

# 3. Trimodal Test
true_trimodal <- function(x) (1/3) * dmvnorm(x, mean1_trimodal, sigma1_trimodal) + (1/3) * dmvnorm(x, mean2_trimodal, sigma2_trimodal) + (1/3) * dmvnorm(x, mean3_trimodal, sigma3_trimodal)
perform_density_test(df_trimodal, true_trimodal, "Trimodal Mixture")


# ##############################################################################
# ---- 1D Density Tests ----
# ##############################################################################


cat("\n\n\n#################################################")
cat("\n--- STARTING 1D DENSITY ESTIMATION TESTS ---")
cat("\n#################################################\n")


# ---- 1D Unimodal (Gaussian) Data Generation ----
mean_unimodal_1d = 0
sd_unimodal_1d = 0.75
X_0_Gauss_1d = rnorm(nObs, mean = mean_unimodal_1d, sd = sd_unimodal_1d)
df_gauss_1d <- data.frame(V1 = X_0_Gauss_1d)

# ---- 1D Bimodal Data Generation ----
mean1_bimodal_1d = -1.5; sd1_bimodal_1d = 0.6
mean2_bimodal_1d = 1.5; sd2_bimodal_1d = 0.4
X_0_Bimodal_1d = c(rnorm(n1_bimodal, mean = mean1_bimodal_1d, sd = sd1_bimodal_1d), 
                   rnorm(n2_bimodal, mean = mean2_bimodal_1d, sd = sd2_bimodal_1d))
df_bimodal_1d <- data.frame(V1 = X_0_Bimodal_1d)

# ---- 1D Trimodal Data Generation ----
mean1_trimodal_1d = -2.0; sd1_trimodal_1d = 0.4
mean2_trimodal_1d = 2.0; sd2_trimodal_1d = 0.3
mean3_trimodal_1d = 0; sd3_trimodal_1d = 0.5
X_0_Trimodal_1d = c(rnorm(n1_trimodal, mean = mean1_trimodal_1d, sd = sd1_trimodal_1d), 
                    rnorm(n2_trimodal, mean = mean2_trimodal_1d, sd = sd2_trimodal_1d),
                    rnorm(n3_trimodal, mean = mean3_trimodal_1d, sd = sd3_trimodal_1d))
df_trimodal_1d <- data.frame(V1 = X_0_Trimodal_1d)

# ---- Reusable 1D Test and Plotting Function ----
plot_1d_comparison <- function(analysis_result, true_density_func, data_points, title) {
    x_grid <- analysis_result$x
    est_density <- analysis_result$density_results$estimator
    true_density <- true_density_func(x_grid)
    
    # Suppress dev.new() window for automated runs if needed, but keep for interactive
    dev.new()
    
    ylim <- range(c(est_density, true_density))
    
    plot(x_grid, est_density, type = 'l', col = 'blue', ylim = ylim,
         xlab = "Value", ylab = "Density", main = title)
    lines(x_grid, true_density, col = 'red', lty = 2)
    rug(data_points)
    legend("topright", legend = c("Estimated", "True"), col = c("blue", "red"), lty = 1:2, bg="white")
    grid()
}

perform_density_test_1d <- function(df, true_density_func, dist_name) {
    
    cat(paste("\n\n================================================\n"))
    cat(paste("--- Testing 1D", dist_name, "Distribution ---\n"))
    cat(paste("================================================\n"))

    # --- Fixed Bandwidth Test ---
    cat("\n--- Fixed Bandwidth Analysis ---\n")
    analysis_fixed <- runDensityAnalysis(df = df, var_cols = colnames(df), nEval = nEval,
                                         adaptive = FALSE, show_plots = FALSE)
    
    true_density_fixed <- true_density_func(analysis_fixed$x)
    z_diff_fixed <- true_density_fixed - analysis_fixed$density_results$estimator
    
    cat("\n--- Error Summary (Fixed Bandwidth) ---\n")
    cat(paste("Max Absolute Error:", max(abs(z_diff_fixed)), "\n"))
    cat(paste("Mean Absolute Error:", mean(abs(z_diff_fixed)), "\n"))
    plot_1d_comparison(analysis_fixed, true_density_func, df[[1]], paste("Comparison:", dist_name, "(Fixed)"))

    # --- Adaptive Bandwidth Test ---
    cat("\n--- Adaptive Bandwidth Analysis ---\n")
    analysis_adaptive <- runDensityAnalysis(df = df, var_cols = colnames(df), nEval = nEval,
                                            adaptive = TRUE, alpha = alpha, show_plots = FALSE)

    true_density_adaptive <- true_density_func(analysis_adaptive$x)
    z_diff_adaptive <- true_density_adaptive - analysis_adaptive$density_results$estimator

    cat("\n--- Error Summary (Adaptive Bandwidth) ---\n")
    cat(paste("Max Absolute Error:", max(abs(z_diff_adaptive)), "\n"))
    cat(paste("Mean Absolute Error:", mean(abs(z_diff_adaptive)), "\n"))
    plot_1d_comparison(analysis_adaptive, true_density_func, df[[1]], paste("Comparison:", dist_name, "(Adaptive)"))
}

# ---- Execute 1D Tests ----

# 1. 1D Unimodal Test
true_unimodal_1d <- function(x) dnorm(x, mean_unimodal_1d, sd_unimodal_1d)
perform_density_test_1d(df_gauss_1d, true_unimodal_1d, "Unimodal Gaussian")

# 2. 1D Bimodal Test
true_bimodal_1d <- function(x) 0.5 * dnorm(x, mean1_bimodal_1d, sd1_bimodal_1d) + 0.5 * dnorm(x, mean2_bimodal_1d, sd2_bimodal_1d)
perform_density_test_1d(df_bimodal_1d, true_bimodal_1d, "Bimodal Mixture")

# 3. 1D Trimodal Test
true_trimodal_1d <- function(x) (1/3) * dnorm(x, mean1_trimodal_1d, sd1_trimodal_1d) + 
                                (1/3) * dnorm(x, mean2_trimodal_1d, sd2_trimodal_1d) + 
                                (1/3) * dnorm(x, mean3_trimodal_1d, sd3_trimodal_1d)
perform_density_test_1d(df_trimodal_1d, true_trimodal_1d, "Trimodal Mixture")

