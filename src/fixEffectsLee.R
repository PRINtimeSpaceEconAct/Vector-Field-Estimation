# Clear workspace and load dependencies
rm(list = ls())
DEBUG = FALSE
source("src/libs/loadLib.R")
library(dplyr)
library(plm)
library(scales)
library(fields)
library(sm)

# Parameters ----
nObs = 200  # Number of cross-sectional units
nT = 20     # Number of time periods for transitions
nEval = 2500  # Number of evaluation points for VF estimation
output_dir <- "test_pics_NonLinearSingleWell"

lengthArrows <- 0.1

# Rotation parameters
theta <- pi/3  # Rotation angle in radians (pi/3 = 60 degrees)
rotation_strength <- 1  # Overall scaling of the rotation


# Data generation ----
set.seed(1)  # For reproducibility

# Generate fixed effects (FE) - unit-specific heterogeneity  
#alpha_i = mvrnorm(nObs, mu=c(0,0), Sigma=0.01*diag(2))
#alpha_i[,1] = alpha_i[,1] - sum(alpha_i[,1])/nObs  # Center to sum to zero
#alpha_i[,2] = alpha_i[,2] - sum(alpha_i[,2])/nObs
alpha_i = matrix(0,nrow=nObs,ncol=2)


# Generate time effects (TE) - common time-varying factors
#gamma_t = mvrnorm(nT, mu=c(0,0), Sigma=0.01*diag(2))
#gamma_t[,1] = gamma_t[,1] - sum(gamma_t[,1])/nT  # Center to sum to zero
#gamma_t[,2] = gamma_t[,2] - sum(gamma_t[,2])/nT
gamma_t = matrix(0,nrow=nT,ncol=2)


# Initial conditions with FE and TE
X0 = mvrnorm(nObs, mu=c(0,0),Sigma = 0.5*diag(2)) + alpha_i + array(rep(gamma_t[1,],each=nObs),dim=c(nObs,2))
X = array(NA,dim = c(nObs,2,nT+1))
X[,,1] = X0

# example 1 - double well ----
VF <- function(X){
    # X = (x,y)
    # U(X) = y^4 - y^2 + x^2  (swapped from x^4 - x^2 + y^2)
    # VF(X) = -grad U(X) = -(2x, 4y^3 - 2y)
    return( -0.1*c(2*X[1], 4*X[2]^3 - 2*X[2]) )
}

# example 2 - non-linear single well ----
# VF <- function(X){
#     # X = (x,y)
#     # U(X) = x^4 + y^4
#     # VF(X) = -grad U(X) = -(4x^3, 4y^3)
#     return( -0.05*c(4*X[1]^3, 4*X[2]^3) )
# }

#example 3 - parameterized rotation around origin ----
# VF <- function(X){
#     # X = (x,y)
#     # Standard 2D rotation of vector (x,y) by angle theta
#     # VF(x,y) = rotation_strength * (x*cos(θ) - y*sin(θ), x*sin(θ) + y*cos(θ))
# 
#     x <- X[1]
#     y <- X[2]
# 
#     vf_x <- rotation_strength * (cos(theta) * x - sin(theta) * y)
#     vf_y <- rotation_strength * (sin(theta) * x + cos(theta) * y)
# 
#     return(c(vf_x, vf_y))
# }

VF_jacobian <- function(X){
    # X = (x,y)
    # VF(x,y) = -0.1 * (2x, 4y^3 - 2y)
    # J_11 = ∂(vf_x)/∂x = -0.2
    # J_12 = ∂(vf_x)/∂y = 0
    # J_21 = ∂(vf_y)/∂x = 0
    # J_22 = ∂(vf_y)/∂y = -0.1 * (12y^2 - 2)
    J_11 <- -0.2
    J_12 <- 0
    J_21 <- 0
    J_22 <- -0.1 * (12 * X[2]^2 - 2)
    return(matrix(c(J_11, J_12, J_21, J_22), nrow=2, byrow=TRUE))
}

# Generate panel data with FE, TE, and vector field dynamics ----
noise = array(rnorm(nObs*nT*2, sd = 0.01), dim=c(nObs,2,nT))
for (t in 1:nT){
    # X_{t+1} = X_t + VF(X_t) + FE + TE + noise
    vf_values <- t(apply(X[,,t, drop=FALSE], 1, VF))
    X[,,t+1] = X[,,t] + 
               vf_values +  # Vector field dynamics
               alpha_i +                   # Fixed effects  
               array(rep(gamma_t[t,], each=nObs), dim=c(nObs,2)) +  # Time effects
               noise[,,t]                  # Random noise
}

# Convert array data to long format data frame for panel analysis
# Note: We include all time periods (1 to nT+1) since the VF estimation 
# uses differences between consecutive periods
data_list <- list()
for (i in 1:nObs) {
  for (t in 1:(nT+1)) {
    data_list <- append(data_list, list(data.frame(
      id = paste0("unit_", i),
      time = t,
      X1 = X[i, 1, t],
      X2 = X[i, 2, t]
    )))
  }
}

# Combine all data
panel_data_df <- do.call(rbind, data_list)

# --- Data Filtering Step ---
# Identify individuals (units) that have non-finite values (NaN, Inf) in their history
ids_with_non_finite <- unique(panel_data_df$id[!is.finite(panel_data_df$X1) | !is.finite(panel_data_df$X2)])

if (length(ids_with_non_finite) > 0) {
  cat(sprintf("\n--- Data Filtering ---\n"))
  cat(sprintf("Found %d individuals with non-finite values. Removing their full histories from the dataset.\n", length(ids_with_non_finite)))
  
  # Filter out the problematic individuals
  panel_data_df <- panel_data_df[!(panel_data_df$id %in% ids_with_non_finite), ]
  
  cat(sprintf("Remaining individuals: %d\n\n", length(unique(panel_data_df$id))))
}


# Create panel data structure using plm  
panel_data <- pdata.frame(panel_data_df, index = c("id", "time"))



# Plot the true vector field for reference ----
cat("\n--- Generating True Vector Field Reference Plot ---\n")

# Create evaluation points (same method as in runPanelVFAnalysis)
X_unrolled <- aperm(X, c(1, 3, 2))
dim(X_unrolled) <- c(nObs * (nT+1), 2)  # X has nT+1 time periods
X_unrolled_finite <- X_unrolled[is.finite(X_unrolled[,1]) & is.finite(X_unrolled[,2]),]
x_eval_true <- defineEvalPoints(X_unrolled_finite, nEval)

# Compute true vector field on evaluation points
VF_true_eval <- t(apply(x_eval_true, 1, VF))

# Prepare data for plotting (using same settings as plotPanelVFAnalysis)
timeInterval <- 1
component_names <- c("X1", "X2")
rescale <- FALSE  # Following the same setting as in the loop

# Create the plot
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
svg(file.path(output_dir, "True_VectorField.svg"), width = 8, height = 8)
par(mar = c(4, 4, 3, 2))

# Plot the true vector field
plot(x_eval_true, type = "n",
     xlab = component_names[1],
     ylab = component_names[2],
     main = "True Vector Field - Double Well Potential (X/Y Swapped)")

if(rescale) {
  abline(h = 1, lty = 3)
  abline(v = 1, lty = 3)
} else {
  abline(h = 0, lty = 3)
  abline(v = 0, lty = 3)
}

# Draw vector field arrows
arrows(x_eval_true[, 1], x_eval_true[, 2],
       x_eval_true[, 1] + lengthArrows * VF_true_eval[, 1],
       x_eval_true[, 2] + lengthArrows * VF_true_eval[, 2],
       angle = 15, col = "blue", length = 0.05)

# Add observed data points
# Reshape X to have all points in a single 2D matrix for density calculation
points_for_density <- matrix(aperm(X, c(1, 3, 2)), ncol = dim(X)[2])
# points(points_for_density[, 1], points_for_density[, 2], pch = 19, col = "red", cex = 0.5)

# Add density contours if sm package is available
if (requireNamespace("sm", quietly = TRUE)) {
  points_for_density_complete <- points_for_density[complete.cases(points_for_density), ]
  if(nrow(points_for_density_complete) > 1) {
    tryCatch({
      # Let sm.density create its own evaluation grid, which is compatible with contour
      est.dens <- sm::sm.density(points_for_density_complete, display = "none")
      contour(est.dens$eval.points[, 1], est.dens$eval.points[, 2], est.dens$estimate,
              add = TRUE, col = "purple")
    }, error = function(e) {
      # Skip density contours if they fail
      cat("Error during density contour plotting:\n")
      print(e)
    })
  }
}

dev.off()

cat(sprintf("True vector field plot saved to: %s\n", file.path(output_dir, "True_VectorField.svg")))
cat("Blue arrows: True double well vector field (X/Y swapped)\n")
cat("Red points: Observed data points (all time steps)\n")
cat("Purple contours: Data density (all points)\n\n")

## Panel Vector Field Analysis (Modern Approach) ----
# Test different kernel types and parameters following VFLifeExpectancyFE.R structure
kernels_to_test <- c("gauss")
ridge_params_to_test <- c(1e-6)
rcond_thresholds_to_test <- c(1e-4)

total_iterations <- length(kernels_to_test) * length(ridge_params_to_test) * length(rcond_thresholds_to_test)
current_iteration <- 0

for (kernel in kernels_to_test) {
  for (ridge in ridge_params_to_test) {
    for (rcond_thresh in rcond_thresholds_to_test) {
      current_iteration <- current_iteration + 1
      
      kernel_type_str <- tools::toTitleCase(kernel)
      ridge_str <- if (ridge == 0) "0" else gsub("e-0", "e-", formatC(ridge, format = "e", digits = 0))
      rcond_str <- if (rcond_thresh == 0) "0" else gsub("e-0", "e-", formatC(rcond_thresh, format = "e", digits = 0))
      filename <- sprintf("Synthetic_%s_ridge_%s_cond%s.svg", kernel_type_str, ridge_str, rcond_str)
      
      cat(sprintf("\n--- [ %d / %d ] Processing: %s ---\n", current_iteration, total_iterations, filename))

      # Run panel VF analysis using the modern approach
      analysis_results <- runPanelVFAnalysis(panel_data,
                                             var_cols = c("X1", "X2"),
                                             id_col = "id",
                                             time_col = "time",
                                             nEval = nEval,
                                             FE = TRUE,
                                             TE = TRUE,
                                             uniform_weights = FALSE,
                                             kernel.type = kernel,
                                            #method.h = "silverman",
                                            # estimation_method = "LL",
                                            #  method.h = NULL,
                                              h=0.6,
                                             chunk_size = 1024,
                                             bootstrap_B = 1, 
                                             ridge_param = ridge,
                                             rcond_threshold = rcond_thresh,
                                             adaptive = TRUE)
      
      # Plot results using the modern approach
      plotPanelVFAnalysis(analysis_results,
                          timeInterval = 1,
                          rescale = FALSE,
                          rescale_ref_index = NULL,
                          years = 1:(nT+1),
                          label_names = paste0("unit_", 1:nObs),
                          component_names = c("X1", "X2"),
                          save_plots = TRUE,
                          save_path = output_dir,
                          lengthArrows = lengthArrows)
      
      # evaluation points
      x_eval <- analysis_results$x
      
      # --- Comparison of VF derivatives ---
      if (analysis_results$FE || analysis_results$TE) {
          
          VF_der_hat <- analysis_results$VF_derivatives_hat
          grad_F1_hat <- VF_der_hat$d1
          grad_F2_hat <- VF_der_hat$d2
          
          # Compute true derivatives
          VF_jac_true <- t(apply(x_eval, 1, VF_jacobian))
          grad_F1_true <- VF_jac_true[, c(1, 3)]
          grad_F2_true <- VF_jac_true[, c(2, 4)]
          
          # Arrow plot for grad(F1)
          svg(file.path(output_dir, paste0("VF_GradF1_Comparison_", filename)), width = 12, height = 6)
          par(mfrow=c(1,2), mar=c(4, 4, 3, 2))
          
          # Left plot: True grad(F1)
          plot(x_eval, type = "n", main = "True grad(F1)", xlab="X1", ylab="X2")
          arrows(x_eval[, 1], x_eval[, 2],
                 x_eval[, 1] + lengthArrows * grad_F1_true[, 1],
                 x_eval[, 2] + lengthArrows * grad_F1_true[, 2],
                 angle = 15, col = "blue", length = 0.05)
          
          # Right plot: Estimated grad(F1)
          plot(x_eval, type = "n", main = "Estimated grad(F1)", xlab="X1", ylab="X2")
          arrows(x_eval[, 1], x_eval[, 2],
                 x_eval[, 1] + lengthArrows * grad_F1_hat[, 1],
                 x_eval[, 2] + lengthArrows * grad_F1_hat[, 2],
                 angle = 15, col = "blue", length = 0.05)
          dev.off()
          
          # Arrow plot for grad(F2)
          svg(file.path(output_dir, paste0("VF_GradF2_Comparison_", filename)), width = 12, height = 6)
          par(mfrow=c(1,2), mar=c(4, 4, 3, 2))
          
          # Left plot: True grad(F2)
          plot(x_eval, type = "n", main = "True grad(F2)", xlab="X1", ylab="X2")
          arrows(x_eval[, 1], x_eval[, 2],
                 x_eval[, 1] + lengthArrows * grad_F2_true[, 1],
                 x_eval[, 2] + lengthArrows * grad_F2_true[, 2],
                 angle = 15, col = "blue", length = 0.05)
          
          # Right plot: Estimated grad(F2)
          plot(x_eval, type = "n", main = "Estimated grad(F2)", xlab="X1", ylab="X2")
          arrows(x_eval[, 1], x_eval[, 2],
                 x_eval[, 1] + lengthArrows * grad_F2_hat[, 1],
                 x_eval[, 2] + lengthArrows * grad_F2_hat[, 2],
                 angle = 15, col = "blue", length = 0.05)
          dev.off()
      }
      
      # --- Comparison Analysis: True vs Estimated ---
      
      # Ensure output directory exists
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      
      # 1) Vector field comparison on evaluation points
      VF_hat <- analysis_results$VF_hat
      
      # Compute true vector field on evaluation points
      VF_true <- t(apply(x_eval, 1, VF))
      
      # Calculate vector field differences
      VF_diff <- VF_hat - VF_true
      
      # Calculate norm (magnitude) of the vector difference
      VF_norm_diff <- sqrt(VF_diff[,1]^2 + VF_diff[,2]^2)
      
      # Calculate relative norm difference
      VF_norm_true <- sqrt(VF_true[,1]^2 + VF_true[,2]^2)
      VF_norm_rel_diff <- VF_norm_diff / VF_norm_true
      
      # Handle special values for absolute differences
      VF_norm_diff[is.nan(VF_norm_diff) | is.infinite(VF_norm_diff)] <- NA
      VF_norm_diff_log <- log10(VF_norm_diff)
      VF_norm_diff_log[is.infinite(VF_norm_diff_log)] <- NA # Replace -Inf with NA
      
      # Handle special values for relative differences, similar to testVF.R
      VF_norm_rel_diff[is.na(VF_norm_rel_diff) | !is.finite(VF_norm_rel_diff)] <- NA
      VF_norm_rel_diff_log <- log10(VF_norm_rel_diff)
      VF_norm_rel_diff_log[is.infinite(VF_norm_rel_diff_log)] <- NA # Replace -Inf with NA
      
      # Create heatmaps
      # The evaluation points form a regular grid
      grid_size <- round(sqrt(nEval))
      actual_nEval <- grid_size^2
      
      # Extract grid coordinates for proper heatmap visualization
      x_unique <- unique(x_eval[,1])
      y_unique <- unique(x_eval[,2])
      x_sorted <- sort(x_unique)
      y_sorted <- sort(y_unique)
      
      # Reshape norm differences into grid format
      VF_norm_diff_grid <- matrix(VF_norm_diff_log[1:actual_nEval], nrow = grid_size, ncol = grid_size, byrow = TRUE)
      VF_norm_rel_diff_grid <- matrix(VF_norm_rel_diff_log[1:actual_nEval], nrow = grid_size, ncol = grid_size, byrow = TRUE)
      
      # Create a single SVG with two plots
      svg(file.path(output_dir, paste0("VF_NormDiff_", filename)), width = 14, height = 6)
      
      # Use layout to give more space to the right plot
      layout(matrix(c(1, 2), nrow = 1), widths = c(1, 1.4))
      par(mar = c(4, 4, 3, 3))
      
      # 1) Absolute difference heatmap
      image.plot(x_sorted, y_sorted, VF_norm_diff_grid, 
                 main = "Absolute Difference (log10)",
                 xlab = "X1", ylab = "X2")
      
      # 2) Relative difference heatmap with percentage scale
      zlim_rel <- range(VF_norm_rel_diff_grid, na.rm = TRUE)
      rel_percent_vals <- 100 * VF_norm_rel_diff
      rel_percent_vals <- rel_percent_vals[is.finite(rel_percent_vals) & rel_percent_vals > 0]
      min_rel_p <- if (length(rel_percent_vals)) min(rel_percent_vals) else 0.1
      max_rel_p <- if (length(rel_percent_vals)) max(rel_percent_vals) else 100
      if (!is.finite(min_rel_p) || min_rel_p <= 0) min_rel_p <- 0.1
      if (!is.finite(max_rel_p) || max_rel_p <= min_rel_p) max_rel_p <- min_rel_p * 10
      candidate_perc <- c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)
      breaks_perc <- candidate_perc[candidate_perc >= min_rel_p & candidate_perc <= max_rel_p]
      if (length(breaks_perc) < 3) {
        breaks_perc <- unique(sort(c(min_rel_p, breaks_perc, max_rel_p)))
      }
      breaks_log <- log10(breaks_perc / 100)
      labels_rel <- paste0(formatC(breaks_perc, format = "fg", digits = 3), "%")
      
      # Increase bottom margin for horizontal legend on this panel only
      par(mar = c(6, 4, 3, 3))
      image.plot(x_sorted, y_sorted, VF_norm_rel_diff_grid, 
                 main = "Relative Difference (% of |V_true|)",
                 xlab = "X1", ylab = "X2",
                 zlim = zlim_rel,
                 horizontal = TRUE,
                 legend.shrink = 0.95,
                 legend.width = 1.6,
                 legend.mar = 4,
                 axis.args = list(at = breaks_log, labels = labels_rel, cex.axis = 0.9),
                 legend.args = list(text = "Relative error (%)", side = 1, line = 2.2, cex = 0.9))
       
       dev.off()
       
       # 2) Fixed and time effects comparison
      FE_hat <- analysis_results$FE_hat
      TE_hat <- analysis_results$TE_hat
      
      # Calculate differences for fixed and time effects
      FE_diff <- if (!is.null(FE_hat)) FE_hat - alpha_i else NULL
      TE_diff <- if (!is.null(TE_hat)) TE_hat - gamma_t else NULL
      
      # Calculate norms for FE and TE differences
      FE_norm_diff <- if (!is.null(FE_diff)) {
        norm_vals <- sqrt(FE_diff[,1]^2 + FE_diff[,2]^2)
        norm_vals[is.nan(norm_vals) | is.infinite(norm_vals)] <- NA
        norm_vals
      } else NULL
      
      TE_norm_diff <- if (!is.null(TE_diff)) {
        norm_vals <- sqrt(TE_diff[,1]^2 + TE_diff[,2]^2)
        norm_vals[is.nan(norm_vals) | is.infinite(norm_vals)] <- NA
        norm_vals
      } else NULL
      
      # Compute summary statistics (with NaN/NA handling)
      stats_summary <- list(
        VF_norm_differences = list(
          n_valid = sum(!is.na(VF_norm_diff)),
          n_total = length(VF_norm_diff),
          mean = mean(VF_norm_diff, na.rm = TRUE),
          sd = sd(VF_norm_diff, na.rm = TRUE),
          min = min(VF_norm_diff, na.rm = TRUE),
          max = max(VF_norm_diff, na.rm = TRUE),
          rmse = sqrt(mean(VF_norm_diff^2, na.rm = TRUE)),
          mae = mean(VF_norm_diff, na.rm = TRUE)  # MAE = mean for norm since it's already positive
        ),
        FE_norm_differences = if (!is.null(FE_norm_diff)) list(
          n_valid = sum(!is.na(FE_norm_diff)),
          n_total = length(FE_norm_diff),
          mean = mean(FE_norm_diff, na.rm = TRUE),
          sd = sd(FE_norm_diff, na.rm = TRUE),
          min = min(FE_norm_diff, na.rm = TRUE),
          max = max(FE_norm_diff, na.rm = TRUE),
          rmse = sqrt(mean(FE_norm_diff^2, na.rm = TRUE)),
          mae = mean(FE_norm_diff, na.rm = TRUE)  # MAE = mean for norm since it's already positive
        ) else NULL,
        TE_norm_differences = if (!is.null(TE_norm_diff)) list(
          n_valid = sum(!is.na(TE_norm_diff)),
          n_total = length(TE_norm_diff),
          mean = mean(TE_norm_diff, na.rm = TRUE),
          sd = sd(TE_norm_diff, na.rm = TRUE),
          min = min(TE_norm_diff, na.rm = TRUE),
          max = max(TE_norm_diff, na.rm = TRUE),
          rmse = sqrt(mean(TE_norm_diff^2, na.rm = TRUE)),
          mae = mean(TE_norm_diff, na.rm = TRUE)  # MAE = mean for norm since it's already positive
        ) else NULL,
        parameters = list(
          kernel = kernel,
          ridge_param = ridge,
          rcond_threshold = rcond_thresh,
          nObs = nObs,
          nT = nT,
          nEval = nEval
        )
      )
      
      # # Save summary statistics to file in a more readable format
      # stats_filename <- file.path(output_dir, paste0("Summary_Stats_", gsub("\\.svg$", ".txt", filename)))
      # 
      # # Create a nicely formatted text output
      # cat("=== VECTOR FIELD ESTIMATION COMPARISON STATISTICS ===\n\n", file = stats_filename)
      # cat("Parameters:\n", file = stats_filename, append = TRUE)
      # cat(sprintf("  Kernel: %s\n", kernel), file = stats_filename, append = TRUE)
      # cat(sprintf("  Ridge parameter: %g\n", ridge), file = stats_filename, append = TRUE)
      # cat(sprintf("  Condition threshold: %g\n", rcond_thresh), file = stats_filename, append = TRUE)
      # cat(sprintf("  Number of observations: %d\n", nObs), file = stats_filename, append = TRUE)
      # cat(sprintf("  Number of time periods: %d\n", nT), file = stats_filename, append = TRUE)
      # cat(sprintf("  Number of evaluation points: %d\n", nEval), file = stats_filename, append = TRUE)
      # cat("\n", file = stats_filename, append = TRUE)
      # 
      # # Vector Field Statistics
      # cat("VECTOR FIELD NORM DIFFERENCES ||Estimated - True||:\n", file = stats_filename, append = TRUE)
      # cat(sprintf("  Valid observations: %d / %d (%.1f%%)\n",
      #             stats_summary$VF_norm_differences$n_valid,
      #             stats_summary$VF_norm_differences$n_total,
      #             100 * stats_summary$VF_norm_differences$n_valid / stats_summary$VF_norm_differences$n_total), 
      #     file = stats_filename, append = TRUE)
      # cat(sprintf("  Mean: %12.8f | SD: %12.8f | RMSE: %12.8f | MAE: %12.8f\n",
      #             stats_summary$VF_norm_differences$mean,
      #             stats_summary$VF_norm_differences$sd,
      #             stats_summary$VF_norm_differences$rmse,
      #             stats_summary$VF_norm_differences$mae), file = stats_filename, append = TRUE)
      # cat(sprintf("  Min: %12.8f | Max: %12.8f\n",
      #             stats_summary$VF_norm_differences$min,
      #             stats_summary$VF_norm_differences$max), file = stats_filename, append = TRUE)
      # 
      # # Add warning if many invalid values
      # if (stats_summary$VF_norm_differences$n_valid < 0.9 * stats_summary$VF_norm_differences$n_total) {
      #   cat("  WARNING: Significant number of NaN/Inf values detected!\n", file = stats_filename, append = TRUE)
      # }
      # cat("\n", file = stats_filename, append = TRUE)
      # 
      # # Fixed Effects Statistics (if available)
      # if (!is.null(stats_summary$FE_norm_differences)) {
      #   cat("FIXED EFFECTS NORM DIFFERENCES ||Estimated - True||:\n", file = stats_filename, append = TRUE)
      #   cat(sprintf("  Valid observations: %d / %d (%.1f%%)\n",
      #               stats_summary$FE_norm_differences$n_valid,
      #               stats_summary$FE_norm_differences$n_total,
      #               100 * stats_summary$FE_norm_differences$n_valid / stats_summary$FE_norm_differences$n_total), 
      #       file = stats_filename, append = TRUE)
      #   cat(sprintf("  Mean: %12.8f | SD: %12.8f | RMSE: %12.8f | MAE: %12.8f\n",
      #               stats_summary$FE_norm_differences$mean,
      #               stats_summary$FE_norm_differences$sd,
      #               stats_summary$FE_norm_differences$rmse,
      #               stats_summary$FE_norm_differences$mae), file = stats_filename, append = TRUE)
      #   cat(sprintf("  Min: %12.8f | Max: %12.8f\n",
      #               stats_summary$FE_norm_differences$min,
      #               stats_summary$FE_norm_differences$max), file = stats_filename, append = TRUE)
      #   
      #   # Add warning if many invalid values
      #   if (stats_summary$FE_norm_differences$n_valid < 0.9 * stats_summary$FE_norm_differences$n_total) {
      #     cat("  WARNING: Significant number of NaN/Inf values detected!\n", file = stats_filename, append = TRUE)
      #   }
      #   cat("\n", file = stats_filename, append = TRUE)
      # } else {
      #   cat("FIXED EFFECTS: Not estimated\n\n", file = stats_filename, append = TRUE)
      # }
      # 
      # # Time Effects Statistics (if available)
      # if (!is.null(stats_summary$TE_norm_differences)) {
      #   cat("TIME EFFECTS NORM DIFFERENCES ||Estimated - True||:\n", file = stats_filename, append = TRUE)
      #   cat(sprintf("  Valid observations: %d / %d (%.1f%%)\n",
      #               stats_summary$TE_norm_differences$n_valid,
      #               stats_summary$TE_norm_differences$n_total,
      #               100 * stats_summary$TE_norm_differences$n_valid / stats_summary$TE_norm_differences$n_total), 
      #       file = stats_filename, append = TRUE)
      #   cat(sprintf("  Mean: %12.8f | SD: %12.8f | RMSE: %12.8f | MAE: %12.8f\n",
      #               stats_summary$TE_norm_differences$mean,
      #               stats_summary$TE_norm_differences$sd,
      #               stats_summary$TE_norm_differences$rmse,
      #               stats_summary$TE_norm_differences$mae), file = stats_filename, append = TRUE)
      #   cat(sprintf("  Min: %12.8f | Max: %12.8f\n",
      #               stats_summary$TE_norm_differences$min,
      #               stats_summary$TE_norm_differences$max), file = stats_filename, append = TRUE)
      #   
      #   # Add warning if many invalid values
      #   if (stats_summary$TE_norm_differences$n_valid < 0.9 * stats_summary$TE_norm_differences$n_total) {
      #     cat("  WARNING: Significant number of NaN/Inf values detected!\n", file = stats_filename, append = TRUE)
      #   }
      #   cat("\n", file = stats_filename, append = TRUE)
      # } else {
      #   cat("TIME EFFECTS: Not estimated\n\n", file = stats_filename, append = TRUE)
      # }
      # 
      # cat(sprintf("Generated at: %s\n", Sys.time()), file = stats_filename, append = TRUE)
      
      cat(sprintf("  - Norm heatmap saved to: VF_NormDiff_%s\n", filename))
      # cat(sprintf("  - Statistics saved to: %s\n", basename(stats_filename)))
      
      # Console warnings for NaN/Inf values
      if (stats_summary$VF_norm_differences$n_valid < 0.9 * stats_summary$VF_norm_differences$n_total) {
        cat(sprintf("  - WARNING: VF has %d/%d invalid values (NaN/Inf)\n", 
                    stats_summary$VF_norm_differences$n_total - stats_summary$VF_norm_differences$n_valid,
                    stats_summary$VF_norm_differences$n_total))
      }
      if (!is.null(stats_summary$FE_norm_differences) && 
          stats_summary$FE_norm_differences$n_valid < 0.9 * stats_summary$FE_norm_differences$n_total) {
        cat(sprintf("  - WARNING: FE has %d/%d invalid values (NaN/Inf)\n", 
                    stats_summary$FE_norm_differences$n_total - stats_summary$FE_norm_differences$n_valid,
                    stats_summary$FE_norm_differences$n_total))
      }
      if (!is.null(stats_summary$TE_norm_differences) && 
          stats_summary$TE_norm_differences$n_valid < 0.9 * stats_summary$TE_norm_differences$n_total) {
        cat(sprintf("  - WARNING: TE has %d/%d invalid values (NaN/Inf)\n", 
                    stats_summary$TE_norm_differences$n_total - stats_summary$TE_norm_differences$n_valid,
                    stats_summary$TE_norm_differences$n_total))
      }
                          
    }
  }
}





