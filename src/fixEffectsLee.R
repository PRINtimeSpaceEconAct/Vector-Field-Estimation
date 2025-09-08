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
nObs = 500  # Number of cross-sectional units
nT = 20     # Number of time periods for transitions
nEval = 1024  # Number of evaluation points for VF estimation
output_dir <- "test_pics_NonLinearSingleWell"

# Data generation ----
set.seed(3)  # For reproducibility

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
X0 = mvrnorm(nObs, mu=c(0,0),Sigma = 0.25*diag(2)) + alpha_i + array(rep(gamma_t[1,],each=nObs),dim=c(nObs,2))
X = array(NA,dim = c(nObs,2,nT+1))
X[,,1] = X0

# example 1 - double well ----
# VF <- function(X){
#     # X = (x,y)
#     # U(X) = x^4 - x^2 + y^2
#     # VF(X) = -grad U(X) = -(4x^3 - 2x, 2y)
#     return( -0.05*c(4*X[1]^3 - 4*X[1], 3*X[2]) )
# }

# example 2 - non-linear single well ----
VF <- function(X){
    # X = (x,y)
    # U(X) = x^4 + y^4
    # VF(X) = -grad U(X) = -(4x^3, 4y^3)
    return( -0.05*c(4*X[1]^3, 4*X[2]^3) )
}

# example 3 - simple counterclockwise rotation around origin ----
# VF <- function(X){
#     # X = (x,y)
#     # Simple counterclockwise rotation around origin
#     # VF(x,y) = (-y, x) gives counterclockwise rotation
#     
#     return( 0.1 * c(-X[2], X[1]) )
# }

# Generate panel data with FE, TE, and vector field dynamics ----
noise = array(rnorm(nObs*nT*2, sd = 0.01), dim=c(nObs,2,nT))
for (t in 1:nT){
    # X_{t+1} = X_t + VF(X_t) + FE + TE + noise
    X[,,t+1] = X[,,t] + 
               t(apply(X[,,t], 1, VF)) +  # Vector field dynamics
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

# Create panel data structure using plm  
panel_data <- pdata.frame(panel_data_df, index = c("id", "time"))



# Plot the true vector field for reference ----
cat("\n--- Generating True Vector Field Reference Plot ---\n")

# Create evaluation points (same method as in runPanelVFAnalysis)
X_unrolled <- aperm(X, c(1, 3, 2))
dim(X_unrolled) <- c(nObs * (nT+1), 2)  # X has nT+1 time periods
x_eval_true <- defineEvalPoints(X_unrolled, nEval)

# Compute true vector field on evaluation points
VF_true_eval <- t(apply(x_eval_true, 1, VF))

# Prepare data for plotting (using same settings as plotPanelVFAnalysis)
timeInterval <- 1
component_names <- c("X1", "X2")
lengthArrows <- 1
rescale <- FALSE  # Following the same setting as in the loop

# Create the plot
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
svg(file.path(output_dir, "True_VectorField.svg"), width = 8, height = 8)
par(mar = c(4, 4, 3, 2))

# Plot the true vector field
plot(x_eval_true, type = "n",
     xlab = component_names[1],
     ylab = component_names[2],
     main = "True Vector Field - Non-Linear Single Well")

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
X0_plot <- X[,,2]  # Final conditions
points(X0_plot[, 1], X0_plot[, 2], pch = 19, col = "red", cex = 0.5)

# Add density contours if sm package is available
if (requireNamespace("sm", quietly = TRUE)) {
  X0_complete <- X0_plot[complete.cases(X0_plot), ]
  if(nrow(X0_complete) > 1) {
    tryCatch({
      est.dens <- sm::sm.density(X0_complete, display = "none")
      contour(est.dens$eval.points[, 1], est.dens$eval.points[, 2], est.dens$estimate,
              add = TRUE, col = "purple")
    }, error = function(e) {
      # Skip density contours if they fail
    })
  }
}

dev.off()

cat(sprintf("True vector field plot saved to: %s\n", file.path(output_dir, "True_VectorField.svg")))
cat("Blue arrows: True non-linear single well vector field\n")
cat("Red points: Initial conditions\n")
cat("Purple contours: Data density\n\n")

## Panel Vector Field Analysis (Modern Approach) ----
# Test different kernel types and parameters following VFLifeExpectancyFE.R structure
kernels_to_test <- c("epa")
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
                                             uniform_weights = TRUE,
                                             kernel.type = kernel,
                                             method.h = "silverman",
                                             chunk_size = 1024,
                                             bootstrap_B = 1, 
                                             ridge_param = ridge,
                                             rcond_threshold = rcond_thresh)
      
      # Plot results using the modern approach
      plotPanelVFAnalysis(analysis_results,
                          timeInterval = 1,
                          rescale = FALSE,
                          rescale_ref_index = NULL,
                          years = 1:(nT+1),
                          label_names = paste0("unit_", 1:nObs),
                          component_names = c("X1", "X2"),
                                                                    save_plots = TRUE,
                                          save_path = output_dir)
      
      # --- Comparison Analysis: True vs Estimated ---
      
      # Ensure output directory exists
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      
      # 1) Vector field comparison on evaluation points
      x_eval <- analysis_results$x
      VF_hat <- analysis_results$VF_hat
      
      # Compute true vector field on evaluation points
      VF_true <- t(apply(x_eval, 1, VF))
      
      # Calculate vector field differences
      VF_diff <- VF_hat - VF_true
      
      # Calculate norm (magnitude) of the vector difference
      VF_norm_diff <- sqrt(VF_diff[,1]^2 + VF_diff[,2]^2)
      
      # Remove NaN and infinite values from norm differences
      VF_norm_diff[is.nan(VF_norm_diff) | is.infinite(VF_norm_diff)] <- NA
      
      # Create heatmap of norm differences
      # The evaluation points form a regular grid: round(sqrt(nEval)) x round(sqrt(nEval))
      grid_size <- round(sqrt(nEval))
      actual_nEval <- grid_size^2
      
      # Extract grid coordinates for proper heatmap visualization
      x_unique <- unique(x_eval[,1])
      y_unique <- unique(x_eval[,2])
      x_sorted <- sort(x_unique)
      y_sorted <- sort(y_unique)
      
      # Reshape norm differences into grid format for visualization
      VF_norm_diff_grid <- matrix(VF_norm_diff[1:actual_nEval], nrow = grid_size, ncol = grid_size, byrow = TRUE)
      
      # Check for problematic values in the grid
      n_na_grid <- sum(is.na(VF_norm_diff_grid))
      if (n_na_grid > 0) {
        cat(sprintf("  - Note: %d/%d grid points contain NA values\n", n_na_grid, actual_nEval))
      }
      
      svg(file.path(output_dir, paste0("VF_NormDiff_", filename)), width = 8, height = 8)
      par(mar = c(4, 4, 3, 2))
      
      # Single heatmap showing norm of vector difference
      image(x_sorted, y_sorted, VF_norm_diff_grid, 
            main = "Vector Field Norm Difference ||Estimated - True||",
            xlab = "X1", ylab = "X2", col = heat.colors(100, rev = TRUE),
            useRaster = TRUE)
      
      # Add contour lines only if there's sufficient variation in the data
      points_range <- range(VF_norm_diff_grid, na.rm = TRUE)
      if (!is.na(points_range[1]) && !is.na(points_range[2]) && 
          points_range[2] > points_range[1] && 
          (points_range[2] - points_range[1]) > .Machine$double.eps) {
        tryCatch({
          contour(x_sorted, y_sorted, VF_norm_diff_grid, add = TRUE, 
                  nlevels = 6, col = "black", lwd = 0.5)
        }, error = function(e) {
          # If contour fails, just skip it silently
          cat("  - Warning: Could not add contour lines to heatmap\n")
        })
      }
      
      # Add range information
      mtext(sprintf("Range: [%.6f, %.6f]", points_range[1], points_range[2]), 
            side = 1, line = 2.5, cex = 0.8, col = "darkblue")
      
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
      
      # Save summary statistics to file in a more readable format
      stats_filename <- file.path(output_dir, paste0("Summary_Stats_", gsub("\\.svg$", ".txt", filename)))
      
      # Create a nicely formatted text output
      cat("=== VECTOR FIELD ESTIMATION COMPARISON STATISTICS ===\n\n", file = stats_filename)
      cat("Parameters:\n", file = stats_filename, append = TRUE)
      cat(sprintf("  Kernel: %s\n", kernel), file = stats_filename, append = TRUE)
      cat(sprintf("  Ridge parameter: %g\n", ridge), file = stats_filename, append = TRUE)
      cat(sprintf("  Condition threshold: %g\n", rcond_thresh), file = stats_filename, append = TRUE)
      cat(sprintf("  Number of observations: %d\n", nObs), file = stats_filename, append = TRUE)
      cat(sprintf("  Number of time periods: %d\n", nT), file = stats_filename, append = TRUE)
      cat(sprintf("  Number of evaluation points: %d\n", nEval), file = stats_filename, append = TRUE)
      cat("\n", file = stats_filename, append = TRUE)
      
      # Vector Field Statistics
      cat("VECTOR FIELD NORM DIFFERENCES ||Estimated - True||:\n", file = stats_filename, append = TRUE)
      cat(sprintf("  Valid observations: %d / %d (%.1f%%)\n",
                  stats_summary$VF_norm_differences$n_valid,
                  stats_summary$VF_norm_differences$n_total,
                  100 * stats_summary$VF_norm_differences$n_valid / stats_summary$VF_norm_differences$n_total), 
          file = stats_filename, append = TRUE)
      cat(sprintf("  Mean: %12.8f | SD: %12.8f | RMSE: %12.8f | MAE: %12.8f\n",
                  stats_summary$VF_norm_differences$mean,
                  stats_summary$VF_norm_differences$sd,
                  stats_summary$VF_norm_differences$rmse,
                  stats_summary$VF_norm_differences$mae), file = stats_filename, append = TRUE)
      cat(sprintf("  Min: %12.8f | Max: %12.8f\n",
                  stats_summary$VF_norm_differences$min,
                  stats_summary$VF_norm_differences$max), file = stats_filename, append = TRUE)
      
      # Add warning if many invalid values
      if (stats_summary$VF_norm_differences$n_valid < 0.9 * stats_summary$VF_norm_differences$n_total) {
        cat("  WARNING: Significant number of NaN/Inf values detected!\n", file = stats_filename, append = TRUE)
      }
      cat("\n", file = stats_filename, append = TRUE)
      
      # Fixed Effects Statistics (if available)
      if (!is.null(stats_summary$FE_norm_differences)) {
        cat("FIXED EFFECTS NORM DIFFERENCES ||Estimated - True||:\n", file = stats_filename, append = TRUE)
        cat(sprintf("  Valid observations: %d / %d (%.1f%%)\n",
                    stats_summary$FE_norm_differences$n_valid,
                    stats_summary$FE_norm_differences$n_total,
                    100 * stats_summary$FE_norm_differences$n_valid / stats_summary$FE_norm_differences$n_total), 
            file = stats_filename, append = TRUE)
        cat(sprintf("  Mean: %12.8f | SD: %12.8f | RMSE: %12.8f | MAE: %12.8f\n",
                    stats_summary$FE_norm_differences$mean,
                    stats_summary$FE_norm_differences$sd,
                    stats_summary$FE_norm_differences$rmse,
                    stats_summary$FE_norm_differences$mae), file = stats_filename, append = TRUE)
        cat(sprintf("  Min: %12.8f | Max: %12.8f\n",
                    stats_summary$FE_norm_differences$min,
                    stats_summary$FE_norm_differences$max), file = stats_filename, append = TRUE)
        
        # Add warning if many invalid values
        if (stats_summary$FE_norm_differences$n_valid < 0.9 * stats_summary$FE_norm_differences$n_total) {
          cat("  WARNING: Significant number of NaN/Inf values detected!\n", file = stats_filename, append = TRUE)
        }
        cat("\n", file = stats_filename, append = TRUE)
      } else {
        cat("FIXED EFFECTS: Not estimated\n\n", file = stats_filename, append = TRUE)
      }
      
      # Time Effects Statistics (if available)
      if (!is.null(stats_summary$TE_norm_differences)) {
        cat("TIME EFFECTS NORM DIFFERENCES ||Estimated - True||:\n", file = stats_filename, append = TRUE)
        cat(sprintf("  Valid observations: %d / %d (%.1f%%)\n",
                    stats_summary$TE_norm_differences$n_valid,
                    stats_summary$TE_norm_differences$n_total,
                    100 * stats_summary$TE_norm_differences$n_valid / stats_summary$TE_norm_differences$n_total), 
            file = stats_filename, append = TRUE)
        cat(sprintf("  Mean: %12.8f | SD: %12.8f | RMSE: %12.8f | MAE: %12.8f\n",
                    stats_summary$TE_norm_differences$mean,
                    stats_summary$TE_norm_differences$sd,
                    stats_summary$TE_norm_differences$rmse,
                    stats_summary$TE_norm_differences$mae), file = stats_filename, append = TRUE)
        cat(sprintf("  Min: %12.8f | Max: %12.8f\n",
                    stats_summary$TE_norm_differences$min,
                    stats_summary$TE_norm_differences$max), file = stats_filename, append = TRUE)
        
        # Add warning if many invalid values
        if (stats_summary$TE_norm_differences$n_valid < 0.9 * stats_summary$TE_norm_differences$n_total) {
          cat("  WARNING: Significant number of NaN/Inf values detected!\n", file = stats_filename, append = TRUE)
        }
        cat("\n", file = stats_filename, append = TRUE)
      } else {
        cat("TIME EFFECTS: Not estimated\n\n", file = stats_filename, append = TRUE)
      }
      
      cat(sprintf("Generated at: %s\n", Sys.time()), file = stats_filename, append = TRUE)
      
      cat(sprintf("  - Norm heatmap saved to: VF_NormDiff_%s\n", filename))
      cat(sprintf("  - Statistics saved to: %s\n", basename(stats_filename)))
      
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





