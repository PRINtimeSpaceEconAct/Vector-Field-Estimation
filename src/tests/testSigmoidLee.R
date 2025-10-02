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
nObs = 100  # Number of cross-sectional units
nT = 20     # Number of time periods for transitions
nEval = 2500  # Number of evaluation points for VF estimation
output_dir <- "test_pics_Sigmoid"

lengthArrows <- 0.1

# Data generation ----
set.seed(1)  # For reproducibility

# Generate X_i,t from N(2, 4^2) as specified in the paper
X_flat <- mvrnorm(n = nObs * nT, mu = c(2, 2), Sigma = 16 * diag(2))
X <- array(X_flat, dim = c(nObs, nT, 2))
X <- aperm(X, c(1, 3, 2)) # Reshape to nObs x 2 x nT

# Generate fixed effects (FE) based on Lee (2022) specification
# eta_i = 0.5 * X_bar_i + v_i
# Calculate time-average X_bar_i for each unit i from the generated X
X_bar_i <- apply(X, c(1, 2), mean)

# Generate v_i from U(-0.5, 0.5)
v_i <- matrix(runif(nObs * 2, -0.5, 0.5), nrow = nObs, ncol = 2)

# Calculate fixed effects alpha_i (eta_i in the paper)
alpha_i <- 0.5 * X_bar_i + v_i

# Generate time effects (TE) - common time-varying factors (not in paper, set to 0)
gamma_t = matrix(0,nrow=nT,ncol=2)


# Sigmoid vector field from Lee (2022) ----
VF <- function(X){
    # m(x) = (1 + exp(-x))^-1
    # VF(X) acts on each component
    return( 1 / (1 + exp(-X)) )
}


VF_jacobian <- function(X){
    # m(x) = (1 + exp(-x))^-1
    # m'(x) = exp(-x) / (1 + exp(-x))^2
    # Jacobian is a diagonal matrix with m'(X_i) on the diagonal
    m_prime_x1 <- exp(-X[1]) / (1 + exp(-X[1]))^2
    m_prime_x2 <- exp(-X[2]) / (1 + exp(-X[2]))^2
    
    J_11 <- m_prime_x1
    J_12 <- 0
    J_21 <- 0
    J_22 <- m_prime_x2
    return(matrix(c(J_11, J_12, J_21, J_22), nrow=2, byrow=TRUE))
}

# Generate exogenous Y = VF(X) + FE + TE + noise ----
# u_i,t is drawn from N(0, 1)
noise = array(rnorm(nObs*nT*2, sd = 1), dim=c(nObs,2,nT))
Y_exog <- array(NA, dim = c(nObs, 2, nT))
for (t in 1:nT) {
    # Y is based on X at time t. X has dimensions (nObs, 2, nT)
    # VF acts on each row of X[,,t]
    Y_exog[,,t] <- t(apply(X[,,t, drop=FALSE], 1, VF)) +
                   alpha_i +
                   array(rep(gamma_t[t,], each=nObs), dim=c(nObs,2)) +
                   noise[,,t]
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
X_unrolled <- aperm(X[,,1:nT, drop=FALSE], c(1, 3, 2)) # Use nT periods
dim(X_unrolled) <- c(nObs * nT, 2)
x_eval <- defineEvalPoints(X_unrolled, nEval)

# Compute true vector field on evaluation points
VF_true_eval <- t(apply(x_eval, 1, VF))

# Prepare data for plotting (using same settings as plotPanelVFAnalysis)
timeInterval <- 1
component_names <- c("X1", "X2")
rescale <- FALSE  # Following the same setting as in the loop

# Create the plot
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
svg(file.path(output_dir, "True_VectorField.svg"), width = 8, height = 8)
par(mar = c(4, 4, 3, 2))

# Plot the true vector field
plot(x_eval, type = "n",
     xlab = component_names[1],
     ylab = component_names[2],
     main = "True Vector Field - Sigmoid")

if(rescale) {
  abline(h = 1, lty = 3)
  abline(v = 1, lty = 3)
} else {
  abline(h = 0, lty = 3)
  abline(v = 0, lty = 3)
}

# Draw vector field arrows
arrows(x_eval[, 1], x_eval[, 2],
       x_eval[, 1] + lengthArrows * VF_true_eval[, 1],
       x_eval[, 2] + lengthArrows * VF_true_eval[, 2],
       angle = 15, col = "blue", length = 0.05)

# Add observed data points
X0_plot <- X[,,1] # Initial conditions
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
cat("Blue arrows: True sigmoid vector field\n")
cat("Red points: Initial conditions\n")
cat("Purple contours: Data density\n\n")

## Panel Vector Field Analysis (Modern Approach) ----
# Test parameters from the paper
kernels_to_test <- c("epa")
c_h_to_test <- c(0.5, 1.0, 1.5, 2.0)
ridge <- 1e-6 # Keep ridge/rcond constant as they are not varied in the paper
rcond_thresh <- 1e-4

# Create evaluation points once
X_unrolled <- aperm(X, c(1, 3, 2)) # Use nT periods
dim(X_unrolled) <- c(nObs * nT, 2)
x_eval <- defineEvalPoints(X_unrolled, nEval)

# Calculate bandwidth component from data variance. Paper uses {var(X_i,t)}^(1/2).
# We use the root of the mean variance of the components.
sd_X <- sqrt(mean(apply(X_unrolled, 2, var)))

total_iterations <- length(kernels_to_test) * length(c_h_to_test)
current_iteration <- 0

for (kernel in kernels_to_test) {
  for (c_h in c_h_to_test) {
      current_iteration <- current_iteration + 1
    
    # Calculate bandwidth h based on the paper's formula
    h <- c_h * sd_X * (nObs * nT)^(-1/7)
      
      kernel_type_str <- tools::toTitleCase(kernel)
    ch_str <- gsub("\\.", "_", format(c_h, nsmall = 1))
    filename_base <- sprintf("Sigmoid_%s_ch_%s", kernel_type_str, ch_str)
    
    cat(sprintf("\n--- [ %d / %d ] Processing: %s (h=%.4f) ---\n", current_iteration, total_iterations, filename_base, h))

    # Run panel VF analysis using the exogenous approach
    panel_vf_results <- estimate_panel_vf_exogenous(
        X = X, # Use nT periods for regressors
        Y = Y_exog,               # The exogenous outcome
        x = x_eval,
        nEval = nEval,
        FE = TRUE,
        TE = TRUE,
        uniform_weights = FALSE,
        kernel.type = kernel,
        method.h = NULL, # Provide h manually
        h = h,
        chunk_size = 1024,
        gc = FALSE,
        ridge_param = ridge,
        rcond_threshold = rcond_thresh
    )
    
    # --- Plotting: True vs Estimated Vector Fields ---
    
      # Estimated VF from results
      VF_hat <- panel_vf_results$estimator
      
      # Compute true vector field on the same evaluation points
      VF_true <- t(apply(x_eval, 1, VF))
      
      # Calculate relative difference for heatmap
      VF_diff <- VF_hat - VF_true
      VF_norm_diff <- sqrt(VF_diff[,1]^2 + VF_diff[,2]^2)
      VF_norm_true <- sqrt(VF_true[,1]^2 + VF_true[,2]^2)
      VF_norm_rel_diff <- VF_norm_diff / VF_norm_true
      VF_norm_rel_diff[is.na(VF_norm_rel_diff) | !is.finite(VF_norm_rel_diff)] <- NA
      VF_norm_rel_diff_log <- log10(VF_norm_rel_diff)
      VF_norm_rel_diff_log[is.infinite(VF_norm_rel_diff_log)] <- NA
      
      grid_size <- round(sqrt(nEval))
      actual_nEval <- grid_size^2
      x_sorted <- sort(unique(x_eval[,1]))
      y_sorted <- sort(unique(x_eval[,2]))
      VF_norm_rel_diff_grid <- matrix(VF_norm_rel_diff_log[1:actual_nEval], nrow = grid_size, ncol = grid_size, byrow = TRUE)
      
      # Estimate data density for contour plots
      X_unrolled_complete <- X_unrolled[complete.cases(X_unrolled), ]
      est.dens <- NULL
      if (requireNamespace("sm", quietly = TRUE) && nrow(X_unrolled_complete) > 1) {
          tryCatch({
              est.dens <- sm::sm.density(X_unrolled_complete, display = "none")
          }, error = function(e) {
              cat("Error during density estimation:\n")
              print(e)
          })
      }
      
      # Create a single SVG with three plots
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      svg_filename <- paste0("VF_Comparison_", filename_base, ".svg")
      svg(file.path(output_dir, svg_filename), width = 18, height = 6)
      
      layout(matrix(c(1, 2, 3), nrow = 1), widths = c(1, 1, 1.4))
      par(mar = c(4, 4, 3, 2))
      
      # Plot 1: True Vector Field
      plot(x_eval, type = "n",
                 xlab = "X1", ylab = "X2",
           main = "True Vector Field")
      arrows(x_eval[, 1], x_eval[, 2],
             x_eval[, 1] + lengthArrows * VF_true[, 1],
             x_eval[, 2] + lengthArrows * VF_true[, 2],
             angle = 15, col = "blue", length = 0.05)
      if (!is.null(est.dens)) {
          contour(est.dens$eval.points[, 1], est.dens$eval.points[, 2], est.dens$estimate,
                  add = TRUE, col = "purple")
      }
      
      # Plot 2: Estimated Vector Field
      plot(x_eval, type = "n",
           xlab = "X1", ylab = "X2",
           main = "Estimated Vector Field")
      arrows(x_eval[, 1], x_eval[, 2],
             x_eval[, 1] + lengthArrows * VF_hat[, 1],
             x_eval[, 2] + lengthArrows * VF_hat[, 2],
             angle = 15, col = "red", length = 0.05)
      if (!is.null(est.dens)) {
          contour(est.dens$eval.points[, 1], est.dens$eval.points[, 2], est.dens$estimate,
                  add = TRUE, col = "purple")
      }
      
      # Plot 3: Relative Difference Heatmap
      par(mar = c(4, 4, 3, 3)) # Adjust margin for color bar
      image.plot(x_sorted, y_sorted, VF_norm_rel_diff_grid, 
                 main = "Relative Difference (log10)",
                 xlab = "X1", ylab = "X2")
      
      dev.off()
      
      cat(sprintf("  - Vector field comparison plot saved to: %s\n", svg_filename))
                          
  }
}
