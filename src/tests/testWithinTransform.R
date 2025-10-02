# Visualize the effect of within_transform on simple one-step VF dynamics

rm(list = ls())
DEBUG <- FALSE

source("src/libs/loadLib.R")

# Parameters
output_dir <- file.path("test_pics_withinTransform")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

grid_side <- 20         # number of points per axis (grid_side x grid_side)
x_lim <- c(-2, 2)
y_lim <- c(-2, 2)
dt <- 1                  # single step size

# Within-transform toggles
FE <- TRUE
TE <- TRUE
uniform_weights <- FALSE  # set to TRUE for uniform weights, FALSE for kernel weights

# Arrow style
arrow_head_len <- 0.05
point_cex <- 0.58

# Vector field VF as in fixEffectsLee.R (double well example)
VF <- function(X) {
    # X = (x, y)
    # U(X) = y^4 - y^2 + x^2
    # VF(X) = -grad U(X) = -(2x, 4y^3 - 2y)
    -0.1 * c(2 * X[1], 4 * X[2]^3 - 2 * X[2])
}

# Build equispaced grid
x_seq <- seq(x_lim[1], x_lim[2], length.out = grid_side)
y_seq <- seq(y_lim[1], y_lim[2], length.out = grid_side)
grid_xy <- as.matrix(expand.grid(x_seq, y_seq))
colnames(grid_xy) <- c("X1", "X2")

# Geometric minimal h for square, equispaced grid (interior centers only)
s <- x_seq[2] - x_seq[1]  # common spacing on both axes
h_4_interior_geom <- s
h_8_interior_geom <- sqrt(2) * s
cat(sprintf("Geometric h (interior centers) for >=4 neighbors: %.10f\n", h_4_interior_geom))
cat(sprintf("Geometric h (interior centers) for >=8 neighbors: %.10f\n", h_8_interior_geom))

nObs <- nrow(grid_xy)
nT <- 3  # need at least 3 to avoid dimension drop inside within_transform

# Create panel array X: [nObs, 2, nT]
X <- array(NA_real_, dim = c(nObs, 2, nT))
X[, , 1] <- grid_xy

# Two steps of dynamics using VF
VF_step1 <- t(apply(X[, , 1], 1, VF))  # nObs x 2
X[, , 2] <- X[, , 1] + dt * VF_step1
VF_step2 <- t(apply(X[, , 2], 1, VF))  # nObs x 2
X[, , 3] <- X[, , 2] + dt * VF_step2

# Use weights evaluated at the initial grid points (ignored if uniform)
x_eval <- X[, , 1]
nEval_chunk <- nrow(x_eval)

# Apply within_transform to both X0 and Y = X1 - X0
h_used <- h_8_interior_geom
if (uniform_weights) {
    wt_label <- "uniform weights"
    title_suffix <- sprintf(" (%s)", wt_label)
    file_tag <- "uniform"
} else {
    wt_label <- "kernel weights"
    title_suffix <- sprintf(" (h=%.6g, %s)", h_used, wt_label)
    h_numeric_tag <- gsub("[^0-9a-zA-Z.-]", "_", sprintf("%.6g", h_used))
    file_tag <- paste0("h_", h_numeric_tag)
}

Filtered <- within_transform(X,
                             FE = FE,
                             TE = TE,
                             uniform_weights = uniform_weights,
                             nEval_chunk = nEval_chunk,
                             h = h_used,
                             x = x_eval)

# Extract transformed points (X0_star) and transformed arrows (Y_star)
# Rebuild unrolled arrays to [time, unit, ...] and take t = 1 across all units
tCount <- nT - 1
nEval_used <- dim(Filtered$X0_star_unrolled)[3]
X0_star_arr <- array(Filtered$X0_star_unrolled, dim = c(tCount, nObs, 2, nEval_used))
X0_star <- matrix(NA_real_, nObs, 2)
for (i in seq_len(nObs)) {
    X0_star[i, ] <- X0_star_arr[1, i, , i]
}
Y1_star_arr <- array(Filtered$Y1_star_unrolled, dim = c(tCount, nObs, nEval_used))
Y2_star_arr <- array(Filtered$Y2_star_unrolled, dim = c(tCount, nObs, nEval_used))
Y_star <- matrix(NA_real_, nObs, 2)
for (i in seq_len(nObs)) {
    Y_star[i, ] <- c(Y1_star_arr[1, i, i], Y2_star_arr[1, i, i])
}

# End points after transform
X1_star <- X0_star + Y_star
Y_original <- X[, , 2] - X[, , 1]

# Combined Plot: Original vs. Transformed
combined_file <- file.path(output_dir, sprintf("WithinTransform_Combined_VF_Arrows_%s.svg", file_tag))
svg(combined_file, width = 16, height = 8)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))

# Plot 1: Original
plot(X[, , 1], type = "n", xlab = "X1", ylab = "X2", main = "Original grid: one VF step")
abline(h = 0, v = 0, lty = 3)
arrows(X[, 1, 1], X[, 2, 1], X[, 1, 2], X[, 2, 2], angle = 15, length = arrow_head_len, col = "blue")

# Plot 2: Transformed
# Compute limits using both X0_star and X1_star so arrows are fully visible
xlim_star <- range(c(X0_star[, 1], X1_star[, 1]))
ylim_star <- range(c(X0_star[, 2], X1_star[, 2]))
plot(X0_star, type = "n", xlab = "X1 (transformed)", ylab = "X2 (transformed)",
     main = paste0("Within-transformed grid: transformed one-step arrows", title_suffix),
     xlim = xlim_star, ylim = ylim_star)
abline(h = 0, v = 0, lty = 3)
arrows(X0_star[, 1], X0_star[, 2], X1_star[, 1], X1_star[, 2], angle = 15, length = arrow_head_len, col = "red")

dev.off()

cat("Saved: ", combined_file, "\n", sep = "")

# Plot 3: Side-by-side comparison of point and VF shifts
point_file <- file.path(output_dir, sprintf("WithinTransform_Point_And_VF_Shift_%s.svg", file_tag))
svg(point_file, width = 16, height = 8)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))

# Left Plot: Position shift
plot(rbind(X[, , 1], X0_star), type = "n", xlab = "X1", ylab = "X2", main = paste0("Position Transformation (X vs X*)", title_suffix))
abline(h = 0, v = 0, lty = 3)
points(X[, , 1], pch = 19, col = "black", cex = point_cex)
points(X0_star, pch = 19, col = "red", cex = point_cex)
legend("topright", legend = c("X (Original)", "X* (Transformed)"), col = c("black", "red"), pch = 19)

# Right Plot: VF step shift
plot(rbind(Y_original, Y_star), type = "n", xlab = "X1", ylab = "X2", main = paste0("VF Step Transformation (Y vs Y*)", title_suffix))
abline(h = 0, v = 0, lty = 3)
points(Y_original, pch = 17, col = "blue", cex = point_cex)
points(Y_star, pch = 17, col = "magenta", cex = point_cex)
legend("topright", legend = c("Y (Original)", "Y* (Transformed)"), col = c("blue", "magenta"), pch = 17)

dev.off()

cat("Saved: ", point_file, "\n", sep = "")


