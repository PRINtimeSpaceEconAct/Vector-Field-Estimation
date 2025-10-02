# ============================================================================
# Vector Field Estimation Test Script
# ============================================================================
# This script tests vector field estimation methods on synthetic data.
# 
# What it does:
# 1. Generates synthetic 2D data following a known vector field (configurable)
# 2. Estimates the vector field using kernel methods (LL or NW, adaptive or fixed)
# 3. Optionally performs bootstrap for significance testing
# 4. Compares estimated field with true field and VAR-based linear estimates
# 5. Produces visualization plots including:
#    - True vector field
#    - Estimated vector field (all arrows)
#    - Estimated vector field (significant arrows only, if bootstrap enabled)
#    - Error vector field
#    - Absolute and relative error heatmaps (log scale)
#    - VAR-based estimates and their errors
#
# Results saved to: test_pics_VF/ directory as PDF files with names based on settings
#   Format: {VF_TYPE}_{ESTIMATOR_TYPE}_{Adaptive/Fixed}_*.pdf
# ============================================================================

# Clear workspace and load dependencies
setwd("~/Library/CloudStorage/OneDrive-UniversityofPisa/timeSpaceEvolutionEcAct/RVF/R code/Vector Field Estimation/")
rm(list = ls())
DEBUG = FALSE
source("src/libs/loadLib.R")
library(fields)
library(latex2exp)

# ============================================================================
# USER CONFIGURATION
# ============================================================================

# Choose vector field: "double_well", "single_well", or "rotation"
VF_TYPE = "double_well"

# Choose estimator: "LL" (Local Linear) or "NW" (Nadaraya-Watson)
ESTIMATOR_TYPE = "LL"

# Use adaptive bandwidth?
USE_ADAPTIVE = TRUE

# Bootstrap settings
RUN_BOOTSTRAP = FALSE  # Set to TRUE to run bootstrap (can be slow)
N_BOOTSTRAP = 100      # Number of bootstrap replicates

# ============================================================================

# Create output directory if it doesn't exist ----
output_dir = "src/tests/test_pics_VF"
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
}

# parameters ----
nObs = 1000
nEval = 2500

# data Generation ----
set.seed(1)
X0 = mvrnorm(nObs, mu=c(0,0),Sigma = 0.5*diag(2))

# Define vector field based on user choice ----
if (VF_TYPE == "double_well") {
    # example 1 - double well
    VF <- function(X){
        # X = (x,y)
        # U(X) = x^4 - x^2 + y^2
        # VF(X) = -grad U(X) = -(4x^3 - 2x, 2y)
        return( -0.1*c(4*X[1]^3 - 2*X[1], 2*X[2]) )
    }
} else if (VF_TYPE == "single_well") {
    # example 2 -- single well
    VF <- function(X){
        # X = (x,y)
        # U(X) = x^2 + y^2
        # VF(X) = -grad U(X) = -(2x, 2y)
        return( -0.01*c(2*X[1], 2*X[2]) )
    }
} else if (VF_TYPE == "rotation") {
    # example 3 -- rotation
    M = matrix(c(cos(pi/2), -sin(pi/2), sin(pi/2), cos(pi/2)),nrow=2,ncol=2)
    VF <- function(X){
        # X = (x,y), theta = pi/2
        return (M %*% X)
    }
} else {
    stop("Invalid VF_TYPE. Choose 'double_well', 'single_well', or 'rotation'.")
}

# Create filename prefix based on settings ----
adaptive_str = if(USE_ADAPTIVE) "Adaptive" else "Fixed"
filename_prefix = paste0(VF_TYPE, "_", ESTIMATOR_TYPE, "_", adaptive_str)

# apply VF
X1 = X0 + t(apply(X0, 1, VF)) +  matrix(rnorm(2*nObs),nrow=nObs) %*% matrix(c(0.01,0.005,0.005,0.02),nrow=2)

# eval points
X_unrolled <- rbind(X0, X1)
x <- defineEvalPoints(X_unrolled, nEval)

# stima ----

# Choose estimator based on user configuration
if (ESTIMATOR_TYPE == "LL" && USE_ADAPTIVE) {
    est_field = LLfieldAdaptive(X0, X1, x=x, kernel.type="epa", method.h = "silverman",
                                chunk_size=1000, gc=TRUE, hOpt = TRUE, alphaOpt = TRUE)
} else if (ESTIMATOR_TYPE == "LL" && !USE_ADAPTIVE) {
    est_field = LLfield(X0, X1, x=x, kernel.type="epa", method.h = "silverman",
                       chunk_size=1000, gc=TRUE)
} else if (ESTIMATOR_TYPE == "NW" && USE_ADAPTIVE) {
    est_field = NWfieldAdaptive(X0, X1, x=x, kernel.type="epa", method.h = "silverman",
                                chunk_size=1000, gc=TRUE, hOpt = TRUE, alphaOpt = TRUE)
} else if (ESTIMATOR_TYPE == "NW" && !USE_ADAPTIVE) {
    est_field = NWfield(X0, X1, x=x, kernel.type="epa", method.h = "silverman",
                       chunk_size=1000, gc=TRUE)
} else {
    stop("Invalid ESTIMATOR_TYPE. Choose 'LL' or 'NW'.")
}

cat("Using:", VF_TYPE, "vector field with", ESTIMATOR_TYPE, 
    if(USE_ADAPTIVE) "adaptive" else "fixed", "bandwidth\n")

# Bootstrap for significance testing ----
if (RUN_BOOTSTRAP) {
    cat("\nRunning bootstrap for significance testing...\n")
    est_field_bootstrap = bootstrapKernelFieldErrors(est_field, B = N_BOOTSTRAP, 
                                                      chunk_size = nrow(est_field$x))
    signifBoot = significanceBootstrap(est_field$estimator, est_field_bootstrap, p_crit = 0.05)
    cat("Bootstrap completed:", sum(signifBoot), "out of", length(signifBoot), "points are significant\n")
} else {
    signifBoot = NULL
    cat("\nBootstrap skipped (set RUN_BOOTSTRAP = TRUE to enable)\n")
}

# plot ----
# plot campo vero ----
dev.new()
op <- par(family = "mono") #Possible families: "mono", "Helvetica","Palatino" or "Times"
VFx = t(apply(x, 1, VF))
lengthArrows=0.1
plot(x, type = "n", xlab = TeX(r'($X_1$)'), ylab=TeX(r'($X_2$)'), main = "")
arrows(x[,1],x[,2],x[,1]+lengthArrows*VFx[,1],x[,2]+lengthArrows*VFx[,2],angle=15,col="black",length=0.05)
abline(h=0)
abline(v=0)
dev.copy2pdf(file=paste0(output_dir, "/", filename_prefix, "_campoVero.pdf"))
# par(op)

## plot campo stimato ----
dev.new()
op <- par(family = "mono") #Possible families: "mono", "Helvetica","Palatino" or "Times" 
signifVFest = significanceVF(est_field)
lengthArrows=0.1
plot(x, type = "n", xlab = TeX(r'($X_1$)'), ylab=TeX(r'($X_2$)'), main = "Tutte frecce")
arrows(est_field$x[,1], est_field$x[,2],
       est_field$x[,1] + lengthArrows*est_field$estimator[,1], 
       est_field$x[,2] + lengthArrows*est_field$estimator[,2],
       length = 0.05, angle = 15, col = "blue")
abline(h=0)
abline(v=0)
dev.copy2pdf(file=paste0(output_dir, "/", filename_prefix, "_campoStimato.pdf"))
par(op)

## plot campo stimato significativo ----
if (RUN_BOOTSTRAP && !is.null(signifBoot)) {
    signifInd = which(signifBoot)
    if (length(signifInd) > 0) {
        dev.new()
        op <- par(family = "mono") #Possible families: "mono", "Helvetica","Palatino" or "Times" 
        lengthArrows=0.1
        plot(x, type = "n", xlab = TeX(r'($X_1$)'), ylab=TeX(r'($X_2$)'), main = "Solo significative")
        arrows(est_field$x[signifInd,1], est_field$x[signifInd,2],
               est_field$x[signifInd,1] + lengthArrows*est_field$estimator[signifInd,1], 
               est_field$x[signifInd,2] + lengthArrows*est_field$estimator[signifInd,2],
               length = 0.05, angle = 15, col = "blue")
        abline(h=0)
        abline(v=0)
        dev.copy2pdf(file=paste0(output_dir, "/", filename_prefix, "_campoStimatoSignif.pdf"))
        par(op)
    } else {
        cat("\nNo significant points found to plot.\n")
    }
}


## plot errore ----
VFx = t(apply(x, 1, VF))
plot(est_field$x, type = "n", xlab = "X", ylab = "Y", main = "Error Vector Field")
arrows(est_field$x[,1], est_field$x[,2],
       est_field$x[,1] + est_field$estimator[,1] - VFx[,1],
       est_field$x[,2] + est_field$estimator[,2] - VFx[,2],
       length = 0.05, angle = 15, col = "red")

## image errore assoluto ----
errorNorm = sqrt((est_field$estimator[,1] - VFx[,1])^2 + (est_field$estimator[,2] - VFx[,2])^2)
image.plot(x = unique(x[,1]), y = unique(x[,2]), z = matrix(log10(errorNorm), nrow=sqrt(nEval), ncol=sqrt(nEval)),xlab="x",ylab="y",main="error norm (log10)")

## image errore relativo ----
VFx = t(apply(x, 1, VF))
dev.new()
op <- par(family = "mono") #Possible families: "mono", "Helvetica","Palatino" or "Times" 
errorNormRel = sqrt((est_field$estimator[,1] - VFx[,1])^2 + (est_field$estimator[,2] - VFx[,2])^2)/sqrt(VFx[,1]^2 + VFx[,2]^2)
errorNormRel[is.na(errorNormRel)] = 1; errorNormRel[!is.finite(errorNormRel)] = 1; errorNormRel[errorNormRel==0] = 1
image.plot(x = unique(x[,1]), y = unique(x[,2]), z = matrix(log10(errorNormRel), nrow=sqrt(nEval), ncol=sqrt(nEval)),xlab=TeX(r'($X_1$)'),ylab=TeX(r'($X_2$)'),main="")
abline(h=0)
abline(v=0)
dev.copy2pdf(file=paste0(output_dir, "/", filename_prefix, "_LogNormRel.pdf"))
par(op)


# stima VAR ----
data = data.frame(X0 = X0[,1], Y0 = X0[,2], X1 = X1[,1], Y1 = X1[,2])
lm1 = lm(X1 ~ X0 + Y0 + 0, data = data)
lm2 = lm(Y1 ~ X0 + Y0 + 0, data = data)
A = matrix(c(coef(lm1)[1:2],coef(lm2)[1:2]),nrow=2,byrow=TRUE)

VFVAR = function(x){
    return((A-diag(2)) %*% x)
}
VFVARx = t(apply(x, 1, VFVAR))

dev.new()
op <- par(family = "mono") #Possible families: "mono", "Helvetica","Palatino" or "Times" 
plot(x, type = "n", xlab = "X1", ylab ="X2", main = "")
arrows(x[,1], x[,2],
       x[,1] +  VFVARx[,1], 
       x[,2] +  VFVARx[,2],
       length = 0.05, angle = 15, col = "blue")
abline(h=0)
abline(v=0)
dev.copy2pdf(file=paste0(output_dir, "/", VF_TYPE, "_VAR_campoStimato.pdf"))
par(op)

## image errore relativo ----
dev.new()
op <- par(family = "mono") #Possible families: "mono", "Helvetica","Palatino" or "Times" 
VFx = t(apply(x, 1, VF))
errorNormRel = sqrt((VFVARx[,1] - VFx[,1])^2 + (VFVARx[,2] - VFx[,2])^2)/sqrt(VFx[,1]^2 + VFx[,2]^2)
image.plot(x = unique(x[,1]), y = unique(x[,2]), z = matrix(log10(errorNormRel), nrow=sqrt(nEval), ncol=sqrt(nEval)),xlab=TeX(r'($X_1$)'),ylab=TeX(r'($X_2$)'),main="")
abline(h=0)
abline(v=0)
dev.copy2pdf(file=paste0(output_dir, "/", VF_TYPE, "_VAR_LogNormRel.pdf"))
par(op)


