source("src/libs/panel.R")

#' Performs statistical significance testing for vector field estimations
#' 
#' @param est An object containing vector field estimation results
#' @param p_crit Critical p-value threshold for significance testing (default: 0.05)
#' 
#' @return A list containing:
#'   \item{Var}{Array of variance matrices for the estimations (2 x 2 x n, where n is number of evaluation points)}
#'   \item{ChiSquare_stat}{Vector of chi-square test statistics (length n)}
#'   \item{p_values}{Vector of p-values for the significance tests (length n)}
#'   \item{signif}{Binary vector of significance indicators (1 = significant, 0 = not significant, length n)}
significanceVF <- function(est,p_crit = 0.05){
    
    # Extract parameters from the estimation object
    h = est$h
    X0 = est$X0
    X1 = est$X1
    
    n = nrow(X0)
    Y = X1-X0
    Yhat = interp2d(X0,est$x,est$estimator)
    eps = Y - Yhat    
    
    # Calculate covariance matrix for homoskedastic errors
    Sigma = cov(eps, use = "complete.obs")

    # Set kernel constant based on kernel type
    if (est$kernel.type == "epa"){ k = 3/5 
    } else if (est$kernel == "gauss") { k = 1/(2*sqrt(pi)) }
    
    # Initialize arrays for results
    Var = array(data = NA, dim = c(2,2,nrow(x)))
    ChiSquare_stat = rep(NA, nrow(x))     
    p_values = rep(NA, nrow(x))
    signif = rep(NA, nrow(x))
    
    # Calculate variance, chi-square statistic, and significance for each estimation point
    for (i in 1:nrow(x)){
        Var[,,i] = (k^2 * Sigma / (est$density[i] * n*h^2)) 
        ChiSquare_stat[i] = t(est$estimator[i,]) %*% solve(Var[,,i]) %*% est$estimator[i,]
        p_values[i] = 1-pchisq(ChiSquare_stat[i], 2)
        p_values[i] = ifelse( is.na(p_values[i]), 1 , p_values[i])
        signif[i] = ifelse(p_values[i] < p_crit, 1, 0)
    }
    
    return(listN(Var,ChiSquare_stat,p_values,signif))
}

#' Performs bootstrap resampling to estimate errors in kernel-based vector field estimations
#' 
#' @param result An object containing vector field estimation results
#' @param B Number of bootstrap replicates (default: 500)
#' @param chunk_size Chunk size for processing (default: number of evaluation points)
#' 
#' @return A 3D array of bootstrapped field estimators with dimensions [nEval x 2 x B],
#'         where nEval is the number of evaluation points, 2 is for x and y components,
#'         and B is the number of bootstrap replicates
bootstrapKernelFieldErrors <- function(result, B = 500, chunk_size = nrow(result$x)) {
    # 1. Extract common parameters from result
    X0 = result$X0
    X1 = result$X1
    x_eval = result$x
    h_opt = result$h
    kernel.type = result$kernel.type
    method.h = result$method.h # Used for pilot density in adaptive case if h is not fixed
    type.est.base = sub("_adaptive", "", result$type.est) # Get "NW" or "LL"

    # 2. Determine if adaptive and get alpha
    # Check if lambda exists and is a vector (indicating adaptive bandwidths per observation)
    is_adaptive = !is.null(result$lambda) && is.vector(result$lambda) && length(result$lambda) > 1
    alpha_opt = if (is_adaptive && !is.null(result$alpha)) result$alpha else NULL

    # 3. Select the appropriate compute function based on type.est.base
    # Assumes these functions are available in the calling environment
    if (type.est.base == "NW") {
        computeFieldComponents <- computeNWFieldComponents
    } else if (type.est.base == "LL") {
        computeFieldComponents <- computeLLFieldComponents
    } else {
        stop("Unsupported estimation type for bootstrap: ", result$type.est)
    }

    # 4. Prepare data and storage
    Nobs = nrow(X0)
    Y = X1 - X0
    Y1 = Y[, 1]
    Y2 = Y[, 2]
    nEvalPoints = nrow(x_eval)
    estimators_array = array(NA, dim = c(nEvalPoints, 2, B))

    # 5. Setup progress bar
    # Assumes createCustomProgressBar is available
    show_progress = !exists("DEBUG") || !DEBUG
    if (show_progress) {
        cat("Running bootstrap with", B, "replicates (Type: ", result$type.est, ")...\n", sep="")
        pb = createCustomProgressBar(min = 0, max = B)
    } else {
        cat("Starting bootstrap with", B, "replicates (Type: ", result$type.est, ")...\n", sep="")
    }

    # 6. Bootstrap loop
    for (b in 1:B) {
        # Resample
        idx_star = sample(1:Nobs, Nobs, replace = TRUE)
        X0_star = X0[idx_star, ]
        Y1_star = Y1[idx_star]
        Y2_star = Y2[idx_star]

        if (show_progress) {
            pb$update(b, paste("Bootstrap iteration", b, "/", B))
        }

        # Calculate lambda for the resampled data if adaptive
        current_lambda = NULL
        if (is_adaptive) {
            # Need getLocalBandwidth function here
            # Assuming getLocalBandwidth is available in the environment
            # Use h_opt for pilot density calculation to maintain consistency with original estimation
            current_lambda = getLocalBandwidth(X0_star, kernel.type=kernel.type, D=NULL,
                                     method.h=NULL, # Explicitly pass h_opt
                                     h=h_opt,
                                     sparse=FALSE, # Assuming defaults consistent with NW/LL field functions
                                     gc=FALSE,
                                     chunk_size=Nobs, # Use Nobs for density calculation over resampled data
                                     alpha = alpha_opt)
        }

        # Re-estimate using the selected compute function
        # Pass method.h = NULL because h = h_opt is fixed for the main estimation step
        # Pass chunk_size for the estimation step itself
        est_star = computeFieldComponents(X0 = X0_star, Y1 = Y1_star, Y2 = Y2_star, x = x_eval,
            nEval = nEvalPoints, kernel.type = kernel.type, D = NULL, method.h = NULL,
            h = h_opt, lambda = current_lambda, sparse = FALSE, gc = FALSE, chunk_size = chunk_size)

        # Store results
        estimators_array[, , b] = est_star$estimator
    }

    # 7. Clean up and return
    if (show_progress) {
        pb$close()
        cat("Bootstrap completed.\n")
    } else {
        cat("\nBootstrap completed.\n")
    }

    # Assumes listN is available
    return(estimators_array)
}


#' Performs bootstrap resampling for panel data vector field estimations.
#'
#' @param result An object from `estimate_panel_vf` containing estimation results and parameters.
#' @param B Number of bootstrap replicates (default: 500).
#'
#' @return A 3D array of bootstrapped field estimators with dimensions [nEval x 2 x B].
bootstrapPanelVF <- function(result, B = 500) {
    # 1. Extract parameters from the result object
    X <- result$X
    x_eval <- result$x
    nEval <- result$nEval
    FE <- result$FE
    TE <- result$TE
    uniform_weights <- result$uniform_weights
    kernel.type <- result$kernel.type
    method.h <- result$method.h
    chunk_size <- result$chunk_size
    sparse <- result$sparse
    gc <- result$gc
    
    X0_raw <- result$X0_raw
    Y1 <- result$Y1_unrolled
    Y2 <- result$Y2_unrolled

    # 2. Prepare data and storage
    dims <- dim(X)
    nObs <- dims[1]
    nT <- dims[3]
    nEvalPoints <- nrow(x_eval)
    estimators_array <- array(NA, dim = c(nEvalPoints, 2, B))
    FE_array <- array(NA, dim = c(nObs, 2, B))
    TE_array <- array(NA, dim = c(nT-1, 2, B))

    # 3. Setup progress bar
    show_progress <- !exists("DEBUG") || !DEBUG
    if (show_progress) {
        cat("Running bootstrap with", B, "replicates for panel VF...\n")
        pb <- createCustomProgressBar(min = 0, max = B)
    } else {
        cat("Starting bootstrap with", B, "replicates for panel VF...\n")
    }

    # 4. Bootstrap loop
    for (b in 1:B) {
        # Resample individuals (cross-sectional bootstrap)
        idx_star <- sample(1:nObs, nObs, replace = TRUE)
        X_star <- X[idx_star, , , drop = FALSE]

        if (show_progress) {
            pb$update(b, paste("Bootstrap iteration", b, "/", B))
        }

        # Re-estimate using the bootstrapped data
        # All original estimation parameters are passed through
        est_star <- estimate_panel_vf(
            X = X_star,
            x = x_eval,
            nEval = nEvalPoints,
            FE = FE,
            TE = TE,
            uniform_weights = uniform_weights,
            kernel.type = kernel.type,
            method.h = method.h,
            chunk_size = chunk_size,
            sparse = sparse,
            gc = gc
        )
        effects <- get_effects(est_star, X_obs = X0_raw, FE = TRUE, TE = TRUE)
        alpha_i_hat <- effects$alpha_i
        gamma_t_hat <- effects$gamma_t

        # Store results
        estimators_array[, , b] <- est_star$estimator
        FE_array[, , b] <- alpha_i_hat
        TE_array[, , b] <- gamma_t_hat
    }

    # 5. Clean up and return
    if (show_progress) {
        pb$close()
        cat("Bootstrap completed.\n")
    } else {
        cat("\nBootstrap completed.\n")
    }

    return(listN(estimators_array, FE_array, TE_array))
}

#' Performs significance testing based on bootstrap samples for vector field estimations
#' 
#' @param result An object containing vector field estimation results
#' @param bootstrapSamples Array of bootstrap samples from bootstrapKernelFieldErrors (nEval x 2 x B)
#' @param p_crit Critical p-value threshold for significance testing (default: 0.05)
#' 
#' @return A logical vector indicating significance for each evaluation point (length nEval),
#'         where TRUE means the point has a significant vector field direction
significanceBootstrap <- function(result, bootstrapSamples, p_crit=0.05){
   
    nEval = nrow(result)
    estConfInt = matrix(NA, nrow=nEval, ncol=2)
    
    # Identify where we can calculate directional vectors (non-zero and non-NaN values)
    whereEst = (result[,1]!=0 & result[,2]!=0) & (!is.nan(result[,1]) & !is.nan(result[,2]))
    whereEstInd = which(whereEst)
    
    # Calculate significance based on bootstrap distribution
    estConfInt[whereEstInd,1] = sapply(whereEstInd, function(i) max(mean(bootstrapSamples[i,1,]<0,na.rm=T),mean(bootstrapSamples[i,1,]>0,na.rm=T)) )
    estConfInt[whereEstInd,2] = sapply(whereEstInd, function(i) max(mean(bootstrapSamples[i,2,]<0,na.rm=T),mean(bootstrapSamples[i,2,]>0,na.rm=T)) )

    # If all bootstrap samples are NA, set to non-significant
    estConfInt[is.na(estConfInt)] = FALSE
    
    # Consider significant only if at least one component is significant
    signifEst = rowSums(estConfInt>(1-p_crit)) >= 1
    
    return(signifEst) 
}


