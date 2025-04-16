significanceVF <- function(est,X0,X1,alpha = 0.05){
    

    h = est$h
    
    n = nrow(X0)
    Y = X1-X0
    Yhat = interp2d(X0,est$x,est$estimator)
    eps = Y - Yhat    
    
    # homoskedastic errors
    Sigma = cov(eps, use = "complete.obs")

    if (est$kernel.type == "epa"){ k = 3/5 
    } else if (est$kernel == "gauss") { k = 1/(2*sqrt(pi)) }
    
    
    Var = array(data = NA, dim = c(2,2,nrow(x)))
    ChiSquare_stat = rep(NA, nrow(x))     
    p_values = rep(NA, nrow(x))
    signif = rep(NA, nrow(x))
    for (i in 1:nrow(x)){
        Var[,,i] = (k^2 * Sigma / (est$density[i] * n*h^2)) 
        ChiSquare_stat[i] = t(est$estimator[i,]) %*% solve(Var[,,i]) %*% est$estimator[i,]
        p_values[i] = 1-pchisq(ChiSquare_stat[i], 2)
        p_values[i] = ifelse( is.na(p_values[i]), 1 , p_values[i])
        signif[i] = ifelse(p_values[i] < alpha, 1, 0)
    }
    
    
    return(listN(Var,ChiSquare_stat,p_values,signif))
}

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






