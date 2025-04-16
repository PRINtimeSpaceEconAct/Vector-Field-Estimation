LLregression <- function(X, Y, x=NULL, nEval=2500, kernel.type="gauss", D=NULL, 
                        method.h=NULL, h=NULL, lambda = NULL, 
                        sparse=FALSE, gc=FALSE, chunk_size=nrow(x)) {

    type.est = "LL"
    resultEst = kernelMethod(X=X,x=x,nEval=nEval,kernel.type=kernel.type,D=D,
                             method.h=method.h,h=h,lambda=lambda,
                             sparse=sparse,gc=gc,chunk_size=chunk_size,type.est=type.est, Y=Y)
    
    return(resultEst)
}

# Adaptive bandwidth version of LL regression
LLregressionAdaptive <- function(X, Y, x=NULL, nEval=2500, kernel.type="gauss", D=NULL,
                                method.h=NULL, h=NULL, lambda=NULL,
                                sparse=FALSE, gc=FALSE, chunk_size=nrow(x), alpha = 0.5) {
    
    # Get adaptive bandwidths using the density-based approach
    lambda = getLocalBandwidth(X, kernel.type=kernel.type, D=D, method.h=method.h, h=h,
                             sparse=sparse, gc=gc, chunk_size=chunk_size, alpha = alpha)
    
    # Apply LL regression with adaptive bandwidths
    est = LLregression(X, Y, x=x, nEval=nEval, kernel.type=kernel.type, D=D, 
                      method.h=method.h, h=h, lambda=lambda,
                      sparse=sparse, gc=gc, chunk_size=chunk_size)
    est$lambda = lambda
    return(est)
}

# Helper function for core LL field component estimation
computeLLFieldComponents <- function(X0, Y1, Y2, x, nEval, kernel.type, D, 
                             method.h, h, lambda, sparse, gc, chunk_size) {
    
    est1 = LLregression(X0, Y1, x=x, nEval=nEval, kernel.type=kernel.type, D=D,
                      method.h=method.h, h=h, lambda=lambda,
                      sparse=sparse, gc=gc, chunk_size=chunk_size)
    
    # Use the evaluation points from the first estimate for the second
    x_eval = est1$x
    est2 = LLregression(X0, Y2, x=x_eval, nEval=nEval, kernel.type=kernel.type, D=D,
                      method.h=method.h, h=h, lambda=lambda,
                      sparse=sparse, gc=gc, chunk_size=chunk_size)
    
    estimator = cbind(est1$estimator, est2$estimator)
    density = est1$density
    h = est1$h
    method.h = est1$method.h
    kernel.type = est1$kernel.type
    lambda = est1$lambda
    Hkk_values = est1$Hkk_values # Get the full vector of Hkk values
    # Return relevant components needed later
    return(listN(estimator, x_eval, density, h, method.h, kernel.type, lambda, Hkk_values))
}


LLfield <- function(X0, X1, x=NULL, nEval=2500, kernel.type="gauss", D=NULL,
                   method.h=NULL, h=NULL, lambda=NULL,
                   sparse=FALSE, gc=FALSE, chunk_size=nrow(x),
                   hOpt = FALSE, nGridh = 10) {

    # Calculate differences for vector field
    Y1 = X1[,1] - X0[,1]
    Y2 = X1[,2] - X0[,2]
    
    # Initialize optimization variables
    optVars = initializeOptimizationVariables()

    # Setup parameters and validate inputs - validation now happens in this function
    optParams = setupOptimizationParameters(X0=X0, kernel.type=kernel.type,
                                           hOpt=hOpt, nGridh=nGridh, h=h, method.h=method.h,
                                           lambda=lambda, isAdaptive=FALSE)

    # Get validated/determined parameters
    h = optParams$h
    method.h = optParams$method.h
    optVars$best_h = h  # Initialize best_h with determined h

    if (hOpt) {
        optVars$hGrid = optParams$hGrid # Store grid if optimizing
        optVars$debugArrays = debugInitArrays(optParams$nGridh_actual) # Use actual grid size
        optVars$AICc_values = array(NA, dim = optParams$nGridh_actual)

        # Create progress bar if not in debug mode
        show_progress = !exists("DEBUG") || !DEBUG
        if (show_progress) {
            cat("Optimizing bandwidth parameter h...\n")
            pb = createCustomProgressBar(min = 0, max = optParams$nGridh_actual)
        }
        
        # --- Optimization Loop h ---
        for (i in 1:optParams$nGridh_actual) {
            debugPrint("Computing h %d/%d", i, optParams$nGridh_actual)
            hi = optParams$hGrid[i]

            # Update progress bar with current parameter value
            if (show_progress) {
                pb$update(i, formatOptimParams(h = hi))
            }

            # Use full chunk_size? Check logic
            est_components_i = computeLLFieldComponents(X0, Y1, Y2, x=X0, nEval=optParams$Nobs,
                                                        kernel.type=kernel.type, D=D,
                                                        method.h=method.h, h=hi, lambda=NULL,
                                                        sparse=sparse, gc=gc, chunk_size=chunk_size) 

            X1Hat_i = X0 + est_components_i$estimator

            metrics_i = calculateAICc(X0=X0, X1=X1, X1Hat=X1Hat_i, 
                                          lambda=NULL, Nobs=optParams$Nobs,
                                          kernelFunction=optParams$kernelFunction,
                                          method="LL", Hkk_values=est_components_i$Hkk_values)
            
            optVars$AICc_values[i] = metrics_i$AICc

            optVars$debugArrays = debugStoreValues(optVars$debugArrays, i, NULL, # Use NULL for 1D
                                                  AICc=metrics_i$AICc,
                                                  RSS=metrics_i$RSS,
                                                  trH=metrics_i$trH,
                                                  freedom=metrics_i$freedom)

            if (!is.na(metrics_i$AICc) && metrics_i$AICc < optVars$best_AICc) {
                optVars$best_AICc = metrics_i$AICc
                optVars$best_h = hi
            }
        }
        # --- End Optimization Loop h ---
        
        # Close progress bar if shown
        if (show_progress) { pb$close() }

        h = optVars$best_h # Update h with the best found value
        printOptimizationResults(hGrid = optVars$hGrid,
                                 h = h,
                                 arrays = optVars$debugArrays,
                                 alpha = NULL,
                                 alphaGrid = NULL)
    } # end if hOpt 

    # Final estimation using the chosen h (either provided, method-derived, or optimized)
    final_x = if (is.null(x)) X0 else x
    final_chunk_size = if (is.null(x)) nrow(X0) else chunk_size
    final_est_components = computeLLFieldComponents(X0, Y1, Y2, x=final_x, nEval=nEval,
                                                   kernel.type=kernel.type, D=D,
                                                   method.h=method.h, h=h, lambda=NULL,
                                                   sparse=sparse, gc=gc, chunk_size=final_chunk_size)
    # Construct result object
    result = list(x = final_est_components$x_eval,
                   X0 = X0, X1 = X1,
                   estimator = final_est_components$estimator,
                   type.est = "LL",
                   density = final_est_components$density,
                   Hkk_values = final_est_components$Hkk_values,
                   kernel.type = final_est_components$kernel.type,
                   h = final_est_components$h, # Use h from final components
                   method.h = final_est_components$method.h,
                   lambda = NULL)

    # Add optimization results if performed
    if (hOpt){
        result$AICc_values = optVars$AICc_values
        result$hGrid = optVars$hGrid
    }
    
    return(result)
}

LLfieldAdaptive <- function(X0, X1, x=NULL, nEval=2500, kernel.type="gauss", D=NULL,
                           method.h=NULL, h=NULL, lambda=NULL,
                           sparse=FALSE, gc=FALSE, chunk_size=nrow(x), alpha = 0.5,
                           hOpt = FALSE, nGridh = 10, alphaOpt = FALSE, nGridAlpha = 5) {
    
    # Calculate differences for vector field
    Y1 = X1[,1] - X0[,1]
    Y2 = X1[,2] - X0[,2]
    
    # Initialize optimization variables
    optVars = initializeOptimizationVariables()

    # Setup parameters and validate inputs - validation now happens in this function
    optParams = setupOptimizationParameters(X0=X0, kernel.type=kernel.type,
                                           hOpt=hOpt, nGridh=nGridh, h=h, method.h=method.h,
                                           alphaOpt=alphaOpt, nGridAlpha=nGridAlpha, alpha=alpha,
                                           lambda=lambda, isAdaptive=TRUE)

    # Get validated/determined parameters
    current_h = optParams$h
    current_method_h = optParams$method.h
    current_alpha = optParams$alpha
    
    # Initialize best values based on fixed values if not optimizing
    optVars$best_h = if(hOpt) NULL else current_h
    optVars$best_alpha = if(alphaOpt) NULL else current_alpha

    if (hOpt || alphaOpt) {
        # Store grids from setup
        optVars$hGrid = optParams$hGrid
        optVars$alphaGrid = optParams$alphaGrid

        # Initialize debug arrays and AICc storage based on optimization dimensions
        if (alphaOpt) {
            optVars$debugArrays = debugInitArrays(optParams$nGridh_actual, optParams$nGridAlpha_actual)
            optVars$AICc_values = array(NA, dim = c(optParams$nGridh_actual, optParams$nGridAlpha_actual))
        } else {
            optVars$debugArrays = debugInitArrays(optParams$nGridh_actual) # 1D
            optVars$AICc_values = array(NA, dim = optParams$nGridh_actual)
        }

        # Setup progress bar if not in debug mode
        show_progress = !exists("DEBUG") || !DEBUG
        total_iterations = optParams$nGridh_actual
        if (alphaOpt) total_iterations = total_iterations * optParams$nGridAlpha_actual
        
        if (show_progress) {
            if (hOpt && alphaOpt) { print("Optimizing bandwidth parameter h and alpha...")
            } else if (hOpt) { print("Optimizing bandwidth parameter h...")
            } else { print("Optimizing alpha parameter...") }
            pb = createCustomProgressBar(min = 0, max = total_iterations)
            progress_counter = 0
        }

        # --- Optimization Loop ---
        for (i in 1:optParams$nGridh_actual) {
            hi = optParams$hGrid[i] # Use h from the grid (single value if hOpt=F)
            debugPrint("Computing h %d/%d (h=%.4f)", i, optParams$nGridh_actual, hi)

            # Pilot density calculation (only depends on h)
            # Pass method.h=NULL since h is explicitly given by hi
            pilotDensity = densityEst2d(X0, x=X0, nEval=optParams$Nobs, kernel.type=kernel.type,
                                       D=D, method.h=NULL, h=hi, sparse=sparse,
                                       gc=gc, chunk_size=chunk_size)

            # Calculate geometric mean directly without any checks
            g = exp(mean(log(pilotDensity$estimator)))

            for (j in 1:optParams$nGridAlpha_actual) {
                alpha_j = optParams$alphaGrid[j] # Use alpha from the grid (single value if alphaOpt=F)
                if (alphaOpt) debugPrint("  Computing alpha %d/%d (alpha=%.4f)", j, optParams$nGridAlpha_actual, alpha_j)

                # Update progress bar with current parameter values
                if (show_progress) {
                    progress_counter = progress_counter + 1
                    if (alphaOpt) {
                        pb$update(progress_counter, formatOptimParams(h = hi, alpha = alpha_j))
                    } else {
                        pb$update(progress_counter, formatOptimParams(h = hi))
                    }
                }

                # Calculate lambda directly
                lambda_ij = (pilotDensity$estimator / g)^(-alpha_j)

                # Pass method.h=NULL here, h=hi is explicit
                est_components_ij = computeLLFieldComponents(X0, Y1, Y2, x=X0, nEval=optParams$Nobs,
                                                             kernel.type=kernel.type, D=D,
                                                             method.h=NULL, h=hi, lambda=lambda_ij,
                                                             sparse=sparse, gc=gc, chunk_size=chunk_size)
                    
                X1Hat_ij = X0 + est_components_ij$estimator
                metrics_ij = calculateAICc(X0=X0, X1=X1, X1Hat=X1Hat_ij, 
                                                lambda=lambda_ij, Nobs=optParams$Nobs,
                                                kernelFunction=optParams$kernelFunction,
                                                method="LL", Hkk_values=est_components_ij$Hkk_values)

                # Store results and update best
                aic_index_1 = i
                aic_index_2 = if(alphaOpt) j else NULL

                if (alphaOpt) {
                   optVars$AICc_values[aic_index_1, aic_index_2] = metrics_ij$AICc
                } else {
                   optVars$AICc_values[aic_index_1] = metrics_ij$AICc # Index directly for 1D
                }

                 optVars$debugArrays = debugStoreValues(optVars$debugArrays, aic_index_1, aic_index_2,
                                                       AICc=metrics_ij$AICc,
                                                       RSS=metrics_ij$RSS,
                                                       trH=metrics_ij$trH,
                                                       freedom=metrics_ij$freedom)

                if (!is.na(metrics_ij$AICc) && metrics_ij$AICc < optVars$best_AICc) {
                    optVars$best_AICc = metrics_ij$AICc
                    optVars$best_h = hi
                    optVars$best_alpha = alpha_j
                    optVars$best_lambda = lambda_ij # Store the lambda corresponding to best h/alpha
                }
            } # End alpha loop
        } # End h loop
        # --- End Optimization Loop ---
        
        # Close progress bar if shown
        if (show_progress) { pb$close() }

        # Set the final parameters based on optimization results
        current_h = optVars$best_h
        current_alpha = optVars$best_alpha
        final_lambda = optVars$best_lambda # Use the lambda found for the best parameters

        
        # Print optimization results
        printOptimizationResults(hGrid = if(hOpt) optVars$hGrid else NULL,
                                 h = current_h,
                                 arrays = optVars$debugArrays,
                                 alpha = current_alpha,
                                 alphaGrid = if(alphaOpt) optVars$alphaGrid else NULL)
        
    } else {
        # No optimization: Calculate lambda directly from provided/determined h and alpha
        # current_h and current_alpha were set after calling setupOptimizationParameters
        # Pass method.h=NULL as h is explicitly known
        final_lambda = getLocalBandwidth(X0, kernel.type=kernel.type, D=D, method.h=NULL, h=current_h,
                            sparse=sparse, gc=gc, chunk_size=chunk_size, alpha = current_alpha)
    }

    # Final estimation using the chosen/optimized parameters
    final_x = if (is.null(x)) X0 else x
    final_chunk_size = if (is.null(x)) nrow(X0) else chunk_size

    # Pass method.h=NULL as h is fixed
    final_est_components = computeLLFieldComponents(X0, Y1, Y2, x=final_x, nEval=nEval,
                                                   kernel.type=kernel.type, D=D,
                                                   method.h=NULL, h=current_h, lambda=final_lambda,
                                                   sparse=sparse, gc=gc, chunk_size=final_chunk_size)
    # Construct result object
    result = list(x = final_est_components$x_eval,
                   X0 = X0, X1 = X1,
                   estimator = final_est_components$estimator,
                   type.est = "LL_adaptive",
                   density = final_est_components$density,
                   Hkk_values = final_est_components$Hkk_values,
                   kernel.type = final_est_components$kernel.type,
                   h = final_est_components$h, # Global h used
                   method.h = current_method_h, # Report method used for h
                   lambda = final_lambda, # Final lambda used (could be vector)
                   alpha = current_alpha) # Alpha used for the final lambda


    # Add optimization results if any optimization was performed
    if (hOpt || alphaOpt){
        result$AICc_values = optVars$AICc_values
        if(hOpt) result$hGrid = optVars$hGrid
        # result$alpha is already added above
        if(alphaOpt) result$alphaGrid = optVars$alphaGrid
    }

    return(result)
}