#' Performs the computationally intensive steps of the panel vector field analysis.
#'
#' This function runs the estimation and bootstrap procedures for panel vector
#' field analysis. By default, it displays plots of the results during execution.
#' This allows immediate visualization of the estimated vector field and effects.
#' The function also returns an object containing all results for further analysis
#' or customized plotting using `plotPanelVFAnalysis`.
#'
#' @param panel_df A data frame containing the panel data.
#' @param var_cols A character vector of length 2 specifying the variable columns (e.g., c("log_GDP", "log_LE")).
#' @param id_col A character string specifying the country/region identifier column.
#' @param time_col A character string specifying the time column.
#' @param nEval Integer, total number of evaluation grid points (default 2500).
#' @param FE Logical, include fixed effects in within transform (default TRUE).
#' @param TE Logical, include time effects in within transform (default TRUE).
#' @param uniform_weights Logical, use uniform weights for within transform (default TRUE).
#' @param kernel.type Kernel type for estimation (default "epa").
#' @param method.h Bandwidth rule (default "silverman").
#' @param chunk_size Integer chunk size for compute functions (default 512).
#' @param bootstrap_B Integer, number of bootstrap replicates (default 100). If NULL or <= 1, bootstrap is skipped.
#' @param p_crit Numeric, significance level for directional significance (default 0.05).
#' @param show_plots Logical, display plots during execution (default TRUE).
#' @param component_names A character vector of length 2 with names for the x and y components (default c("X1", "X2")).
#' @param lengthArrows Numeric, scaling factor for arrow lengths in vector field plot (default 1).
#' @param rescale Logical, rescale axes and field relative to a reference time slice (default TRUE).
#' @param rescale_ref_index Integer in [1, nT], reference time slice for rescaling (default is last time slice).
#' @param years Optional numeric/integer vector of length nT to label the TE plot's x-axis (default 1:nT).
#' @param label_names Optional character vector of length nObs for FE scatter labels (default NULL -> no labels).
#'
#' @return A list containing all inputs, results, and intermediate objects.
#' @export
runPanelVFAnalysis <- function(panel_df,
                               var_cols,
                               id_col,
                               time_col,
                               nEval = 2500,
                               FE = TRUE,
                               TE = TRUE,
                               uniform_weights = TRUE,
                               kernel.type = "epa",
                               method.h = "silverman",
                               chunk_size = 512,
                               bootstrap_B = 100,
                               p_crit = 0.05,
                               estimation_method = "LL",
                               adaptive = FALSE,
                               hOpt = FALSE,
                               alphaOpt = FALSE,
                               h = NULL,
                               alpha = 0.5,
                               ridge_param = 1e-5,
                               rcond_threshold = 1e-3,
                               show_plots = TRUE,
                               component_names = c("X1", "X2"),
                               lengthArrows = 1,
                               rescale = TRUE,
                               rescale_ref_index = NULL,
                               years = NULL,
                               label_names = NULL) {

    is_two_df_case <- is.list(panel_df) && length(panel_df) == 2 && all(sapply(panel_df, is.data.frame))

    if (is_two_df_case) {
        if (FE || TE) {
            warning("Two dataframes provided. Forcing FE=FALSE and TE=FALSE for this mode.")
        }
        FE <- FALSE
        TE <- FALSE

        two_df_data <- handle_two_df_input(panel_df, var_cols)
        X0 <- two_df_data$X0
        X1 <- two_df_data$X1
        X <- two_df_data$X
        dims <- two_df_data$dims
        nObs <- two_df_data$nObs
        nT <- two_df_data$nT
        X_unrolled <- two_df_data$X_unrolled

    } else {
        # Convert panel dataframe to 3D array
        X <- panel_df_to_array(panel_df = panel_df, 
                                 var_cols = var_cols, 
                                 id_col = id_col, 
                                 time_col = time_col)
        # --- Basic checks ---
        dims <- dim(X)
        nObs <- dims[1]
        nT <- dims[3]
        X_unrolled <- aperm(X, c(1, 3, 2))
        dim(X_unrolled) <- c(nObs * nT, 2)
    }


    cat("Panel vector field analysis started.\n")

    # Initialize variables for optimization grids (used in summary)
    h_grid_opt <- NULL
    alpha_grid_opt <- NULL

    # 1) Evaluation grid
    x <- defineEvalPoints(X_unrolled, nEval)

    # --- Estimation ---
    if (FE || TE) {
        # 2) Panel Estimation
        if (adaptive) {
            # If hOpt or alphaOpt is TRUE, first run LLfieldAdaptive to find optimal parameters
            if (hOpt || alphaOpt) {
                cat("Running LLfieldAdaptive to find optimal h and/or alpha...\n")
                
                # Prepare data for LLfieldAdaptive (same as non-panel case)
                X0 <- aperm(X[, , 1:(nT - 1), drop=FALSE], c(1, 3, 2))
                dim(X0) <- c(nObs * (nT - 1), 2)
                X1 <- aperm(X[, , 2:nT, drop=FALSE], c(1, 3, 2))
                dim(X1) <- c(nObs * (nT - 1), 2)
                
                # Run LLfieldAdaptive to optimize h and/or alpha
                opt_results <- LLfieldAdaptive(
                    X0, X1, 
                    x = x, 
                    nEval = nEval, 
                    kernel.type = kernel.type, 
                    method.h = method.h, 
                    h = h, 
                    chunk_size = chunk_size, 
                    hOpt = hOpt, 
                    alphaOpt = alphaOpt, 
                    alpha = alpha
                )
                
                # Extract optimal parameters
                h_opt <- opt_results$h
                alpha_opt <- opt_results$alpha
                
                # Extract grids for later use in summary
                h_grid_opt <- if(hOpt) opt_results$hGrid else NULL
                alpha_grid_opt <- if(alphaOpt) opt_results$alphaGrid else NULL
                
                cat("Optimal h:", round(h_opt, 4), "\n")
                cat("Optimal alpha:", round(alpha_opt, 4), "\n")
                cat("Running estimate_panel_vf_adaptive with optimal parameters...\n")
                
                # Now run the panel estimation with optimal parameters
                panel_vf_results <- estimate_panel_vf_adaptive(
                    X = X,
                    x = x,
                    nEval = nEval,
                    FE = FE,
                    TE = TE,
                    uniform_weights = uniform_weights,
                    kernel.type = kernel.type,
                    method.h = method.h,
                    h = h_opt,
                    chunk_size = chunk_size,
                    gc = FALSE,
                    alpha = alpha_opt,
                    ridge_param = ridge_param,
                    rcond_threshold = rcond_threshold
                )
            } else {
                # Standard adaptive estimation without optimization
                panel_vf_results <- estimate_panel_vf_adaptive(
                    X = X,
                    x = x,
                    nEval = nEval,
                    FE = FE,
                    TE = TE,
                    uniform_weights = uniform_weights,
                    kernel.type = kernel.type,
                    method.h = method.h,
                    h = h,
                    chunk_size = chunk_size,
                    gc = FALSE,
                    alpha = alpha,
                    ridge_param = ridge_param,
                    rcond_threshold = rcond_threshold
                )
            }
        } else {
            panel_vf_results <- estimate_panel_vf(
                X = X,
                x = x,
                nEval = nEval,
                FE = FE,
                TE = TE,
                uniform_weights = uniform_weights,
                kernel.type = kernel.type,
                method.h = method.h,
                h = h,
                chunk_size = chunk_size,
                gc = FALSE,
                ridge_param = ridge_param,
                rcond_threshold = rcond_threshold
            )
        }
        VF_hat <- panel_vf_results$estimator
        X0_raw <- panel_vf_results$X0_raw

    } else {
        # 2) Standard VF Estimation (no FE/TE)
        
        # Prepare data for standard VF functions
        if (!is_two_df_case) {
            X0 <- aperm(X[, , 1:(nT - 1), drop=FALSE], c(1, 3, 2))
            dim(X0) <- c(nObs * (nT - 1), 2)
            X1 <- aperm(X[, , 2:nT, drop=FALSE], c(1, 3, 2))
            dim(X1) <- c(nObs * (nT - 1), 2)
        }

        if (estimation_method == "LL") {
            if (adaptive) {
                vf_results <- LLfieldAdaptive(X0, X1, x = x, nEval = nEval, kernel.type = kernel.type, method.h = method.h, h = h, chunk_size = chunk_size, hOpt = hOpt, alphaOpt = alphaOpt, alpha = alpha)
            } else {
                vf_results <- LLfield(X0, X1, x = x, nEval = nEval, kernel.type = kernel.type, method.h = method.h, h = h, chunk_size = chunk_size, hOpt = hOpt)
            }
        } else if (estimation_method == "NW") {
            if (adaptive) {
                vf_results <- NWfieldAdaptive(X0, X1, x = x, nEval = nEval, kernel.type = kernel.type, method.h = method.h, h = h, chunk_size = chunk_size, hOpt = hOpt, alphaOpt = alphaOpt, alpha = alpha)
            } else {
                vf_results <- NWfield(X0, X1, x = x, nEval = nEval, kernel.type = kernel.type, method.h = method.h, h = h, chunk_size = chunk_size, hOpt = hOpt)
            }
        } else {
            stop("Invalid estimation_method. Choose 'LL' or 'NW'.")
        }
        
        # Extract optimization grids for summary (if available)
        h_grid_opt <- if(hOpt && !is.null(vf_results$hGrid)) vf_results$hGrid else NULL
        alpha_grid_opt <- if(alphaOpt && !is.null(vf_results$alphaGrid)) vf_results$alphaGrid else NULL
        
        panel_vf_results <- vf_results # To maintain compatibility with bootstrap
        VF_hat <- vf_results$estimator
        X0_raw <- vf_results$X0
    }

    VF_derivatives_hat <- NULL
    if (FE || TE) {
        VF_derivatives_hat <- list(
            d1 = panel_vf_results$derivative_estimator_1$estimator,
            d2 = panel_vf_results$derivative_estimator_2$estimator
        )
    }

    # 3) Bootstrap
    run_bootstrap <- !is.null(bootstrap_B) && bootstrap_B > 1

    if (run_bootstrap) {
        if (FE || TE) {
            bootstrap_samples <- bootstrapPanelVF(panel_vf_results, B = bootstrap_B)
        } else {
            # For non-panel models, we need a different bootstrap function
            # and to structure its output to match what the rest of the script expects.
            estimators_array <- bootstrapKernelFieldErrors(panel_vf_results, B = bootstrap_B)
            bootstrap_samples <- list(
                estimators_array = estimators_array,
                FE_array = NULL,
                TE_array = NULL
            )
        }
    } else {
        # Initialize with NULL to prevent errors in later steps
        bootstrap_samples <- list(estimators_array = NULL, FE_array = NULL, TE_array = NULL)
    }

    # 4) Effects and significance
    
    if (FE || TE) {
        effects <- get_effects(panel_vf_results, X_obs = X0_raw, FE = TRUE, TE = TRUE)
        FE_hat <- effects$alpha_i
        TE_hat <- effects$gamma_t
        if (run_bootstrap) {
            FE_signif <- significanceBootstrap(FE_hat, bootstrap_samples$FE_array, p_crit = p_crit)
            TE_signif <- significanceBootstrap(TE_hat, bootstrap_samples$TE_array, p_crit = p_crit)
        } else {
            FE_signif <- rep(TRUE, nrow(FE_hat))
            TE_signif <- rep(TRUE, nrow(TE_hat))
        }
    } else {
        # No FE/TE to compute or test
        FE_hat <- NULL
        TE_hat <- NULL
        FE_signif <- NULL
        TE_signif <- NULL
    }

    if (run_bootstrap) {
        VF_signif <- significanceBootstrap(VF_hat, bootstrap_samples$estimators_array, p_crit = p_crit)
    } else {
        VF_signif <- rep(TRUE, nrow(VF_hat))
    }
    
    cat("Analysis completed.\n")

    # --- Print summary ---
    cat("\n--- Analysis Summary ---\n")
    cat("Estimation Method:", estimation_method, if(adaptive) "(Adaptive)" else "", "\n")
    cat("Kernel Type:", kernel.type, "\n")
    
    if (hOpt) {
        if (!is.null(h_grid_opt)) {
            h_range <- range(h_grid_opt, na.rm = TRUE)
            cat("Bandwidth (h) optimized over [", round(h_range[1], 4), ", ", round(h_range[2], 4), "]\n", sep = "")
            cat("  - Optimal h found:", round(panel_vf_results$h, 4), "\n")
        } else {
            cat("Bandwidth (h) optimized\n")
            cat("  - Optimal h found:", round(panel_vf_results$h, 4), "\n")
        }
    } else if (is.null(h)) {
        cat("Bandwidth (h) Method:", method.h, "\n")
        cat("  - Estimated h:", round(panel_vf_results$h, 4), "\n")
    } else {
        cat("Bandwidth (h) Provided:", round(h, 4), "\n")
    }

    if (adaptive) {
        if (alphaOpt) {
            if (!is.null(alpha_grid_opt)) {
                alpha_range <- range(alpha_grid_opt, na.rm = TRUE)
                cat("Alpha optimized over [", round(alpha_range[1], 4), ", ", round(alpha_range[2], 4), "]\n", sep = "")
                cat("  - Optimal alpha found:", round(panel_vf_results$alpha, 4), "\n")
            } else {
                cat("Alpha optimized\n")
                cat("  - Optimal alpha found:", round(panel_vf_results$alpha, 4), "\n")
            }
        } else {
            cat("Alpha:", round(panel_vf_results$alpha, 4), "\n")
        }
    }
    
    cat("Fixed Effects (FE):", FE, "\n")
    cat("Time Effects (TE):", TE, "\n")
    cat("Bootstrap Replicates:", if(run_bootstrap) bootstrap_B else 0, "\n")
    cat("Chunk Size:", chunk_size, "\n")
    cat("------------------------\n\n")

    # --- Plotting ---
    if (show_plots) {
        cat("Generating plots...\n")
        
        # Validate component_names
        if (length(component_names) != 2 || !is.character(component_names)) {
            warning("`component_names` must be a character vector of length 2. Using default names 'X1', 'X2'.")
            component_names <- c("X1", "X2")
        }
        
        # Set up years if needed for TE plots
        if (is.null(years)) {
            years <- 1:nT
        } else if (length(years) != nT) {
            warning("`years` length does not match time dimension of X. Using 1:nT instead.")
            years <- 1:nT
        }
        
        # Prepare plot data (filter by significance and optionally rescale)
        VF_hat_plot <- VF_hat
        VF_hat_plot[!VF_signif, ] <- 0
        x_plot <- x
        X0_raw_plot <- X0_raw
        
        if (rescale) {
            if (is.null(rescale_ref_index)) rescale_ref_index <- nT
            if (rescale_ref_index < 1 || rescale_ref_index > nT) {
                warning("rescale_ref_index out of bounds. Skipping rescale.")
            } else {
                ref1 <- stats::median(X[, 1, rescale_ref_index], na.rm = TRUE)
                ref2  <- stats::median(X[, 2, rescale_ref_index], na.rm = TRUE)
                if (ref1 == 0 || ref2 == 0 || is.na(ref1) || is.na(ref2)) {
                    warning("Rescale reference medians invalid; skipping rescale.")
                } else {
                    x_plot[, 1] <- x_plot[, 1] / ref1
                    x_plot[, 2] <- x_plot[, 2] / ref2
                    VF_hat_plot[, 1] <- VF_hat_plot[, 1] / ref1
                    VF_hat_plot[, 2] <- VF_hat_plot[, 2] / ref2
                    X0_raw_plot[, 1] <- X0_raw_plot[, 1] / ref1
                    X0_raw_plot[, 2] <- X0_raw_plot[, 2] / ref2
                }
            }
        }
        
        old_par <- graphics::par(no.readonly = TRUE)
        on.exit(graphics::par(old_par), add = TRUE)
        
        # 1) Vector Field Plot
        graphics::plot(x_plot, type = "n",
                      xlab = component_names[1],
                      ylab = component_names[2],
                      main = "Estimated Vector Field")
        if(rescale) {
            graphics::abline(h = 1, lty = 3)
            graphics::abline(v = 1, lty = 3)
        } else {
            graphics::abline(h = 0, lty = 3)
            graphics::abline(v = 0, lty = 3)
        }
        graphics::arrows(x_plot[, 1], x_plot[, 2],
                        x_plot[, 1] + lengthArrows * VF_hat_plot[, 1],
                        x_plot[, 2] + lengthArrows * VF_hat_plot[, 2],
                        angle = 15, col = "black", length = 0.05)
        
        # Add density contours if sm package is available
        if (requireNamespace("sm", quietly = TRUE)) {
            X0_raw_plot_complete <- X0_raw_plot[complete.cases(X0_raw_plot), ]
            if(nrow(X0_raw_plot_complete) > 1) {
                est.dens <- sm::sm.density(X0_raw_plot_complete, display = "none")
                graphics::contour(est.dens$eval.points[, 1], est.dens$eval.points[, 2], est.dens$estimate,
                                 add = TRUE, col = "purple")
            }
        }
        
        # 2) Fixed Effects Plots
        if (!is.null(FE_hat)) {
            cex_vals_fe <- ifelse(FE_signif, 1.0, 0.1)
            x_init1 <- X[, 1, 1]
            x_init2 <- X[, 2, 1]
            
            readline(prompt="Press [enter] to see the next plot")
            graphics::plot(x_init1, FE_hat[, 1], 
                          main = paste("Fixed Effects vs. Initial", component_names[1]), 
                          xlab = paste(component_names[1], "at t0"), 
                          ylab = paste("FE (", component_names[1], ")", sep = ""), 
                          pch = 19, cex = cex_vals_fe)
            graphics::abline(h = 0)
            graphics::grid()
            if (!is.null(label_names) && length(label_names) == nObs) {
                graphics::text(x_init1, FE_hat[, 1], labels = label_names, cex = 0.5, pos = 4, col = "red")
            }
            
            readline(prompt="Press [enter] to see the next plot")
            graphics::plot(x_init2, FE_hat[, 2], 
                          main = paste("Fixed Effects vs. Initial", component_names[2]), 
                          xlab = paste(component_names[2], "at t0"), 
                          ylab = paste("FE (", component_names[2], ")", sep = ""), 
                          pch = 19, cex = cex_vals_fe)
            graphics::abline(h = 0)
            graphics::grid()
            if (!is.null(label_names) && length(label_names) == nObs) {
                graphics::text(x_init2, FE_hat[, 2], labels = label_names, cex = 0.5, pos = 4, col = "red")
            }
        }
        
        # 3) Time Effects Plots with Bootstrap Bands (if bootstrap was run)
        if (!is.null(TE_hat)) {
            cex_vals_te <- ifelse(TE_signif, 1.0, 0.1)
            
            if (run_bootstrap) {
                TE_quantiles <- apply(bootstrap_samples$TE_array, c(1, 2), stats::quantile, 
                                     probs = c(0.025, 0.975), na.rm = TRUE)
                lower_bound <- TE_quantiles[1, , ]
                upper_bound <- TE_quantiles[2, , ]
                
                # Component 1
                readline(prompt="Press [enter] to see the next plot")
                y1_range <- range(c(lower_bound[, 1], upper_bound[, 1], TE_hat[, 1]), na.rm = TRUE)
                graphics::plot(years[2:nT], TE_hat[, 1], 
                              main = paste("Time Effects (", component_names[1], ")", sep = ""), 
                              xlab = "time", 
                              ylab = paste("TE (", component_names[1], ")", sep = ""), 
                              type = "n", ylim = y1_range)
                graphics::polygon(c(years[2:nT], rev(years[2:nT])), 
                                 c(lower_bound[, 1], rev(upper_bound[, 1])), 
                                 col = "grey80", border = NA)
                graphics::grid()
                graphics::abline(h = 0)
                graphics::lines(years[2:nT], TE_hat[, 1], pch = 19, cex = cex_vals_te, type = "b")
                
                # Component 2
                readline(prompt="Press [enter] to see the next plot")
                y2_range <- range(c(lower_bound[, 2], upper_bound[, 2], TE_hat[, 2]), na.rm = TRUE)
                graphics::plot(years[2:nT], TE_hat[, 2], 
                              main = paste("Time Effects (", component_names[2], ")", sep = ""), 
                              xlab = "time", 
                              ylab = paste("TE (", component_names[2], ")", sep = ""), 
                              type = "n", ylim = y2_range)
                graphics::polygon(c(years[2:nT], rev(years[2:nT])), 
                                 c(lower_bound[, 2], rev(upper_bound[, 2])), 
                                 col = "grey80", border = NA)
                graphics::grid()
                graphics::abline(h = 0)
                graphics::lines(years[2:nT], TE_hat[, 2], pch = 19, cex = cex_vals_te, type = "b")
            } else {
                # Without bootstrap, just simple line plots
                readline(prompt="Press [enter] to see the next plot")
                graphics::plot(years[2:nT], TE_hat[, 1], 
                              main = paste("Time Effects (", component_names[1], ")", sep = ""), 
                              xlab = "time", 
                              ylab = paste("TE (", component_names[1], ")", sep = ""), 
                              pch = 19, cex = cex_vals_te, type = "b")
                graphics::abline(h = 0)
                graphics::grid()
                
                readline(prompt="Press [enter] to see the next plot")
                graphics::plot(years[2:nT], TE_hat[, 2], 
                              main = paste("Time Effects (", component_names[2], ")", sep = ""), 
                              xlab = "time", 
                              ylab = paste("TE (", component_names[2], ")", sep = ""), 
                              pch = 19, cex = cex_vals_te, type = "b")
                graphics::abline(h = 0)
                graphics::grid()
            }
        }
        
        cat("Plotting finished.\n")
    }

    return(listN(
        # Inputs
        panel_df, var_cols, id_col, time_col, X,
        nEval, FE, TE, uniform_weights, kernel.type, method.h, chunk_size,
        bootstrap_B, p_crit, estimation_method, adaptive, hOpt, alphaOpt, h, alpha,
        ridge_param, rcond_threshold,
        # Results
        x, panel_vf_results, bootstrap_samples,
        VF_hat, FE_hat, TE_hat,
        VF_derivatives_hat,
        VF_signif, FE_signif, TE_signif
    ))
}

#' Generates plots from a completed panel vector field analysis.
#'
#' This function takes the results from `runPanelVFAnalysis` and produces
#' visualizations for the vector field, fixed effects, and time effects.
#' It allows for customization of plotting parameters without re-running the
#' expensive computations.
#'
#' @param analysis_results A list object returned by `runPanelVFAnalysis`.
#' @param timeInterval Numeric, step between consecutive panels (used to scale arrows, default 1).
#' @param rescale Logical, rescale axes and field relative to a reference time slice (default TRUE).
#' @param rescale_ref_index Integer in [1, nT], reference time slice for rescaling (default is last time slice).
#' @param years Optional numeric/integer vector of length nT to label the TE plot's x-axis (default 1:nT).
#' @param label_names Optional character vector of length nObs for FE scatter labels (default NULL -> no labels).
#' @param component_names A character vector of length 2 with names for the x and y components.
#' @param save_plots Logical, save the vector field plot to a file (default FALSE).
#' @param save_path Character, path where to save the plot (default "test_pics_VF").
#'
#' @return Invisibly returns NULL. This function is called for its side effect of generating plots.
#' @export
plotPanelVFAnalysis <- function(analysis_results,
                                timeInterval = 1,
                                rescale = TRUE,
                                rescale_ref_index = NULL,
                                years = NULL,
                                label_names = NULL,
                                component_names = c("X1", "X2"),
                                save_plots = FALSE,
                                save_path = "test_pics_VF",
                                lengthArrows = 1) {

    # --- Extract data from results object ---
    X <- analysis_results$X
    x <- analysis_results$x
    VF_hat <- analysis_results$VF_hat
    
    # Handle different result structures for X0_raw
    if (!is.null(analysis_results$panel_vf_results$X0_raw)) {
        X0_raw <- analysis_results$panel_vf_results$X0_raw
    } else if (!is.null(analysis_results$panel_vf_results$X0)) {
        X0_raw <- analysis_results$panel_vf_results$X0
    } else {
        stop("Could not find initial data points (X0_raw or X0) in analysis results.")
    }

    VF_signif <- analysis_results$VF_signif
    
    if (is.null(analysis_results$FE_hat)) {
        FE_hat <- NULL
        FE_signif <- NULL
        FE_bootstrap <- NULL
    } else {
        FE_hat <- analysis_results$FE_hat
        FE_signif <- analysis_results$FE_signif
        FE_bootstrap <- analysis_results$bootstrap_samples$FE_array
    }
    
    if (is.null(analysis_results$TE_hat)) {
        TE_hat <- NULL
        TE_signif <- NULL
        TE_bootstrap <- NULL
    } else {
        TE_hat <- analysis_results$TE_hat
        TE_signif <- analysis_results$TE_signif
        TE_bootstrap <- analysis_results$bootstrap_samples$TE_array
    }

    # --- Basic checks ---
    if (length(component_names) != 2 || !is.character(component_names)) {
        warning("`component_names` must be a character vector of length 2. Using default names 'X1', 'X2'.")
        component_names <- c("X1", "X2")
    }
    
    dims <- dim(X)
    nObs <- dims[1]
    nT <- dims[3]

    if (is.null(years)) {
        years <- 1:nT
    } else if (length(years) != nT) {
        warning("`years` length does not match time dimension of X. Using 1:nT instead.")
        years <- 1:nT
    }

    cat("Preparing plot data...\n")

    # --- Prepare plot data (filter/rescale) ---
    VF_hat_plot <- VF_hat
    VF_hat_plot[!VF_signif, ] <- 0
    x_plot <- x
    X0_raw_plot <- X0_raw

    if (rescale) {
        if (is.null(rescale_ref_index)) rescale_ref_index <- nT
        if (rescale_ref_index < 1 || rescale_ref_index > nT) stop("rescale_ref_index out of bounds.")
        ref1 <- stats::median(X[, 1, rescale_ref_index], na.rm = TRUE)
        ref2  <- stats::median(X[, 2, rescale_ref_index], na.rm = TRUE)
        if (ref1 == 0 || ref2 == 0 || is.na(ref1) || is.na(ref2)) {
            warning("Rescale reference medians invalid; skipping rescale.")
        } else {
            x_plot[, 1] <- x_plot[, 1] / ref1
            x_plot[, 2] <- x_plot[, 2] / ref2
            VF_hat_plot[, 1] <- VF_hat_plot[, 1] / ref1
            VF_hat_plot[, 2] <- VF_hat_plot[, 2] / ref2
            X0_raw_plot[, 1] <- X0_raw_plot[, 1] / ref1
            X0_raw_plot[, 2] <- X0_raw_plot[, 2] / ref2
        }
    }

    cat("Generating plots...\n")

    # Local function to generate the vector field plot
    plot_vf <- function() {
        graphics::plot(x_plot, type = "n",
                       xlab = component_names[1],
                       ylab = component_names[2],
                       main = "Estimated Vector Field")
        if(rescale) {
          graphics::abline(h = 1, lty = 3)
          graphics::abline(v = 1, lty = 3)
        } else {
          graphics::abline(h = 0, lty = 3)
          graphics::abline(v = 0, lty = 3)
        }
        graphics::arrows(x_plot[, 1], x_plot[, 2],
                         x_plot[, 1] + lengthArrows * VF_hat_plot[, 1],
                         x_plot[, 2] + lengthArrows * VF_hat_plot[, 2],
                         angle = 15, col = "black", length = 0.05)
        #graphics::points(X0_raw_plot[, 1], X0_raw_plot[, 2], pch = 19, col = "red", cex = 0.5)
        if (requireNamespace("sm", quietly = TRUE)) {
            X0_raw_plot_complete <- X0_raw_plot[complete.cases(X0_raw_plot), ]
            if(nrow(X0_raw_plot_complete) > 1) {
                est.dens <- sm::sm.density(X0_raw_plot_complete, display = "none")
                graphics::contour(est.dens$eval.points[, 1], est.dens$eval.points[, 2], est.dens$estimate,
                                  add = TRUE, col = "purple")
            }
        }
    }

    # --- Plotting ---
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par))
    
    if (save_plots) {
        # Save Vector Field plot and exit
        dir.create(save_path, showWarnings = FALSE)
        
        kernel_type_str <- tools::toTitleCase(analysis_results$kernel.type)
        
        ridge_val <- analysis_results$ridge_param
        ridge_str <- if (ridge_val == 0) "0" else gsub("e-0", "e-", formatC(ridge_val, format = "e", digits = 0))
        
        rcond_val <- analysis_results$rcond_threshold
        rcond_str <- if (rcond_val == 0) "0" else gsub("e-0", "e-", formatC(rcond_val, format = "e", digits = 0))
        
        filename <- sprintf("%s_ridge_%s_cond%s.svg", kernel_type_str, ridge_str, rcond_str)
        filepath <- file.path(save_path, filename)
        
        cat("Saving vector field plot to:", filepath, "\n")
        
        svg(filepath)
        plot_vf()
        dev.off()

    } else {
        # Interactive plotting session

        # Vector Field plot
        plot_vf()
        
        # FE plots
        if (!is.null(FE_hat)) {
            cex_vals_fe <- ifelse(FE_signif, 1.0, 0.1)
            x_init1 <- X[, 1, 1]
            x_init2 <- X[, 2, 1]

            readline(prompt="Press [enter] to see the next plot")
            graphics::plot(x_init1, FE_hat[, 1], main = paste("Fixed Effects vs. Initial", component_names[1]), xlab = paste(component_names[1], "at t0"), ylab = paste("FE (", component_names[1], ")", sep = ""), pch = 19, cex = cex_vals_fe)
            graphics::abline(h = 0)
            graphics::grid()
            if (!is.null(label_names) && length(label_names) == nObs) {
                graphics::text(x_init1, FE_hat[, 1], labels = label_names, cex = 0.5, pos = 4, col = "red")
            }
            
            readline(prompt="Press [enter] to see the next plot")
            graphics::plot(x_init2, FE_hat[, 2], main = paste("Fixed Effects vs. Initial", component_names[2]), xlab = paste(component_names[2], "at t0"), ylab = paste("FE (", component_names[2], ")", sep = ""), pch = 19, cex = cex_vals_fe)
            graphics::abline(h = 0)
            graphics::grid()
            if (!is.null(label_names) && length(label_names) == nObs) {
                graphics::text(x_init2, FE_hat[, 2], labels = label_names, cex = 0.5, pos = 4, col = "red")
            }
            
        }

        # TE plots with bootstrap bands
        if (!is.null(TE_hat)) {
            cex_vals_te <- ifelse(TE_signif, 1.0, 0.1)
            TE_quantiles <- apply(TE_bootstrap, c(1, 2), stats::quantile, probs = c(0.025, 0.975), na.rm = TRUE)
            lower_bound <- TE_quantiles[1, , ]
            upper_bound <- TE_quantiles[2, , ]

            # Component 1
            readline(prompt="Press [enter] to see the next plot")
            y1_range <- range(c(lower_bound[, 1], upper_bound[, 1], TE_hat[, 1]), na.rm = TRUE)
            graphics::plot(years[2:nT], TE_hat[, 1], main = paste("Time Effects (", component_names[1], ")", sep = ""), xlab = "time", ylab = paste("TE (", component_names[1], ")", sep = ""), type = "n", ylim = y1_range)
            graphics::polygon(c(years[2:nT], rev(years[2:nT])), c(lower_bound[, 1], rev(upper_bound[, 1])), col = "grey80", border = NA)
            graphics::grid()
            graphics::abline(h = 0)
            graphics::lines(years[2:nT], TE_hat[, 1], pch = 19, cex = cex_vals_te, type = "b")
            
            # Component 2
            readline(prompt="Press [enter] to see the next plot")
            y2_range <- range(c(lower_bound[, 2], upper_bound[, 2], TE_hat[, 2]), na.rm = TRUE)
            graphics::plot(years[2:nT], TE_hat[, 2], main = paste("Time Effects (", component_names[2], ")", sep = ""), xlab = "time", ylab = paste("TE (", component_names[2], ")", sep = ""), type = "n", ylim = y2_range)
            graphics::polygon(c(years[2:nT], rev(years[2:nT])), c(lower_bound[, 2], rev(upper_bound[, 2])), col = "grey80", border = NA)
            graphics::grid()
            graphics::abline(h = 0)
            graphics::lines(years[2:nT], TE_hat[, 2], pch = 19, cex = cex_vals_te, type = "b")
            
        }
    }

    cat("Plotting finished.\n")
    
    invisible(NULL)
}


#' Performs 2D kernel density estimation on a dataframe and optionally plots the results.
#'
#' This is a wrapper around the `densityEst2d` and `densityEst2dAdaptive` functions.
#' It handles data preparation from a dataframe, supports 1D and 2D variables,
#' and returns a comprehensive results object. For 1D variables, it creates a
#' pseudo-2D representation with a constant second dimension. If `show_plots` or
#' `save_plots` is TRUE, it will also call `plotDensityAnalysis` to generate visualizations.
#'
#' @param df A data frame containing the variable(s).
#' @param var_cols A character vector of length 1 or 2 specifying the variable columns.
#' @param nEval Integer, total number of evaluation grid points (default 2500).
#' @param kernel.type Kernel type for estimation (default "gauss").
#' @param method.h Bandwidth rule (default "silverman").
#' @param h Numeric, bandwidth parameter. If NULL, it's determined by `method.h`.
#' @param adaptive Logical, whether to use adaptive kernel density estimation (default FALSE).
#' @param alpha Numeric, sensitivity parameter for adaptive bandwidth (default 0.5).
#' @param chunk_size Integer chunk size for compute functions. If NULL, defaults to `nEval`.
#' @param component_names A character vector of length 2 with names for the x and y components.
#' @param show_plots Logical, open plotting devices to display plots (default TRUE).
#'
#' @return A list containing all inputs, results, and intermediate objects.
#' @export
runDensityAnalysis <- function(df,
                               var_cols,
                               nEval = 2500,
                               kernel.type = "gauss",
                               method.h = "silverman",
                               h = NULL,
                               adaptive = FALSE,
                               alpha = 0.5,
                               chunk_size = NULL,
                               component_names = NULL,
                               show_plots = TRUE) {
    
    # --- Input validation and data preparation ---
    if (!is.data.frame(df)) {
        stop("`df` must be a data frame.")
    }
    if (!is.character(var_cols) || !length(var_cols) %in% c(1, 2)) {
        stop("`var_cols` must be a character vector of length 1 or 2.")
    }
    if (!all(var_cols %in% colnames(df))) {
        stop("One or more `var_cols` not found in the data frame.")
    }

    is_1d <- length(var_cols) == 1

    X <- df[, var_cols, drop = FALSE]
    X <- X[complete.cases(X), , drop = FALSE]

    if (is_1d) {
        X <- X[[1]]
    } else {
        X <- as.matrix(X)
    }

    nObs <- if(is_1d) length(X) else nrow(X)
    
    if (nObs == 0) {
        stop("No complete observations to perform density estimation.")
    }

    cat("Density analysis started.\n")

    # 1) Evaluation grid
    if (is_1d) {
        x <- defineEvalPoints1d(X, nEval)
    } else {
        x <- defineEvalPoints(X, nEval)
    }
    
    if (is.null(chunk_size)) {
        chunk_size <- nEval
    }

    # 2) Density Estimation
    if (is_1d) {
        if (adaptive) {
            density_results <- densityEst1dAdaptive(X = X, x = x, nEval = nEval, kernel.type = kernel.type,
                                                   method.h = method.h, h = h, alpha = alpha,
                                                   chunk_size = chunk_size, gc = FALSE)
        } else {
            density_results <- densityEst1d(X = X, x = x, nEval = nEval, kernel.type = kernel.type,
                                           method.h = method.h, h = h,
                                           chunk_size = chunk_size, gc = FALSE)
        }
    } else {
        if (adaptive) {
            density_results <- densityEst2dAdaptive(X = X, x = x, nEval = nEval, kernel.type = kernel.type,
                                                   method.h = method.h, h = h, alpha = alpha, 
                                                   chunk_size = chunk_size, gc = FALSE)
        } else {
            density_results <- densityEst2d(X = X, x = x, nEval = nEval, kernel.type = kernel.type,
                                           method.h = method.h, h = h,
                                           chunk_size = chunk_size, gc = FALSE)
        }
    }

    cat("Analysis completed.\n")
    
    # --- Print summary ---
    cat("\n--- Analysis Summary ---\n")
    cat("Estimation Method: Kernel Density", if(adaptive) "(Adaptive)" else "", "\n")
    cat("Kernel Type:", kernel.type, "\n")
    if (is.null(h)) {
        cat("Bandwidth (h) Method:", method.h, "\n")
        cat("  - Estimated h:", round(density_results$h, 4), "\n")
    } else {
        cat("Bandwidth (h) Provided:", round(h, 4), "\n")
    }

    if (adaptive) {
        cat("Alpha:", round(alpha, 4), "\n")
    }
    cat("Input dimensionality:", if(is_1d) "1D" else "2D", "\n")
    cat("------------------------\n\n")


    analysis_results <- listN(
        # Inputs
        df, var_cols, nEval, kernel.type, method.h, h, adaptive, alpha, chunk_size,
        # Data
        X, x, is_1d,
        # Results
        density_results
    )

    if (show_plots) {
        cat("Generating density plot...\n")
        
        # Set up component names
        if (is.null(component_names)) {
            if (is_1d) {
                component_names <- c(var_cols[1], "Density")
            } else {
                component_names <- var_cols
            }
        }
        
        if ((is_1d && length(component_names) != 2) || (!is_1d && length(component_names) != 2)) {
            warning("`component_names` is not valid. Using default names.")
            if(is_1d) component_names <- c(var_cols[1], "Density") else component_names <- var_cols
        }
        
        # Plot
        old_par <- graphics::par(no.readonly = TRUE)
        on.exit(graphics::par(old_par), add = TRUE)
        
        if (is_1d) {
            graphics::plot(x, density_results$estimator, type = 'l',
                          xlab = component_names[1], ylab = component_names[2],
                          main = "Estimated 1D Density")
            graphics::rug(X)
            graphics::grid()
        } else {
            # For 2D, create a contour plot
            nEval_side <- round(sqrt(nEval))
            x1_seq <- x[1:nEval_side, 1]
            x2_seq <- x[seq(1, nrow(x), by = nEval_side), 2]
            density_matrix <- matrix(density_results$estimator, nrow = nEval_side, ncol = nEval_side)
            
            graphics::contour(x1_seq, x2_seq, density_matrix,
                            xlab = component_names[1], ylab = component_names[2],
                            main = "Estimated 2D Density")
            graphics::points(X[, 1], X[, 2], pch = 19, cex = 0.3, col = "red")
            graphics::grid()
        }
        
        cat("Plotting finished.\n")
    }

    return(analysis_results)
}


#' Generates plots from a completed density analysis.
#'
#' This function takes the results from `runDensityAnalysis` and produces
#' visualizations for the estimated density.
#'
#' @param analysis_results A list object returned by `runDensityAnalysis`.
#' @param component_names A character vector of length 2 with names for the x and y components.
#'
#' @return Invisibly returns NULL. This function is called for its side effect of generating plots.
#' @export
plotDensityAnalysis <- function(analysis_results,
                                component_names = NULL) {
    
    # --- Extract data from results object ---
    X <- analysis_results$X
    x <- analysis_results$x
    density_results <- analysis_results$density_results
    is_1d <- analysis_results$is_1d
    var_cols <- analysis_results$var_cols

    if (is.null(component_names)) {
        if (is_1d) {
            component_names <- c(var_cols[1], "Density")
        } else {
            component_names <- var_cols
        }
    }
    
    if ((is_1d && length(component_names) != 2) || (!is_1d && length(component_names) != 2)) {
        warning("`component_names` is not valid. Using default names.")
        if(is_1d) component_names <- c(var_cols[1], "Density") else component_names <- var_cols
    }


    cat("Generating density plot...\n")

    # --- Plotting ---
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par))

    graphics::plot(x, density_results$estimator, type = 'l',
             xlab = component_names[1], ylab = component_names[2],
             main = "Estimated 1D Density")
    graphics::rug(X[, 1])
    graphics::grid()
    
    invisible(NULL)
}


