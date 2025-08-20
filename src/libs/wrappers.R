#' Performs the computationally intensive steps of the panel vector field analysis.
#'
#' This function runs the estimation and bootstrap procedures for panel vector
#' field analysis. It does not produce any plots, but returns an object
#' containing all results needed for plotting. This allows the analysis to be
#' run once, and the plots to be generated or customized separately using
#' `plotPanelVFAnalysis`.
#'
#' @param X 3D numeric array of size nObs x 2 x nT (e.g., log-GDP/log-LE).
#' @param nEval Integer, total number of evaluation grid points (default 2500).
#' @param FE Logical, include fixed effects in within transform (default TRUE).
#' @param TE Logical, include time effects in within transform (default TRUE).
#' @param uniform_weights Logical, use uniform weights for within transform (default TRUE).
#' @param kernel.type Kernel type for estimation (default "epa").
#' @param method.h Bandwidth rule (default "silverman").
#' @param chunk_size Integer chunk size for compute functions (default 512).
#' @param bootstrap_B Integer, number of bootstrap replicates (default 100).
#' @param p_crit Numeric, significance level for directional significance (default 0.05).
#'
#' @return A list containing all estimation and bootstrap results.
#' @export
runPanelVFAnalysis <- function(X,
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
                               alpha = 0.5) {

    # --- Basic checks ---
    if (length(dim(X)) != 3 || dim(X)[2] != 2) {
        stop("X must be a 3D array with second dimension equal to 2 (nObs x 2 x nT).")
    }
    dims <- dim(X)
    nObs <- dims[1]
    nT <- dims[3]

    # --- Top-level progress bar over main stages ---
    show_progress <- !exists("DEBUG") || !DEBUG
    if (show_progress) {
        cat("Starting panel vector field analysis...\n")
        pb <- createCustomProgressBar(min = 0, max = 4)
    } else {
        cat("Panel vector field analysis started.\n")
    }

    # 1) Evaluation grid
    if (show_progress) pb$update(1, "Build evaluation grid")
    X_unrolled <- aperm(X, c(1, 3, 2))
    dim(X_unrolled) <- c(nObs * nT, 2)
    x <- defineEvalPoints(X_unrolled, nEval)

    # --- Estimation ---
    if (FE || TE) {
        # 2) Panel Estimation
        if (show_progress) pb$update(2, "Estimate panel vector field")
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
            sparse = FALSE,
            gc = FALSE
        )
        VF_hat <- panel_vf_results$estimator
        X0_raw <- panel_vf_results$X0_raw

    } else {
        # 2) Standard VF Estimation (no FE/TE)
        if (show_progress) {
            pb$update(2, paste("Estimate VF (", estimation_method, if(adaptive) " adaptive" else "", ")", sep=""))
            if(hOpt || alphaOpt) cat("\n") # Add newline before optimization progress starts
        }
        
        # Prepare data for standard VF functions
        X0 <- aperm(X[, , 1:(nT - 1), drop=FALSE], c(1, 3, 2))
        dim(X0) <- c(nObs * (nT - 1), 2)
        X1 <- aperm(X[, , 2:nT, drop=FALSE], c(1, 3, 2))
        dim(X1) <- c(nObs * (nT - 1), 2)

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
        
        panel_vf_results <- vf_results # To maintain compatibility with bootstrap
        VF_hat <- vf_results$estimator
        X0_raw <- vf_results$X0
    }

    # 3) Bootstrap
    if (show_progress) {
        pb$update(3, paste("Bootstrap B=", bootstrap_B, sep = ""))
        cat("\n") # Add newline before bootstrap progress starts
    }
    if(FE || TE) {
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

    # 4) Effects and significance
    if (show_progress) pb$update(4, "Compute FE/TE and significance")
    
    if (FE || TE) {
        effects <- get_effects(panel_vf_results, X_obs = X0_raw, FE = TRUE, TE = TRUE)
        FE_hat <- effects$alpha_i
        TE_hat <- effects$gamma_t
        FE_signif <- significanceBootstrap(FE_hat, bootstrap_samples$FE_array, p_crit = p_crit)
        TE_signif <- significanceBootstrap(TE_hat, bootstrap_samples$TE_array, p_crit = p_crit)
    } else {
        # No FE/TE to compute or test
        FE_hat <- NULL
        TE_hat <- NULL
        FE_signif <- NULL
        TE_signif <- NULL
    }

    VF_signif <- significanceBootstrap(VF_hat, bootstrap_samples$estimators_array, p_crit = p_crit)
    
    if (show_progress) {
        pb$close()
        cat("Analysis completed.\n")
    } else {
        cat("Analysis completed.\n")
    }

    # --- Print summary ---
    cat("\n--- Analysis Summary ---\n")
    cat("Estimation Method:", estimation_method, if(adaptive) "(Adaptive)" else "", "\n")
    cat("Kernel Type:", kernel.type, "\n")
    
    if (hOpt) {
        h_range <- range(panel_vf_results$hGrid, na.rm = TRUE)
        cat("Bandwidth (h) optimized over [", round(h_range[1], 4), ", ", round(h_range[2], 4), "]\n", sep = "")
        cat("  - Optimal h found:", round(panel_vf_results$h, 4), "\n")
    } else if (is.null(h)) {
        cat("Bandwidth (h) Method:", method.h, "\n")
        cat("  - Estimated h:", round(panel_vf_results$h, 4), "\n")
    } else {
        cat("Bandwidth (h) Provided:", round(h, 4), "\n")
    }

    if (adaptive) {
        if (alphaOpt) {
            alpha_range <- range(panel_vf_results$alphaGrid, na.rm = TRUE)
            cat("Alpha optimized over [", round(alpha_range[1], 4), ", ", round(alpha_range[2], 4), "]\n", sep = "")
            cat("  - Optimal alpha found:", round(panel_vf_results$alpha, 4), "\n")
        } else {
            cat("Alpha:", round(panel_vf_results$alpha, 4), "\n")
        }
    }
    
    cat("Fixed Effects (FE):", FE, "\n")
    cat("Time Effects (TE):", TE, "\n")
    cat("Bootstrap Replicates:", bootstrap_B, "\n")
    cat("Chunk Size:", chunk_size, "\n")
    cat("------------------------\n\n")

    return(listN(
        # Inputs
        X, nEval, FE, TE, uniform_weights, kernel.type, method.h, chunk_size,
        bootstrap_B, p_crit, estimation_method, adaptive, hOpt, alphaOpt, h, alpha,
        # Results
        x, panel_vf_results, bootstrap_samples,
        VF_hat, FE_hat, TE_hat,
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
#' @param out_dir Output directory for saved plots (default "outpics").
#' @param save_plots Logical, save plots as PDFs to `out_dir` (default TRUE).
#' @param show_plots Logical, open plotting devices to display plots (default TRUE).
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
                                out_dir = "outpics",
                                save_plots = TRUE,
                                show_plots = TRUE) {

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

    if (save_plots && !dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
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

    # --- Plotting ---
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par))
    
    # Vector Field plot
    lengthArrows <- (5 / timeInterval) * 1e-1 * 0.5
    if (show_plots) grDevices::dev.new()
    graphics::par(family = "mono")
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
    graphics::points(X0_raw_plot[, 1], X0_raw_plot[, 2], pch = 19, col = "red", cex = 0.5)
    if (requireNamespace("sm", quietly = TRUE)) {
        X0_raw_plot_complete <- X0_raw_plot[complete.cases(X0_raw_plot), ]
        if(nrow(X0_raw_plot_complete) > 1) {
            est.dens <- sm::sm.density(X0_raw_plot_complete, display = "none")
            graphics::contour(est.dens$eval.points[, 1], est.dens$eval.points[, 2], est.dens$estimate,
                              add = TRUE, col = "purple")
        }
    }
    if (save_plots) grDevices::dev.copy2pdf(file = file.path(out_dir, "VF.pdf"), width = 7, height = 7, family = "mono")

    # FE plots
    if (!is.null(FE_hat)) {
        cex_vals_fe <- ifelse(FE_signif, 1.0, 0.1)
        x_init1 <- X[, 1, 1]
        x_init2 <- X[, 2, 1]

        if (show_plots) grDevices::dev.new()
        graphics::par(family = "mono")
        graphics::plot(x_init1, FE_hat[, 1], main = paste("Fixed Effects vs. Initial", component_names[1]), xlab = paste(component_names[1], "at t0"), ylab = paste("FE (", component_names[1], ")", sep = ""), pch = 19, cex = cex_vals_fe)
        graphics::abline(h = 0)
        graphics::grid()
        if (!is.null(label_names) && length(label_names) == nObs) {
            graphics::text(x_init1, FE_hat[, 1], labels = label_names, cex = 0.5, pos = 4, col = "red")
        }
        if (save_plots) grDevices::dev.copy2pdf(file = file.path(out_dir, paste0("FE_vs_initial_", gsub(" ", "_", component_names[1]), ".pdf")), width = 7, height = 7, family = "mono")

        if (show_plots) grDevices::dev.new()
        graphics::par(family = "mono")
        graphics::plot(x_init2, FE_hat[, 2], main = paste("Fixed Effects vs. Initial", component_names[2]), xlab = paste(component_names[2], "at t0"), ylab = paste("FE (", component_names[2], ")", sep = ""), pch = 19, cex = cex_vals_fe)
        graphics::abline(h = 0)
        graphics::grid()
        if (!is.null(label_names) && length(label_names) == nObs) {
            graphics::text(x_init2, FE_hat[, 2], labels = label_names, cex = 0.5, pos = 4, col = "red")
        }
        if (save_plots) grDevices::dev.copy2pdf(file = file.path(out_dir, paste0("FE_vs_initial_", gsub(" ", "_", component_names[2]), ".pdf")), width = 7, height = 7, family = "mono")
    }

    # TE plots with bootstrap bands
    if (!is.null(TE_hat)) {
        cex_vals_te <- ifelse(TE_signif, 1.0, 0.1)
        TE_quantiles <- apply(TE_bootstrap, c(1, 2), stats::quantile, probs = c(0.025, 0.975), na.rm = TRUE)
        lower_bound <- TE_quantiles[1, , ]
        upper_bound <- TE_quantiles[2, , ]

        # Component 1
        if (show_plots) grDevices::dev.new()
        graphics::par(family = "mono")
        y1_range <- range(c(lower_bound[, 1], upper_bound[, 1], TE_hat[, 1]), na.rm = TRUE)
        graphics::plot(years[2:nT], TE_hat[, 1], main = paste("Time Effects (", component_names[1], ")", sep = ""), xlab = "time", ylab = paste("TE (", component_names[1], ")", sep = ""), type = "n", ylim = y1_range)
        graphics::polygon(c(years[2:nT], rev(years[2:nT])), c(lower_bound[, 1], rev(upper_bound[, 1])), col = "grey80", border = NA)
        graphics::grid()
        graphics::abline(h = 0)
        graphics::lines(years[2:nT], TE_hat[, 1], pch = 19, cex = cex_vals_te, type = "b")
        if (save_plots) grDevices::dev.copy2pdf(file = file.path(out_dir, paste0("TE_", gsub(" ", "_", component_names[1]), ".pdf")), width = 7, height = 7, family = "mono")

        # Component 2
        if (show_plots) grDevices::dev.new()
        graphics::par(family = "mono")
        y2_range <- range(c(lower_bound[, 2], upper_bound[, 2], TE_hat[, 2]), na.rm = TRUE)
        graphics::plot(years[2:nT], TE_hat[, 2], main = paste("Time Effects (", component_names[2], ")", sep = ""), xlab = "time", ylab = paste("TE (", component_names[2], ")", sep = ""), type = "n", ylim = y2_range)
        graphics::polygon(c(years[2:nT], rev(years[2:nT])), c(lower_bound[, 2], rev(upper_bound[, 2])), col = "grey80", border = NA)
        graphics::grid()
        graphics::abline(h = 0)
        graphics::lines(years[2:nT], TE_hat[, 2], pch = 19, cex = cex_vals_te, type = "b")
        if (save_plots) grDevices::dev.copy2pdf(file = file.path(out_dir, paste0("TE_", gsub(" ", "_", component_names[2]), ".pdf")), width = 7, height = 7, family = "mono")
    }

    cat("Plotting finished.\n")
    
    invisible(NULL)
}


