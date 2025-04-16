listN <- function(...){
    # automatically give names to list elements = var name
    anonList <- list(...)
    names(anonList) <- as.character(substitute(list(...)))[-1]
    anonList
}

determinant3 <- function(a00, a01, a02, a10, a11, a12, a20, a21, a22){
    return(a00 * (a11 * a22 - a12 * a21) - a01 * (a10 * a22 - a12 * a20) + a02 * (a10 * a21 - a11 * a20))
}

defineEvalPoints <- function(X,nEval){
    xGrid = seq(from=min(X[,1]), to=max(X[,1]), length.out=round(sqrt(nEval)))
    yGrid = seq(from=min(X[,2]), to=max(X[,2]), length.out=round(sqrt(nEval)))
    x = as.matrix(expand.grid(xGrid,yGrid))
    return(x)
}

define_h_method.h <- function(X=X, h=h, method.h=method.h, kernel.type="gauss") {
    # Check for conflicting bandwidth specifications
    if (!is.null(h) && !is.null(method.h)) {
        stop("Cannot specify both h and method.h. Please provide only one.")
    }
    
    # Define constant based on kernel type
    c = if (kernel.type == "epa") 1.77 else 0.96
    
    # If h is provided directly, use it regardless of method.h
    if (!is.null(h)) {
        # Keep h as is
        if (is.null(method.h)) method.h <- "manual" # Indicate h was manually set
    } else if (is.null(method.h)) {
        # Default method when both h and method.h are null
        debugPrint("No h or method.h provided, using default: silverman")
        method.h = "silverman"
        h = c * nrow(X)^(-1/6)
    } else {
        # Select h based on specified method
        debugPrint("Using method.h = %s to determine h", method.h)
        h = switch(method.h,
                   "silverman" = c * nrow(X)^(-1/6),
                   #"sj" = { requireNamespace("KernSmooth", quietly = TRUE); KernSmooth::dpik(X) }, # 1D only
                   # "botev" = botev(X), # 1D only
                   stop(paste("method.h '", method.h, "' not recognized", sep=""))
        )
    }
    
    debugPrint("Determined h = %.4f using method = %s", h, method.h)
    
    return(listN(h, method.h))
}

interp2d <- function(z,x,VF){
    # z = where to interpolate <- (nInterp x 2)
    # x = evaluation points where f is computed <- (nEval x 2)
    # VF = Vector Field to interpolate <- (nEval x 2)
    
    VFz1 = interp1d(z,x,VF[,1])
    VFz2 = interp1d(z,x,VF[,2])
    
    return(cbind(VFz1,VFz2))
}

interp1d <- function(z,x,f){
    # z = where to interpolate <- (nInterp x 2)
    # x = evaluation points where f is computed <- (nEval x 2)
    # f = function to interpolate <- (nEval)
    
    xS = sort(unique(x[,1]))
    yS = sort(unique(x[,2]))
    nEval = length(xS)
    interpZ = interp2(xS,yS,matrix(f,length(xS),length(yS),byrow=T),z[,1],z[,2])
    return(interpZ)
}

# Initialize variables used for optimization loops and storing results
initializeOptimizationVariables <- function() {
  list(
    AICc_values = NULL,
    hGrid = NULL,
    alphaGrid = NULL,
    best_AICc = Inf,
    best_h = NULL,
    best_alpha = NULL,
    best_lambda = NULL,
    debugArrays = list() # Initialize as empty list
  )
}

# Setup common parameters needed for optimization (AICc calculation)
# Returns a list containing necessary parameters like Nobs, detS, grids, kernelFunction, etc.
setupOptimizationParameters <- function(X0, kernel.type,
                                        hOpt, nGridh, h = NULL, method.h = NULL,
                                        alphaOpt = FALSE, nGridAlpha = 5, alpha = NULL,
                                        lambda = NULL, isAdaptive = FALSE) {
  # --- Begin Input Validation ---
  # 1. Validate lambda usage
  if (!isAdaptive && !is.null(lambda)) {
    stop("Non-adaptive methods do not support bandwidths (lambda). Use adaptive method instead.")
  }
  
  # 2. Validate h/hOpt combinations
  if (hOpt && !is.null(h)) {
    stop("Cannot specify h when hOpt is TRUE.")
  }
  
  # 3. Validate method.h with hOpt 
  if (hOpt && !is.null(method.h)) {
    warning("method.h is ignored when hOpt is TRUE.")
    method.h <- NULL  # Prevent it being used later
  }
  
  # 4. Validate alpha/alphaOpt combinations (only relevant for adaptive methods)
  if (isAdaptive) {
    if (!alphaOpt && is.null(alpha)) {
      stop("Must specify alpha or set alphaOpt = TRUE.")
    }
    if (alphaOpt && !is.null(alpha) && alpha != 0.5) {
      warning("Specified 'alpha' is ignored when alphaOpt is TRUE (optimization determines best alpha).")
    }
    
    # 5. Validate lambda with optimization flags
    if (!hOpt && !alphaOpt && !is.null(lambda)) {
      stop("Cannot specify lambda directly when not optimizing h or alpha (it's calculated from h and alpha).")
    }
  }
  # --- End Input Validation ---

  # Continue with existing setup logic
  Nobs = nrow(X0)
  covX = cov(X0)
  detS = det(covX)
  kernelFunction = defineKernel(kernel.type)

  hGrid = NULL
  nGridh_actual = 1 # Default if hOpt is FALSE

  if (hOpt) {
    # Calculate hGrid if optimizing h
    # Use silverman for the pilot range calculation, consistent with previous logic
    list.h.pilot = define_h_method.h(X0, NULL, "silverman", kernel.type)
    hStart = list.h.pilot$h / 10
    hEnd = list.h.pilot$h * 2
    hGrid = exp(log(hStart) + (log(hEnd) - log(hStart)) * (0:(nGridh - 1)) / (nGridh - 1))
    nGridh_actual = nGridh
  } else {
    # If h is fixed (hOpt is FALSE), determine the single h value
    # Use define_h_method.h to handle cases where h is NULL but method.h is not
    list.h.fixed = define_h_method.h(X0, h = h, method.h = method.h, kernel.type = kernel.type)
    hGrid = list.h.fixed$h # hGrid contains the single fixed h value
    h = list.h.fixed$h # Ensure h variable holds the determined value
    method.h = list.h.fixed$method.h # Ensure method.h reflects the used method
    nGridh_actual = 1
  }

  # Determine alpha grid and dimensions (only relevant for adaptive methods)
  alphaGrid <- NULL
  nGridAlpha_actual <- 1 # Default if alphaOpt is FALSE
  
  if (isAdaptive) {
    alpha_start = 0 # Hardcoded, consider making params
    alpha_end = 0.9 # Hardcoded, consider making params

    if (alphaOpt) {
      alphaGrid = alpha_start + (alpha_end - alpha_start) * (0:(nGridAlpha - 1)) / (nGridAlpha - 1)
      nGridAlpha_actual = nGridAlpha
    } else if (!is.null(alpha)) {
      # If alpha is fixed
      alphaGrid = alpha
      # nGridAlpha_actual remains 1
    }
  }

  # Return all potentially useful parameters
  return(listN(Nobs, detS, hGrid, kernelFunction, nGridh_actual, 
               alphaGrid, nGridAlpha_actual, h, method.h, alpha, isAdaptive))
}

# Debug utilities

# Debug print function - only prints if DEBUG is TRUE
debugPrint <- function(message, ...) {
  if (exists("DEBUG") && DEBUG) {
    if (length(list(...)) > 0) {
      print(sprintf(message, ...))
    } else {
      print(message)
    }
  }
}

# Creates arrays for debug information
debugInitArrays <- function(dim1, dim2 = NULL) {
  if (exists("DEBUG") && DEBUG) {
    if (is.null(dim2)) {
      # 1D arrays
      result <- list(
        AICc = array(NA, dim = dim1),
        RSS = array(NA, dim = dim1),
        trH = array(NA, dim = dim1),
        freedom = array(NA, dim = dim1)
      )
    } else {
      # 2D arrays
      result <- list(
        AICc = array(NA, dim = c(dim1, dim2)),
        RSS = array(NA, dim = c(dim1, dim2)),
        trH = array(NA, dim = c(dim1, dim2)),
        freedom = array(NA, dim = c(dim1, dim2))
      )
    }
    return(result)
  } else {
    # Return empty list if not in debug mode
    return(list())
  }
}

# Store debug values in arrays
debugStoreValues <- function(arrays, i, j = NULL, AICc = NA, RSS = NA, trH = NA, freedom = NA) {
  if (exists("DEBUG") && DEBUG && length(arrays) > 0) {
    if (is.null(j)) {
      # 1D case
      if (!is.na(AICc)) arrays$AICc[i] <- AICc
      if (!is.na(RSS)) arrays$RSS[i] <- RSS
      if (!is.na(trH)) arrays$trH[i] <- trH
      if (!is.na(freedom)) arrays$freedom[i] <- freedom
    } else {
      # 2D case
      if (!is.na(AICc)) arrays$AICc[i,j] <- AICc
      if (!is.na(RSS)) arrays$RSS[i,j] <- RSS
      if (!is.na(trH)) arrays$trH[i,j] <- trH
      if (!is.na(freedom)) arrays$freedom[i,j] <- freedom
    }
  }
  return(arrays)
}

# Print optimization results for both 1D and 2D cases
printOptimizationResults <- function(hGrid, h, arrays = NULL, alpha = NULL, alphaGrid = NULL) {
  if (!exists("DEBUG") || !DEBUG) return()
  
  print(paste("Optimal h:", format(h, digits=2, nsmall=2)))
  if (!is.null(alpha)) {
    print(paste("Optimal alpha:", format(alpha, digits=2, nsmall=2)))
  }
  print(paste("hGrid:", paste(format(hGrid, digits=2, nsmall=2), collapse=" ")))
  
  # 2D case - alpha was optimized (alphaGrid exists), regardless of final optimal alpha
  if (!is.null(alphaGrid)) {
    print(paste("alphaGrid:", paste(format(alphaGrid, digits=2, nsmall=2), collapse=" ")))
    
    if (is.null(arrays) || length(arrays) == 0) return()
    
    # Print trH matrix
    cat("\ntrH values for all h and alpha combinations:\n")
    cat("      ") # Space for row labels
    for (j in 1:length(alphaGrid)) {
      cat(sprintf("alpha=%.2f ", alphaGrid[j]))
    }
    cat("\n")
    
    for (i in 1:length(hGrid)) {
      cat(sprintf("h=%.2f ", hGrid[i]))
      for (j in 1:length(alphaGrid)) {
        cat(sprintf("%.2f     ", arrays$trH[i,j]))
      }
      cat("\n")
    }
    
    # Print freedom matrix
    cat("\nfreedom values for all h and alpha combinations:\n")
    cat("      ") # Space for row labels
    for (j in 1:length(alphaGrid)) {
      cat(sprintf("alpha=%.2f ", alphaGrid[j]))
    }
    cat("\n")
    
    for (i in 1:length(hGrid)) {
      cat(sprintf("h=%.2f ", hGrid[i]))
      for (j in 1:length(alphaGrid)) {
        cat(sprintf("%.2f     ", arrays$freedom[i,j]))
      }
      cat("\n")
    }
    
    # Print log(RSS) matrix
    cat("\nlog(RSS) values for all h and alpha combinations:\n")
    cat("      ") # Space for row labels
    for (j in 1:length(alphaGrid)) {
      cat(sprintf("alpha=%.2f ", alphaGrid[j]))
    }
    cat("\n")
    
    for (i in 1:length(hGrid)) {
      cat(sprintf("h=%.2f ", hGrid[i]))
      for (j in 1:length(alphaGrid)) {
        cat(sprintf("%.2f     ", log(arrays$RSS[i,j])))
      }
      cat("\n")
    }
    
    # Print AICc matrix
    cat("\nAICc values for all h and alpha combinations:\n")
    cat("      ") # Space for row labels
    for (j in 1:length(alphaGrid)) {
      cat(sprintf("alpha=%.2f ", alphaGrid[j]))
    }
    cat("\n")
    
    for (i in 1:length(hGrid)) {
      cat(sprintf("h=%.2f ", hGrid[i]))
      for (j in 1:length(alphaGrid)) {
        cat(sprintf("%.2f     ", arrays$AICc[i,j]))
      }
      cat("\n")
    }
    
    # Highlight the optimal values
    cat(sprintf("\nOptimal values: h=%.2f, alpha=%.2f\n", h, alpha))
  } 
  # 1D case - only h was optimized OR alpha was fixed
  else {
    if (is.null(arrays) || length(arrays) == 0) return()
    
    # For 1D arrays (when alphaOpt is FALSE), print in tabular format
    cat("\nValues for all h", ifelse(!is.null(alpha), paste(" with fixed alpha=", format(alpha, digits=2, nsmall=2), ":", sep=""), ":"), "\n")
    cat(sprintf("%-10s %-10s %-10s %-10s %-10s\n", "h", "trH", "freedom", "log(RSS)", "AICc"))
    cat(sprintf("%-10s %-10s %-10s %-10s %-10s\n", "----------", "----------", "----------", "----------", "----------"))
    
    # Print each row
    for (i in 1:length(hGrid)) {
      cat(sprintf("%-10.2f %-10.2f %-10.2f %-10.2f %-10.2f\n",
                  hGrid[i], 
                  ifelse(is.null(arrays$trH), NA, arrays$trH[i]), 
                  ifelse(is.null(arrays$freedom), NA, arrays$freedom[i]), 
                  ifelse(is.null(arrays$RSS), NA, log(arrays$RSS[i])), 
                  ifelse(is.null(arrays$AICc), NA, arrays$AICc[i])))
    }
    
    # Highlight the optimal value
    cat(sprintf("\nOptimal value: h=%.2f\n", h))
  }
}

# Custom progress bar that can display parameter values
createCustomProgressBar <- function(min = 0, max = 1, initial = 0, 
                                   width = getOption("width") - 35) {
  # Initialize environment to store progress bar state
  pb <- new.env(parent = emptyenv())
  pb$min <- min
  pb$max <- max
  pb$width <- width
  pb$current <- initial
  pb$last_drawn <- -1
  pb$start_time <- Sys.time()
  pb$params <- ""
  
  # Function to set current value
  update_progress <- function(value, params = NULL) {
    pb$current <- value
    if (!is.null(params)) {
      pb$params <- params
    }
    
    # Only redraw if significant progress was made or parameters changed
    if (pb$current != pb$last_drawn || !is.null(params)) {
      # Calculate percentage complete
      pct <- (pb$current - pb$min) / (pb$max - pb$min)
      pct <- max(0, min(1, pct))  # Ensure between 0 and 1
      
      # Calculate how many characters to fill
      chars <- floor(pct * pb$width)
      
      # Calculate time elapsed and ETA
      elapsed <- difftime(Sys.time(), pb$start_time, units = "secs")
      if (pct > 0) {
        eta <- as.numeric(elapsed) * (1 - pct) / pct
        eta_text <- sprintf("ETA: %02d:%02d", floor(eta / 60), floor(eta %% 60))
      } else {
        eta_text <- "ETA: --:--"
      }
      
      # Build progress bar string
      bar <- paste0(
        "[", 
        paste(rep("=", chars), collapse = ""),
        ">",
        paste(rep(" ", pb$width - chars), collapse = ""),
        "] ",
        sprintf("%3d%%", floor(pct * 100)),
        " | ", eta_text
      )
      
      # Add parameter information if provided
      if (pb$params != "") {
        # Clear current line and print progress + parameters
        cat("\r", bar, " | ", pb$params, "    ", sep = "")
      } else {
        # Clear current line and print progress only
        cat("\r", bar, "    ", sep = "")
      }
      
      pb$last_drawn <- pb$current
    }
  }
  
  # Function to close progress bar
  close_progress <- function() {
    cat("\n")
  }
  
  # Return functions for usage
  list(
    update = update_progress,
    close = close_progress
  )
}

# Helper function to format optimization loop parameters for display
formatOptimParams <- function(h = NULL, alpha = NULL) {
  params <- ""
  if (!is.null(h)) {
    params <- paste0(params, "h=", format(h, digits=4, nsmall=4))
  }
  if (!is.null(alpha)) {
    if (params != "") params <- paste0(params, ", ")
    params <- paste0(params, "α=", format(alpha, digits=2, nsmall=2))
  }
  return(params)
}


