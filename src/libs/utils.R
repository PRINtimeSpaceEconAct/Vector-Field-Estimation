#' Creates a named list from variables using their names
#' 
#' @param ... Variables to include in the list
#' 
#' @return A list with names matching the variable names
listN <- function(...){
    # automatically give names to list elements = var name
    anonList <- list(...)
    names(anonList) <- as.character(substitute(list(...)))[-1]
    anonList
}

#' Defines a grid of evaluation points over a 2D domain
#' 
#' @param X Matrix of points that define the domain boundaries (nObs x 2)
#' @param nEval Total number of evaluation points to generate
#' 
#' @return A matrix of evaluation points arranged in a grid (nEval x 2)
defineEvalPoints <- function(X,nEval){
    xGrid = seq(from=min(X[,1]), to=max(X[,1]), length.out=round(sqrt(nEval)))
    yGrid = seq(from=min(X[,2]), to=max(X[,2]), length.out=round(sqrt(nEval)))
    x = as.matrix(expand.grid(xGrid,yGrid))
    return(x)
}

#' Determines bandwidth parameter h and method for kernel estimation
#' 
#' @param X Matrix of data points (nObs x 2)
#' @param h Bandwidth parameter (if NULL, determined by method.h)
#' @param method.h Method for bandwidth selection (if NULL and h is NULL, defaults to "silverman")
#' @param kernel.type Type of kernel function to use (default: "gauss")
#' 
#' @return A list containing:
#'   \item{h}{Bandwidth parameter determined}
#'   \item{method.h}{Method used for bandwidth selection}
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

#' Interpolates a 2D vector field at specified points
#' 
#' @param z Matrix of points where to interpolate (nInterp x 2)
#' @param x Matrix of evaluation points where the vector field is known (nEval x 2)
#' @param VF Matrix of vector field values at x (nEval x 2)
#' 
#' @return Matrix of interpolated vector field values at z (nInterp x 2)
interp2d <- function(z,x,VF){
    # z = where to interpolate <- (nInterp x 2)
    # x = evaluation points where f is computed <- (nEval x 2)
    # VF = Vector Field to interpolate <- (nEval x 2)
    
    VFz1 = interp1d(z,x,VF[,1])
    VFz2 = interp1d(z,x,VF[,2])
    
    return(cbind(VFz1,VFz2))
}

#' Interpolates a scalar function at specified points in 2D
#' 
#' @param z Matrix of points where to interpolate (nInterp x 2)
#' @param x Matrix of evaluation points where the function is known (nEval x 2)
#' @param f Vector of function values at x (nEval)
#' 
#' @return Vector of interpolated function values at z (nInterp)
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

#' Initializes variables used for optimization in kernel methods
#' 
#' @return A list of initialized optimization variables:
#'   \item{AICc_values}{NULL, will store AICc values during optimization}
#'   \item{hGrid}{NULL, will store grid of h values}
#'   \item{alphaGrid}{NULL, will store grid of alpha values}
#'   \item{best_AICc}{Inf, will store best AICc value found}
#'   \item{best_h}{NULL, will store best h value found}
#'   \item{best_alpha}{NULL, will store best alpha value found}
#'   \item{best_lambda}{NULL, will store best lambda value found}
#'   \item{debugArrays}{Empty list for debug information}
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

#' Sets up parameters for bandwidth optimization in kernel methods
#' 
#' @param X0 Matrix of initial points (nObs x 2)
#' @param kernel.type Type of kernel function to use
#' @param hOpt Whether to optimize bandwidth h
#' @param nGridh Number of grid points for h optimization
#' @param h Bandwidth parameter (if NULL and not optimizing, determined by method.h)
#' @param method.h Method for bandwidth selection
#' @param alphaOpt Whether to optimize alpha parameter for adaptive methods
#' @param nGridAlpha Number of grid points for alpha optimization
#' @param alpha Sensitivity parameter for adaptive bandwidth
#' @param lambda Vector of local bandwidths (nObs, for adaptive methods)
#' @param isAdaptive Whether using adaptive bandwidth estimation
#' 
#' @return A list containing all parameters needed for optimization:
#'   \item{Nobs}{Number of observations}
#'   \item{detS}{Determinant of covariance matrix}
#'   \item{hGrid}{Grid of h values for optimization (nGridh)}
#'   \item{kernelFunction}{Kernel function to use}
#'   \item{nGridh_actual}{Actual number of h grid points}
#'   \item{alphaGrid}{Grid of alpha values for optimization (nGridAlpha)}
#'   \item{nGridAlpha_actual}{Actual number of alpha grid points}
#'   \item{h}{Current h value}
#'   \item{method.h}{Method used for h selection}
#'   \item{alpha}{Current alpha value}
#'   \item{isAdaptive}{Whether using adaptive bandwidth}
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

#' Calculates Akaike Information Criterion with correction (AICc) for model selection
#' 
#' @param X0 Matrix of initial points (nObs x 2)
#' @param X1 Matrix of terminal points (nObs x 2)
#' @param X1Hat Matrix of predicted terminal points (nObs x 2)
#' @param lambda Vector of local bandwidths (nObs, for adaptive methods)
#' @param Nobs Number of observations
#' @param kernelFunction Kernel function used in estimation
#' @param method Estimation method ("LL" or "NW")
#' @param Hkk_values Hat matrix diagonal values (nObs, required for LL method)
#' @param density Density estimates at evaluation points (nObs, required for NW method)
#' @param h Bandwidth parameter (required for NW method)
#' @param detS Determinant of covariance matrix (required for NW method)
#' 
#' @return A list containing:
#'   \item{AICc}{Akaike Information Criterion with correction}
#'   \item{CovDet}{Determinant of the covariance matrix of residuals}
#'   \item{trH}{Trace of hat matrix}
#'   \item{freedom}{Effective degrees of freedom}
calculateAICc <- function(X0, X1, X1Hat, lambda, Nobs, kernelFunction, method, 
                          Hkk_values = NULL, # Specific to LL
                          density = NULL, h = NULL, detS = NULL # Specific to NW
                          ) {
    
    # Initialize return values for potential early exit during calculation errors
    error_result <- listN(AICc = Inf, CovDet = NA, trH = NA, freedom = NA)

    # --- Calculate tr(H) based on method ---
    trH = NA # Initialize
    if (method == "LL") {
        if (is.null(lambda)) {
            # Original LL formula
            trH = min( kernelFunction(0, 0) * sum(Hkk_values, na.rm = TRUE)/sum(!is.nan(Hkk_values))*Nobs ,Nobs-2)
        } else {
            # Adaptive LL formula 
            trH = min(kernelFunction(0, 0) * sum( Hkk_values / (lambda^2) , na.rm = TRUE)/sum(!is.nan(Hkk_values))*Nobs,Nobs-2)
        }
    } else if (method == "NW") {
       if (is.null(lambda)) {
            # Original NW formula
            trH = min((kernelFunction(0, 0) / (h^2 * Nobs * sqrt(detS))) * sum(1 / density, na.rm = TRUE) / sum(!is.nan(density))*Nobs,Nobs-2)
        } else {
            # Adaptive NW formula
            trH = min((kernelFunction(0, 0) / (h^2 * Nobs * sqrt(detS))) * sum((1 / lambda^2) * (1 / density), na.rm = TRUE) / sum(!is.nan(density))*Nobs,Nobs-2)
       }
    } 

    # --- Calculate degrees of freedom and AICc (common part) ---
    # Use complete.obs for robustness if X1Hat or X1 have NAs
    CovDet = det(cov(X1Hat - X1, use="complete.obs")) 

    # Ensure denominator is not zero or negative
    denominator = 1 - (2 * trH + 2) / (2 * Nobs)
    
    # Check for invalid calculation results (trH, denominator)
    if (is.na(trH) || is.nan(trH) || is.na(denominator) || is.nan(denominator) || denominator <= 1e-9 || is.na(CovDet) || is.nan(CovDet)) {
        freedom = Inf
        AICc = Inf
    } else {
       freedom = (1 + (2 * trH) / (2 * Nobs)) / denominator
       AICc = log(CovDet) + freedom
    }
    
    # Check for NaN/Inf results after calculation
    if(is.nan(AICc) || is.infinite(AICc)) AICc <- Inf
    if(is.nan(freedom) || is.infinite(freedom)) freedom <- Inf

    return(listN(AICc, CovDet, trH, freedom))
}

#' Converts a panel data frame to a 3D array
#'
#' This function reshapes a data frame in long panel format to a 3D array
#' of dimensions (n_individuals, n_variables, n_time_periods). It assumes
#' a balanced panel.
#'
#' @param panel_df A data frame in long format.
#' @param var_cols A character vector of column names for the variables to include in the array.
#' @param id_col A string specifying the name of the individual identifier column.
#' @param time_col A string specifying the name of the time identifier column.
#'
#' @return A 3D array with the panel data, with dimensions corresponding to individuals, variables, and time.
#' @import plm
panel_df_to_array <- function(panel_df, var_cols, id_col, time_col) {
  # Ensure plm package is available
  if (!requireNamespace("plm", quietly = TRUE)) {
    stop("Package 'plm' is required for this function but is not installed.")
  }
  
  # 1. Check: Must have exactly two variables.
  if (length(var_cols) != 2) {
    stop("This function is designed for 2D vector fields and requires exactly two variables in `var_cols`.")
  }
  
  # 2. Check: No missing values in the specified variable columns.
  if (any(is.na(panel_df[, var_cols]))) {
    stop("Missing values (NA/NaN) detected in `var_cols`. Please fill or remove them before conversion.")
  }
  
  # Create pdata.frame for plm checks. This also checks for duplicate (id, time) pairs.
  p_df <- plm::pdata.frame(panel_df, index = c(id_col, time_col), drop.index = FALSE)
  
  # 3. Check: The panel must be balanced using plm.
  if (!plm::is.pbalanced(p_df)) {
    stop("The panel data is unbalanced. Each individual must have an observation for each time point. Consider using `plm::make.pbalanced`.")
  }
  
  # 4. Check: Time intervals must be equally spaced.
  time_points <- sort(unique(panel_df[[time_col]]))
  
  # Attempt to convert time column to numeric if it isn't already
  if (!is.numeric(time_points)) {
    numeric_time_points <- suppressWarnings(as.numeric(as.character(time_points)))
    if (any(is.na(numeric_time_points))) {
      stop("Time column contains non-numeric values and could not be converted for interval checking.")
    }
    time_points <- numeric_time_points
  }
  
  if (length(time_points) > 1) {
    intervals <- diff(time_points)
    if (length(unique(intervals)) > 1) {
      stop("Time intervals are not equally spaced. All time steps must be of the same duration.")
    }
  }
  
  # Get dimensions
  ids <- sort(unique(panel_df[[id_col]]))
  nObs <- length(ids)
  nT <- length(time_points)
  nVars <- length(var_cols)
  
  # Create the array with named dimensions for clarity
  X <- array(NA, 
             dim = c(nObs, nVars, nT), 
             dimnames = list(as.character(ids), var_cols, as.character(time_points)))
  
  # Populate the array efficiently
  for (t_idx in 1:nT) {
    time_val <- time_points[t_idx]
    
    # Subset data for the current time point
    time_slice <- panel_df[panel_df[[time_col]] == time_val, ]
    
    # Order the slice by ID to ensure correct assignment into the array
    time_slice_ordered <- time_slice[match(ids, time_slice[[id_col]]), ]
    
    X[, , t_idx] <- as.matrix(time_slice_ordered[, var_cols])
  }
  
  return(X)
}

#' Processes input when two dataframes are provided for direct X0, X1 comparison.
#'
#' This function is called by `runPanelVFAnalysis` when `panel_df` is a list of
#' two dataframes. It validates the input and prepares the `X0`, `X1`, and a
#' compatible `X` array for the non-panel estimation workflow.
#'
#' @param panel_df A list containing two dataframes of the same dimensions.
#' @param var_cols A character vector of column names for the variables.
#'
#' @return A list containing:
#'   \item{X0}{A matrix of initial states from the first dataframe.}
#'   \item{X1}{A matrix of final states from the second dataframe.}
#'   \item{X}{A placeholder 3D array for compatibility with downstream functions.}
#'   \item{dims}{Dimensions of the placeholder X array.}
#'   \item{nObs}{Number of observations.}
#'   \item{nT}{Number of time points (always 2 in this case).}
#'   \item{X_unrolled}{The X0 matrix, used for defining evaluation points.}
handle_two_df_input <- function(panel_df, var_cols) {
    df0 <- panel_df[[1]]
    df1 <- panel_df[[2]]

    if (!identical(dim(df0), dim(df1))) {
        stop("The two dataframes provided in panel_df must have the same dimensions.")
    }

    if (!all(var_cols %in% colnames(df0)) || !all(var_cols %in% colnames(df1))) {
        stop("var_cols not found in both dataframes.")
    }

    X0 <- as.matrix(df0[, var_cols])
    X1 <- as.matrix(df1[, var_cols])

    # Create a placeholder X for compatibility.
    X <- array(c(X0, X1), dim = c(nrow(X0), 2, 2))
    dims <- dim(X)
    nObs <- dims[1]
    nT <- dims[3]
    X_unrolled <- X0 # For eval points, use the starting points.

    return(listN(X0, X1, X, dims, nObs, nT, X_unrolled))
}

