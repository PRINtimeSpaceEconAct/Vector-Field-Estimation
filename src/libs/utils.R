listN <- function(...){
    # automatically give names to list elements = var name
    anonList <- list(...)
    names(anonList) <- as.character(substitute(list(...)))[-1]
    anonList
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

# Consolidated AICc calculation function (without input validation)
calculateAICc <- function(X0, X1, X1Hat, lambda, Nobs, kernelFunction, method, 
                          Hkk_values = NULL, # Specific to LL
                          density = NULL, h = NULL, detS = NULL # Specific to NW
                          ) {
    
    # Initialize return values for potential early exit during calculation errors
    error_result <- listN(AICc = Inf, RSS = NA, trH = NA, freedom = NA)

    # --- Calculate tr(H) based on method ---
    trH = NA # Initialize
    if (method == "LL") {
        if (is.null(lambda)) {
            # Original LL formula
            trH = kernelFunction(0, 0) * sum(Hkk_values, na.rm = TRUE)
        } else {
            # Adaptive LL formula 
            trH = kernelFunction(0, 0) * sum( Hkk_values / (lambda^2) , na.rm = TRUE)
        }
    } else if (method == "NW") {
       if (is.null(lambda)) {
            # Original NW formula
            trH = (kernelFunction(0, 0) / (h^2 * Nobs * sqrt(detS))) * sum(1 / density, na.rm = TRUE)
        } else {
            # Adaptive NW formula
            trH = (kernelFunction(0, 0) / (h^2 * Nobs * sqrt(detS))) * sum((1 / lambda^2) * (1 / density), na.rm = TRUE)
       }
    } 

    # --- Calculate degrees of freedom and AICc (common part) ---
    # Use complete.obs for robustness if X1Hat or X1 have NAs
    RSS = det(cov(X1Hat - X1, use="complete.obs")) 

    # Ensure denominator is not zero or negative
    denominator = (1 - (2 * trH + 2) / (2 * Nobs)) 
    
    # Check for invalid calculation results (trH, denominator)
    if (is.na(trH) || is.nan(trH) || is.na(denominator) || is.nan(denominator) || denominator <= 1e-9 || is.na(RSS) || is.nan(RSS)) {
        freedom = Inf
        AICc = Inf
    } else {
       freedom = (1 + (2 * trH) / (2 * Nobs)) / denominator
       AICc = log(RSS) + freedom
    }
    
    # Check for NaN/Inf results after calculation
    if(is.nan(AICc) || is.infinite(AICc)) AICc <- Inf
    if(is.nan(freedom) || is.infinite(freedom)) freedom <- Inf

    return(listN(AICc, RSS, trH, freedom))
}


