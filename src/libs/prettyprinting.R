#' Utility functions for debugging and pretty printing optimization results
#' 

#' Prints debug messages only when DEBUG is TRUE
#' 
#' @param message The message to print
#' @param ... Additional arguments for sprintf formatting
debugPrint <- function(message, ...) {
  if (exists("DEBUG") && DEBUG) {
    if (length(list(...)) > 0) {
      print(sprintf(message, ...))
    } else {
      print(message)
    }
  }
}

#' Initializes arrays for storing debug information during optimization
#' 
#' @param dim1 First dimension size (number of h values)
#' @param dim2 Second dimension size (number of alpha values, if applicable)
#' 
#' @return A list of arrays for storing debugging metrics
debugInitArrays <- function(dim1, dim2 = NULL) {
  if (exists("DEBUG") && DEBUG) {
    if (is.null(dim2)) {
      # 1D arrays
      result <- list(
        AICc = array(NA, dim = dim1),
        CovDet = array(NA, dim = dim1),
        trH = array(NA, dim = dim1),
        freedom = array(NA, dim = dim1)
      )
    } else {
      # 2D arrays
      result <- list(
        AICc = array(NA, dim = c(dim1, dim2)),
        CovDet = array(NA, dim = c(dim1, dim2)),
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

#' Stores debug values in arrays during optimization
#' 
#' @param arrays List of arrays to store values in
#' @param i First index (h index)
#' @param j Second index (alpha index, if applicable)
#' @param AICc AICc value to store
#' @param CovDet CovDet value to store
#' @param trH Trace of hat matrix value to store
#' @param freedom Degrees of freedom value to store
#' 
#' @return Updated arrays with new values stored
debugStoreValues <- function(arrays, i, j = NULL, AICc = NA, CovDet = NA, trH = NA, freedom = NA) {
  if (exists("DEBUG") && DEBUG && length(arrays) > 0) {
    if (is.null(j)) {
      # 1D case
      if (!is.na(AICc)) arrays$AICc[i] <- AICc
      if (!is.na(CovDet)) arrays$CovDet[i] <- CovDet
      if (!is.na(trH)) arrays$trH[i] <- trH
      if (!is.na(freedom)) arrays$freedom[i] <- freedom
    } else {
      # 2D case
      if (!is.na(AICc)) arrays$AICc[i,j] <- AICc
      if (!is.na(CovDet)) arrays$CovDet[i,j] <- CovDet
      if (!is.na(trH)) arrays$trH[i,j] <- trH
      if (!is.na(freedom)) arrays$freedom[i,j] <- freedom
    }
  }
  return(arrays)
}

#' Prints optimization results in a formatted way
#' 
#' @param hGrid Vector of h values used in optimization
#' @param h Optimal h value
#' @param arrays List of arrays containing debug information (AICc, CovDet, trH, freedom)
#' @param alpha Optimal alpha value (if applicable)
#' @param alphaGrid Vector of alpha values used in optimization (if applicable)
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
    
    # Print log(CovDet) matrix
    cat("\nlog(CovDet) values for all h and alpha combinations:\n")
    cat("      ") # Space for row labels
    for (j in 1:length(alphaGrid)) {
      cat(sprintf("alpha=%.2f ", alphaGrid[j]))
    }
    cat("\n")
    
    for (i in 1:length(hGrid)) {
      cat(sprintf("h=%.2f ", hGrid[i]))
      for (j in 1:length(alphaGrid)) {
        cat(sprintf("%.2f     ", log(arrays$CovDet[i,j])))
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
    cat(sprintf("%-10s %-10s %-10s %-10s %-10s\n", "h", "trH", "freedom", "log(CovDet)", "AICc"))
    cat(sprintf("%-10s %-10s %-10s %-10s %-10s\n", "----------", "----------", "----------", "----------", "----------"))
    
    # Print each row
    for (i in 1:length(hGrid)) {
      cat(sprintf("%-10.2f %-10.2f %-10.2f %-10.2f %-10.2f\n",
                  hGrid[i], 
                  ifelse(is.null(arrays$trH), NA, arrays$trH[i]), 
                  ifelse(is.null(arrays$freedom), NA, arrays$freedom[i]), 
                  ifelse(is.null(arrays$CovDet), NA, log(arrays$CovDet[i])), 
                  ifelse(is.null(arrays$AICc), NA, arrays$AICc[i])))
    }
    
    # Highlight the optimal value
    cat(sprintf("\nOptimal value: h=%.2f\n", h))
  }
}

#' Creates a custom progress bar with parameter display for optimization
#' 
#' @param min Minimum value of the progress bar (default: 0)
#' @param max Maximum value of the progress bar (default: 1)
#' @param initial Initial value of the progress bar (default: 0)
#' @param width Width of the progress bar in characters (default: console width - 35)
#' 
#' @return A list of functions to update and close the progress bar
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

#' Formats optimization parameters for display in progress bar
#' 
#' @param h Bandwidth parameter value (optional)
#' @param alpha Alpha parameter value for adaptive bandwidth (optional)
#' 
#' @return A formatted string representation of the parameters
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