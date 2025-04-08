NWregression <- function(X, Y, x=NULL, nEval=2500, kernel.type="gauss", D=NULL, 
                        method.h=NULL, h=NULL, lambda = NULL, 
                        sparse=FALSE, gc=FALSE, chunk_size=nrow(x)) {

    type.est = "NW"
    resultEst = kernelMethod(X=X,x=x,nEval=nEval,kernel.type=kernel.type,D=D,
                             method.h=method.h,h=h,lambda=lambda,
                             sparse=sparse,gc=gc,chunk_size=chunk_size,type.est=type.est, Y=Y)
    
    return(resultEst)
}

# Adaptive bandwidth version of NW regression
NWregressionAdaptive <- function(X, Y, x=NULL, nEval=2500, kernel.type="gauss", D=NULL,
                                method.h=NULL, h=NULL, lambda=NULL,
                                sparse=FALSE, gc=FALSE, chunk_size=nrow(x), alpha = 0.5) {
    
    # Get adaptive bandwidths using the density-based approach
    lambda = getLocalBandwidth(X, kernel.type=kernel.type, D=D, method.h=method.h, h=h,
                             sparse=sparse, gc=gc, chunk_size=chunk_size, alpha = alpha)
    
    # Apply NW regression with adaptive bandwidths
    est = NWregression(X, Y, x=x, nEval=nEval, kernel.type=kernel.type, D=D, 
                      method.h=method.h, h=h, lambda=lambda,
                      sparse=sparse, gc=gc, chunk_size=chunk_size)
    est$lambda = lambda
    return(est)
}

NWfield <- function(X0, X1, x=NULL, nEval=2500, kernel.type="gauss", D=NULL,
                    method.h=NULL, h=NULL, lambda=NULL,
                    sparse=FALSE, gc=FALSE, chunk_size=nrow(x),
                    hOpt = FALSE, nGridh = 10) {

    Y1 = X1[,1] - X0[,1]
    Y2 = X1[,2] - X0[,2]
    
    if (hOpt == TRUE){
        Nobs = nrow(X0)
        covX = cov(X0)
        invS = solve(covX)
        detS = det(covX)
        # define grid of h
        list.h = define_h_method.h(X0, NULL ,"silverman", kernel.type)
        hStart = list.h$h/10
        hEnd = list.h$h*2
        hGrid = exp(log(hStart) + (log(hEnd) - log(hStart)) * (0:(nGridh-1))/(nGridh-1))
        # define kernel function
        kernelFunction = defineKernel(kernel.type)
        
        AICc = array(NA,dim = nGridh)
        RSS = array(NA,dim = nGridh)
        trH = array(NA,dim = nGridh)
        if (DEBUG) trH_old = array(NA,dim = nGridh)
        freedom = array(NA,dim = nGridh)
        for (i in 1:nGridh){
            if (DEBUG) print(paste("Computing h ", i, "/", nGridh, sep=""))
            hi = hGrid[i]
            # est1 = NWregression(X0, Y1, x=x, nEval=nEval, kernel.type=kernel.type, D=D, 
            #                     method.h=method.h, h=hi, lambda=lambda,
            #                     sparse=sparse, gc=gc, chunk_size=chunk_size)
            # est2 = NWregression(X0, Y2, x=est1$x, nEval=nEval, kernel.type=kernel.type, D=D, 
            #                     method.h=method.h, h=hi, lambda=lambda,
            #                     sparse=sparse, gc=gc, chunk_size=chunk_size)
            # VFhi = list(x = est1$x, estimator = cbind(est1$estimator, est2$estimator))
            # X1Hat = forecastDiscrete(X0,VFhi,speedFactor=1,nPeriods=1)
            
            
            # stima alla cazzo fatta sulle obs con un forecast solo
            est1 = NWregression(X0, Y1, x=X0, nEval=nEval, kernel.type=kernel.type, D=D, 
                                method.h=method.h, h=hi, lambda=lambda,
                                sparse=sparse, gc=gc, chunk_size=chunk_size)
            est2 = NWregression(X0, Y2, x=est1$x, nEval=nEval, kernel.type=kernel.type, D=D, 
                                method.h=method.h, h=hi, lambda=lambda,
                                sparse=sparse, gc=gc, chunk_size=chunk_size)
            
            X1Hat = X0 + cbind(est1$estimator, est2$estimator)
            
            if (DEBUG) trH_old[i] = (kernelFunction(0,0)/(hi^2 * Nobs * sqrt(detS))) * sum(1/est1$density)
            trH[i] = kernelFunction(0,0) * sum(1/est1$kernel_sum)
            freedom[i] = (1 + (2*trH[i])/(2*Nobs))/(1 - (2*trH[i]+2)/(2*Nobs))
            RSS[i] = det(cov(X1Hat - X1))
            AICc[i] = log(RSS[i]) + freedom[i]
            
            # maxLength = max(sqrt(rowSums(VFhi$estimator)^2),na.rm=T)
            # nPeriods = ceil(10*maxLength)
            # X1Hat = forecastDiscrete(X0,VFhi,speedFactor=1/nPeriods,nPeriods=nPeriods)
            # RSS[i] = sum(rowSums((X1Hat[,,nPeriods] - X1)^2,na.rm=TRUE))
        }

        h = hGrid[which.min(AICc)]
        if (DEBUG) print(paste("Optimal h: ", format(h, digits=2, nsmall=2)))
        if (DEBUG) print(paste("hGrid: ", paste(format(hGrid, digits=2, nsmall=2), collapse=" ")))
        if (DEBUG) print(paste("trH: ", paste(format(trH, digits=2, nsmall=2), collapse=" ")))
        if (DEBUG) print(paste("trH_old: ", paste(format(trH_old, digits=2, nsmall=2), collapse=" ")))
        if (DEBUG) print(paste("freedom: ", paste(format(freedom, digits=2, nsmall=2), collapse=" ")))
        if (DEBUG) print(paste("RSS: ", paste(format(RSS, digits=2, nsmall=2), collapse=" ")))
        if (DEBUG) print(paste("AICc: ", paste(format(AICc, digits=2, nsmall=2), collapse=" ")))
        
    }
    
    est1 = NWregression(X0, Y1, x=x, nEval=nEval, kernel.type=kernel.type, D=D, 
                      method.h=method.h, h=h, lambda=lambda,
                      sparse=sparse, gc=gc, chunk_size=chunk_size)
    est2 = NWregression(X0, Y2, x=est1$x, nEval=nEval, kernel.type=kernel.type, D=D, 
                      method.h=method.h, h=h, lambda=lambda,
                      sparse=sparse, gc=gc, chunk_size=chunk_size)
    
    # Stack est1$estimator and est2$estimator
    estimator = cbind(est1$estimator, est2$estimator)
    x = est1$x
    density = est1$density
    h = est1$h
    method.h = est1$method.h
    kernel.type = est1$kernel.type
    lambda = est1$lambda
    type.est = "NW"
    if (hOpt == TRUE){
        return(listN(x, X0, X1, estimator, type.est, density, kernel.type, h, method.h, lambda, AICc))
    } else {
        return(listN(x, X0, X1, estimator, type.est, density, kernel.type, h, method.h, lambda))
    }
}

NWfieldAdaptive <- function(X0, X1, x=NULL, nEval=2500, kernel.type="gauss", D=NULL,
                            method.h=NULL, h=NULL, lambda=NULL,
                            sparse=FALSE, gc=FALSE, chunk_size=nrow(x), alpha = 0.5,
                            hOpt = FALSE, nGridh = 10, alphaOpt = FALSE, nGridAlpha = 5) {
    
    # Hardcoded alpha range
    alpha_start = 0
    alpha_end = 0.9
    
    # Parameter validation
    if (hOpt == TRUE && !is.null(h)) {
        stop("Cannot specify h when hOpt is TRUE")
    }
    
    if (alphaOpt == TRUE && !is.null(alpha) && alpha != 0.5) {
        stop("Cannot specify alpha when alphaOpt is TRUE")
    }
    
    Y1 = X1[,1] - X0[,1]
    Y2 = X1[,2] - X0[,2]
    
    if (hOpt == TRUE || alphaOpt == TRUE){
        Nobs = nrow(X0)
        covX = cov(X0)
        invS = solve(covX)
        detS = det(covX)
        # define grid of h
        list.h = define_h_method.h(X0, NULL ,"silverman", kernel.type)
        hStart = list.h$h/10
        hEnd = list.h$h*2
        hGrid = exp(log(hStart) + (log(hEnd) - log(hStart)) * (0:(nGridh-1))/(nGridh-1))
        # Define kernel function
        kernelFunction = defineKernel(kernel.type)
        
        # Create arrays and variables based on whether we're optimizing alpha
        if (alphaOpt) {
            # Define grid of alpha
            alphaGrid = alpha_start + (alpha_end - alpha_start) * (0:(nGridAlpha-1))/(nGridAlpha-1)
            
            # Arrays to store results for all combinations of h and alpha
            AICc_all = array(NA, dim = c(nGridh, nGridAlpha))
            RSS_all = array(NA, dim = c(nGridh, nGridAlpha))
            trH_all = array(NA, dim = c(nGridh, nGridAlpha))
            freedom_all = array(NA, dim = c(nGridh, nGridAlpha))
        } else {
            # Arrays to store results for h only
            AICc_all = array(NA, dim = nGridh)
            RSS_all = array(NA, dim = nGridh)
            trH_all = array(NA, dim = nGridh)
            freedom_all = array(NA, dim = nGridh)
        }
        
        # Best values
        best_AICc = Inf
        best_h = NULL
        best_alpha = if(alphaOpt) NULL else alpha
        best_lambda = NULL
        
        for (i in 1:nGridh){
            if (DEBUG) print(paste("Computing h ", i, "/", nGridh, sep=""))
            hi = hGrid[i]
            
            # Compute pilot density estimation once for this h value
            pilotDensity = densityEst2d(X0, x=X0, nEval=Nobs, kernel.type=kernel.type,
                                       D=D, method.h=method.h, h=hi, sparse=sparse,
                                       gc=gc, chunk_size=chunk_size)
            # Compute g once for this density estimation
            g = exp(mean(log(pilotDensity$estimator)))
            
            if (alphaOpt) {
                # Loop over alpha values if alphaOpt is TRUE
                for (j in 1:nGridAlpha) {
                    if (DEBUG) print(paste("  Computing alpha ", j, "/", nGridAlpha, sep=""))
                    alpha_j = alphaGrid[j]
                    
                    # Compute lambda using the formula without calling getLocalBandwidth again
                    lambda_ij = (pilotDensity$estimator/g)^(-alpha_j)
                    
                    # stima alla cazzo fatta sulle obs con un forecast solo
                    est1 = NWregression(X0, Y1, x=X0, nEval=nEval, kernel.type=kernel.type, D=D,
                                        method.h=method.h, h=hi, lambda=lambda_ij,
                                        sparse=sparse, gc=gc, chunk_size=chunk_size)
                    est2 = NWregression(X0, Y2, x=est1$x, nEval=nEval, kernel.type=kernel.type, D=D,
                                        method.h=method.h, h=hi, lambda=lambda_ij,
                                        sparse=sparse, gc=gc, chunk_size=chunk_size)
                    
                    X1Hat = X0 + cbind(est1$estimator, est2$estimator)
                    
                    trH_all[i,j] = (kernelFunction(0,0)/(hi^2 * Nobs * sqrt(detS))) * sum((1/lambda_ij^2)*1/est1$density)
                    #trH_all[i,j] = kernelFunction(0,0) * sum(1/est1$kernel_sum)
                    freedom_all[i,j] = (1 + (2*trH_all[i,j])/(2*Nobs))/(1 - (2*trH_all[i,j]+2)/(2*Nobs))
                    RSS_all[i,j] = det(cov(X1Hat - X1))
                    AICc_all[i,j] = log(RSS_all[i,j]) + freedom_all[i,j]
                    
                    # Check if this is the best combination so far
                    if (AICc_all[i,j] < best_AICc) {
                        best_AICc = AICc_all[i,j]
                        best_h = hi
                        best_alpha = alpha_j
                        best_lambda = lambda_ij
                    }
                }
            } else {
                # If not optimizing alpha, use the provided alpha value
                alpha_j = alpha
                
                # Compute lambda using the formula without calling getLocalBandwidth again
                lambda_ij = (pilotDensity$estimator/g)^(-alpha_j)
                
                # stima alla cazzo fatta sulle obs con un forecast solo
                est1 = NWregression(X0, Y1, x=X0, nEval=nEval, kernel.type=kernel.type, D=D,
                                    method.h=method.h, h=hi, lambda=lambda_ij,
                                    sparse=sparse, gc=gc, chunk_size=chunk_size)
                est2 = NWregression(X0, Y2, x=est1$x, nEval=nEval, kernel.type=kernel.type, D=D,
                                    method.h=method.h, h=hi, lambda=lambda_ij,
                                    sparse=sparse, gc=gc, chunk_size=chunk_size)
                
                X1Hat = X0 + cbind(est1$estimator, est2$estimator)
                
                trH_all[i] = (kernelFunction(0,0)/(hi^2 * Nobs * sqrt(detS))) * sum((1/lambda_ij^2)*1/est1$density)
                #trH_all[i] = kernelFunction(0,0) * sum(1/est1$kernel_sum)
                freedom_all[i] = (1 + (2*trH_all[i])/(2*Nobs))/(1 - (2*trH_all[i]+2)/(2*Nobs))
                RSS_all[i] = det(cov(X1Hat - X1))
                AICc_all[i] = log(RSS_all[i]) + freedom_all[i]
                
                # Check if this is the best combination so far
                if (AICc_all[i] < best_AICc) {
                    best_AICc = AICc_all[i]
                    best_h = hi
                    best_lambda = lambda_ij
                }
            }
        }

        # Set the optimal values
        h = best_h
        alpha = best_alpha
        lambda = best_lambda
        AICc = AICc_all
        
        if (DEBUG) {
            print(paste("Optimal h:", format(h, digits=2, nsmall=2)))
            print(paste("Optimal alpha:", format(alpha, digits=2, nsmall=2)))
            print(paste("hGrid:", paste(format(hGrid, digits=2, nsmall=2), collapse=" ")))
            
            # Only print alphaGrid if alphaOpt is TRUE
            if (alphaOpt) {
                print(paste("alphaGrid:", paste(format(alphaGrid, digits=2, nsmall=2), collapse=" ")))
                
                # For 2D arrays (when alphaOpt is TRUE), print the full 2D matrices
                
                # Print trH matrix
                cat("\ntrH values for all h and alpha combinations:\n")
                cat("      ") # Space for row labels
                for (j in 1:nGridAlpha) {
                    cat(sprintf("alpha=%.2f ", alphaGrid[j]))
                }
                cat("\n")
                
                for (i in 1:nGridh) {
                    cat(sprintf("h=%.2f ", hGrid[i]))
                    for (j in 1:nGridAlpha) {
                        cat(sprintf("%.2f     ", trH_all[i,j]))
                    }
                    cat("\n")
                }
                
                # Print freedom matrix
                cat("\nfreedom values for all h and alpha combinations:\n")
                cat("      ") # Space for row labels
                for (j in 1:nGridAlpha) {
                    cat(sprintf("alpha=%.2f ", alphaGrid[j]))
                }
                cat("\n")
                
                for (i in 1:nGridh) {
                    cat(sprintf("h=%.2f ", hGrid[i]))
                    for (j in 1:nGridAlpha) {
                        cat(sprintf("%.2f     ", freedom_all[i,j]))
                    }
                    cat("\n")
                }
                
                # Print RSS matrix
                cat("\n log(RSS) values for all h and alpha combinations:\n")
                cat("      ") # Space for row labels
                for (j in 1:nGridAlpha) {
                    cat(sprintf("alpha=%.2f ", alphaGrid[j]))
                }
                cat("\n")
                
                for (i in 1:nGridh) {
                    cat(sprintf("h=%.2f ", hGrid[i]))
                    for (j in 1:nGridAlpha) {
                        cat(sprintf("%.2f     ", log(RSS_all[i,j])))
                    }
                    cat("\n")
                }
                
                # Print AICc matrix
                cat("\nAICc values for all h and alpha combinations:\n")
                cat("      ") # Space for row labels
                for (j in 1:nGridAlpha) {
                    cat(sprintf("alpha=%.2f ", alphaGrid[j]))
                }
                cat("\n")
                
                for (i in 1:nGridh) {
                    cat(sprintf("h=%.2f ", hGrid[i]))
                    for (j in 1:nGridAlpha) {
                        cat(sprintf("%.2f     ", AICc_all[i,j]))
                    }
                    cat("\n")
                }
                
                # Highlight the optimal values
                cat(sprintf("\nOptimal values: h=%.2f, alpha=%.2f\n", h, alpha))
            } else {
                # For 1D arrays (when alphaOpt is FALSE), print in tabular format
                
                # Create a header
                cat("\nValues for all h with fixed alpha=", format(alpha, digits=2, nsmall=2), ":\n", sep="")
                cat(sprintf("%-10s %-10s %-10s %-10s %-10s\n", "h", "trH", "freedom", "RSS", "AICc"))
                cat(sprintf("%-10s %-10s %-10s %-10s %-10s\n", "----------", "----------", "----------", "----------", "----------"))
                
                # Print each row
                for (i in 1:nGridh) {
                    cat(sprintf("%-10.2f %-10.2f %-10.2f %-10.2f %-10.2f\n",
                                hGrid[i], trH_all[i], freedom_all[i], RSS_all[i], AICc_all[i]))
                }
                
                # Highlight the optimal value
                cat(sprintf("\nOptimal value: h=%.2f\n", h))
            }
        }
    } else {
        # If not optimizing, use the provided h and alpha
        # Get adaptive bandwidths using the provided h and alpha
        lambda = getLocalBandwidth(X0, kernel.type=kernel.type, D=D, method.h=method.h, h=h,
                                 sparse=sparse, gc=gc, chunk_size=chunk_size, alpha = alpha)
    }
    
    est1 = NWregression(X0, Y1, x=x, nEval=nEval, kernel.type=kernel.type, D=D,
                      method.h=method.h, h=h, lambda=lambda,
                      sparse=sparse, gc=gc, chunk_size=chunk_size)
    est2 = NWregression(X0, Y2, x=est1$x, nEval=nEval, kernel.type=kernel.type, D=D,
                      method.h=method.h, h=h, lambda=lambda,
                      sparse=sparse, gc=gc, chunk_size=chunk_size)
    
    # Stack est1$estimator and est2$estimator
    estimator = cbind(est1$estimator, est2$estimator)
    x = est1$x
    density = est1$density
    h = est1$h
    method.h = est1$method.h
    kernel.type = est1$kernel.type
    type.est = "NW"
    
    # Create result object
    result = listN(x, X0, X1, estimator, type.est, density, kernel.type, h, method.h, lambda)
    
    # Add optimization results
    if (hOpt == TRUE || alphaOpt == TRUE){
        result$AICc = AICc
        result$hGrid = hGrid
        
        # Add alpha-related information only if alphaOpt is TRUE
        if (alphaOpt == TRUE) {
            result$alpha = alpha
            result$alphaGrid = alphaGrid
        } else if (hOpt == TRUE) {
            # If only hOpt is TRUE, still include alpha but not alphaGrid
            result$alpha = alpha
        }
    }
    
    return(result)
}

