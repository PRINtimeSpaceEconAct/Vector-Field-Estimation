NWregression <- function(X, Y, x=NULL, nEval=2500, kernel="epa", D=NULL, 
                        method.h=NULL, h=NULL, lambda = NULL, 
                        sparse=FALSE, gc=FALSE, chunk_size=1024) {

    if (is.null(x)){
        xGrid = seq(from=min(X[,1]), to=max(X[,1]), length.out=round(sqrt(nEval)))
        yGrid = seq(from=min(X[,2]), to=max(X[,2]), length.out=round(sqrt(nEval)))
        x = as.matrix(expand.grid(xGrid,yGrid))
    }
    nEval = nrow(x)
    
    # Compute numerator (sum of Y_i * K_h)
    numerator = densityEst2d(X, x=x, nEval=nEval, kernel=kernel, D=D, 
                            method.h=method.h, h=h, lambda=lambda,
                            sparse=sparse, gc=gc, chunk_size=chunk_size, Y=Y)
    
    # Compute denominator (sum of K_h)
    denominator = densityEst2d(X, x=x, nEval=nEval, kernel=kernel, D=D, 
                              method.h=method.h, h=h, lambda=lambda,
                              sparse=sparse, gc=gc, chunk_size=chunk_size)
    
    # Compute final NW estimate
    NWest = numerator$densityEst / denominator$densityEst
    
    return(listN(x, NWest))
}

# Adaptive bandwidth version of NW regression
NWregressionAdaptive <- function(X, Y, x=NULL, nEval=2500, kernel="epa", D=NULL,
                                method.h=NULL, h=NULL,
                                sparse=FALSE, gc=FALSE, chunk_size=1024, alpha = 0.5) {
    
    # Get adaptive bandwidths using the density-based approach
    lambda = getLocalBandwidth(X, kernel=kernel, D=D, method.h=method.h, h=h,
                             sparse=sparse, gc=gc, chunk_size=chunk_size, alpha = alpha)
    
    # Apply NW regression with adaptive bandwidths
    est = NWregression(X, Y, x=x, nEval=nEval, kernel=kernel, D=D, 
                      method.h=method.h, h=h, lambda=lambda,
                      sparse=sparse, gc=gc, chunk_size=chunk_size)
    
    return(est)
}
