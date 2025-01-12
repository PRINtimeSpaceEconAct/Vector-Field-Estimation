densityEst2d <- function(X, x=NULL, nEval=2500, kernel="epa", D=NULL, 
                        method.h=NULL, h=NULL, lambda = NULL, 
                        sparse=FALSE, gc=FALSE, chunk_size=1024) {
    
    nObs = nrow(X)
    covX = cov(X)
    invS = solve(covX)
    detS = det(covX)

    if (is.null(x)){
        xGrid = seq(from=min(X[,1]), to=max(X[,1]), length.out=round(sqrt(nEval)))
        yGrid = seq(from=min(X[,2]), to=max(X[,2]), length.out=round(sqrt(nEval)))
        x = as.matrix(expand.grid(xGrid,yGrid))
    }
    nEval = nrow(x)

    if (is.null(lambda)){
        lambda = rep(1,nObs)
    }
    
    # If h is provided directly, use it regardless of method.h
    if (!is.null(h)) {
        # Keep h as is
    } else if (is.null(method.h)) {
        # Default method when both h and method.h are null
        method.h = "silverman"
        h =  nrow(X)^(-1/6)
    } else {
        # Select h based on specified method
        h = switch(method.h,
            "silverman" = nrow(X)^(-1/6),
            "sj" = bw.SJ(X),
            "botev" = botev(X),
            stop("method.h not recognized")
        )
    }
    
    print(paste("Using h = ", h, "and method = ", method.h))

    # Initialize empty vector for density estimates
    densityEst = numeric(nEval)
    
    # Process x in chunks
    chunks = split(seq_len(nEval), ceiling(seq_len(nEval)/chunk_size))

    for(chunk in chunks) {
        x_chunk = x[chunk, , drop=FALSE]
        if (is.null(D)){ 
            D_chunk = computeDcomponents(X, x_chunk, sparse=sparse)
            if (gc == TRUE){ gc() }
        }
        
        M = mahalanobis(D_chunk$z1, D_chunk$z2, A=invS, den=h^2 * lambda^2)
        if (gc == TRUE){ gc() }
        
        # Kernel computation
        if (is.null(kernel) || kernel == "epa"){ 
            K = epaKernel(M) 
        } else if (kernel == "gauss"){ 
            K = gaussKernel(M) 
        } else { 
            stop("kernel not recognized") 
        }
        
        densityEst[chunk] = 1/(nObs * h^2 * sqrt(detS)) * colSums(sweep(K,1,lambda^2,"/"))
        if (gc == TRUE){ gc() }
    }
    
    return(listN(x, densityEst))
}

getLocalBandwidth <- function(X, kernel="epa", D=NULL, method.h=NULL, h=NULL,
                            sparse=FALSE, gc=FALSE, chunk_size=1024, alpha = 0.5){
    nObs = nrow(X)
    nEval = nObs

    pilotDensity = densityEst2d(X, x=X, nEval=nEval, kernel=kernel,
                                D=D, method.h=method.h, h=h, sparse=sparse,
                                gc=gc, chunk_size=chunk_size)
    g = exp(mean(pilotDensity$densityEst))
    lambda = (pilotDensity$densityEst/g)^(-alpha)

    return(lambda)
}

densityEst2dAdaptive <- function(X, x=NULL, nEval=2500, kernel="epa", D=NULL,
                                method.h=NULL, h=NULL,
                                sparse=FALSE, gc=FALSE, chunk_size=1024, alpha = 0.5){

    lambda = getLocalBandwidth(X, kernel=kernel, D=D, method.h=method.h, h=h,
                                sparse=sparse, gc=gc, chunk_size=chunk_size, alpha = alpha)

    est = densityEst2d(X, x=x, nEval=nEval, kernel=kernel, D=D, method.h=method.h, h=h,
                                lambda=lambda, sparse=sparse, gc=gc, chunk_size=chunk_size)
    return(est)
}