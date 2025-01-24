kernelMethod <- function(X, x=NULL, nEval=2500, kernel.type="gauss", D=NULL, 
                             method.h=NULL, h=NULL, lambda=NULL, 
                             sparse=FALSE, gc=FALSE, chunk_size=nrow(x), type.est=NULL, Y=NULL) {
    
    nObs = nrow(X)
    covX = cov(X)
    invS = solve(covX)
    sqrtinvS = sqrtm(solve(covX))
    detS = det(covX)
    
    if (is.null(x)) { x = defineEvalPoints(X,nEval) }
    nEval =  nrow(x)
    if (DEBUG) {
        print(paste("nEval: ",nEval))
        print(paste("nObs: ",nObs)) }
    if (is.null(chunk_size)) { chunk_size = nrow(x) }
    if (is.null(lambda)) { lambda = rep(1,nObs) }
    
    if (is.null(type.est)) { stop("type.est not specified") }
    if (is.null(h) | is.null(method.h)) { 
        list.h = define_h_method.h(X,h,method.h, kernel.type)
        h = list.h$h
        method.h = list.h$method.h }
    kernelFunction = defineKernel(kernel.type)
    
    Z = X %*% sqrtinvS
    z = x %*% sqrtinvS
    
    estimator = numeric(nEval)
    chunks = split(seq_len(nEval), ceiling(seq_len(nEval)/chunk_size))
    if (DEBUG) {
        print(paste("Computing ",length(chunks)," chunks"))
        print(paste("Chunk size: ",chunk_size)) }
    
    # start estimate
    for(i in 1:length(chunks)) {
        if (DEBUG) print(paste("Computing chunk ",i,"/",length(chunks),sep=""))
        
        chunk = chunks[[i]]
        z_chunk = z[chunk, , drop=FALSE]
        
        D_chunk = computeDcomponents(Z, z_chunk, sparse=sparse) 
        M = mahalanobis(D_chunk$z1, D_chunk$z2, A=diag(2), den=h^2 * lambda^2)
        # Kernel computation
        K = kernelFunction(M)
        if (gc == TRUE){ gc() }
        
        K_scaled = sweep(K, 1, lambda^2, "/")
        
        # insert bootstrap here 
        
        # switch between density, NW
        switch(type.est,
               "density" = {
                   estimator[chunk] = computeTerms(D_chunk, Y, h, detS, K_scaled, type.est) },
               "NW" = {
                   estimator[chunk] = computeTerms(D_chunk, Y, h, detS, K_scaled, type.est) },
               "LL" = {
                   estimator[chunk] = computeTerms(D_chunk, Y, h, detS, K_scaled, type.est) },
               {
                   stop(paste("Invalid type.est:", type.est, ". Must be either 'density' or 'NW'."))
               }
        )
        if (gc == TRUE){ gc() }
    }
    
    return(listN(x, estimator))
}



computeTerms <- function(distances, Y, h, detS, K_scaled, type.est){
    nObs = dim(distances$z1)[1]
    d1 = distances$z1
    d2 = distances$z2

    switch(type.est,
        "density" = {
            return(1/(nObs * h^2 * sqrt(detS)) * colSums(K_scaled))
        },
        "NW" = {
            K_scaledY = sweep(K_scaled, 1, Y, "*")
            numerator = colSums(K_scaledY)
            denominator = colSums(K_scaled)
            return(numerator/denominator)
        },
        "LL" = {
            K_scaledY = sweep(K_scaled, 1, Y, "*")

            S0 = colSums(K_scaled)
            S10 = colSums(K_scaled * d1)
            S20 = colSums(K_scaled * d2)
            S12 = colSums(K_scaled * d1 * d2)
            S11 = colSums(K_scaled * d1^2)
            S22 = colSums(K_scaled * d2^2)
            T0 = colSums(K_scaledY)
            T1 = colSums(K_scaledY * d1)
            T2 = colSums(K_scaledY * d2)
            
            numerator = determinant3(T0, S10, S20, T1, S11, S12, T2, S12, S22)
            denominator = determinant3(S0, S10, S20, S10, S11, S12, S20, S12, S22)
            return(numerator/denominator)
        }
    )

}

defineEvalPoints <- function(X,nEval){
    xGrid = seq(from=min(X[,1]), to=max(X[,1]), length.out=round(sqrt(nEval)))
    yGrid = seq(from=min(X[,2]), to=max(X[,2]), length.out=round(sqrt(nEval)))
    x = as.matrix(expand.grid(xGrid,yGrid))
    return(x)
}

define_h_method.h <- function(X, h, method.h, kernel.type="gauss") {
    # Check for conflicting bandwidth specifications
    if (!is.null(h) && !is.null(method.h)) {
        stop("Cannot specify both h and method.h. Please provide only one.")
    }
    
    # Define constant based on kernel type
    c = if (kernel.type == "epa") 1.77 else 0.96
    
    # If h is provided directly, use it regardless of method.h
    if (!is.null(h)) {
        # Keep h as is
    } else if (is.null(method.h)) {
        # Default method when both h and method.h are null
        method.h = "silverman"
        h = c * nrow(X)^(-1/6)
    } else {
        # Select h based on specified method
        h = switch(method.h,
                   "silverman" = c * nrow(X)^(-1/6),
                   "sj" = bw.SJ(X),
                   "botev" = botev(X),
                   stop("method.h not recognized")
        )
    }
    
    if (DEBUG) print(paste("Using h = ", h, "and method = ", method.h))
    
    return(listN(h, method.h))
}

defineKernel <- function(kernel.type) {
    if (is.null(kernel.type) || kernel.type == "epa"){ return(epaKernel) }
    else if (kernel.type == "gauss") { return(gaussKernel) } 
    else { stop("kernel not recognized") }
}

getLocalBandwidth <- function(X, kernel.type="gauss", D=NULL, method.h=NULL, h=NULL,
                              sparse=FALSE, gc=FALSE, chunk_size=1024, alpha = 0.5){
    nObs = nrow(X)
    nEval = nObs
    
    pilotDensity = densityEst2d(X, x=X, nEval=nEval, kernel.type=kernel.type,
                                D=D, method.h=method.h, h=h, sparse=sparse,
                                gc=gc, chunk_size=chunk_size)
    g = exp(mean(log(pilotDensity$estimator)))
    lambda = (pilotDensity$estimator/g)^(-alpha)
    
    return(lambda)
}


epaKernel <- function(z){
    # to be computed on z S^-1 z
    2/pi*(1-z)*(z <= 1)
}

gaussKernel <- function(z){
    # to be computed on z S^-1 z
    1/(2*pi)*exp(-z/2)
}

