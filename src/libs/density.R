


densityEst2d <- function(X,x=NULL,nEval=2500,
             kernel="epa",D=NULL,method.h=NULL,h=NULL,sparse=FALSE,gc=FALSE){
    
    nObs = nrow(X)
    covX = cov(X)
    invS = solve(covX)
    detS = det(covX)

    if (is.null(x)){
        xGrid = seq(from=min(X[,1]),to=max(X[,1]),length.out=round(sqrt(nEval)))
        yGrid = seq(from=min(X[,2]),to=max(X[,2]),length.out=round(sqrt(nEval)))
        x = as.matrix(expand.grid(xGrid,yGrid))
    }
    nEval = nrow(x)
    
    if (is.null(D)){ 
        D = computeDcomponents(X,x,sparse=sparse); if (gc == TRUE){ gc() }
    }
    
    # defaults to Sheather & Jones (1991) 
    # if (is.null(h) & is.null(method.h)){ h = 1.77*nrow(X)^(-1/6) } 
    # else if (method.h == "silverman") { h = 1.77*nrow(X)^(-1/6) }
    # else if (method.h == "sj") { h = bw.SJ(X) }
    # else if (method.h == "scott") { h = bw.scott(X) }
    # else if (method.h == "botev") { h = botev(X) }
    # else { stop("method.h not recognized") }
    
    
    M = mahalanobis(D$z1,D$z2,A=invS,den=rep(h^2,nObs)); if (gc == TRUE){ gc() }
    
    # defaults to Epanechnikov kernel
    if (is.null(kernel)){  K = epaKernel(M) }
    else if (kernel == "epa"){ K = epaKernel(M) }
    else if (kernel == "gauss"){ K = gaussKernel(M) }
    else { stop("kernel not recognized") }
    
    densityEst = 1/(nObs *h^2 * sqrt(detS)) * colSums(K)
    
    return(listN(x,densityEst))
}