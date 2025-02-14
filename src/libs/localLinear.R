LLregression <- function(X, Y, x=NULL, nEval=2500, kernel.type="gauss", D=NULL, 
                        method.h=NULL, h=NULL, lambda = NULL, 
                        sparse=FALSE, gc=FALSE, chunk_size=nrow(x)) {

    type.est = "LL"
    resultEst = kernelMethod(X=X,x=x,nEval=nEval,kernel.type=kernel.type,D=D,
                             method.h=method.h,h=h,lambda=lambda,
                             sparse=sparse,gc=gc,chunk_size=chunk_size,type.est=type.est, Y=Y)
    
    return(resultEst)
}

# Adaptive bandwidth version of LL regression
LLregressionAdaptive <- function(X, Y, x=NULL, nEval=2500, kernel.type="gauss", D=NULL,
                                method.h=NULL, h=NULL, lambda=NULL,
                                sparse=FALSE, gc=FALSE, chunk_size=nrow(x), alpha = 0.5) {
    
    # Get adaptive bandwidths using the density-based approach
    lambda = getLocalBandwidth(X, kernel.type=kernel.type, D=D, method.h=method.h, h=h,
                             sparse=sparse, gc=gc, chunk_size=chunk_size, alpha = alpha)
    
    # Apply LL regression with adaptive bandwidths
    est = LLregression(X, Y, x=x, nEval=nEval, kernel.type=kernel.type, D=D, 
                      method.h=method.h, h=h, lambda=lambda,
                      sparse=sparse, gc=gc, chunk_size=chunk_size)
    est$lambda = lambda
    return(est)
}

LLfield <- function(X0, X1, x=NULL, nEval=2500, kernel.type="gauss", D=NULL,
                   method.h=NULL, h=NULL, lambda=NULL,
                   sparse=FALSE, gc=FALSE, chunk_size=nrow(x)) {
    
    Y1 = X1[,1] - X0[,1]
    Y2 = X1[,2] - X0[,2]
    est1 = LLregression(X0, Y1, x=x, nEval=nEval, kernel.type=kernel.type, D=D, 
                       method.h=method.h, h=h, lambda=lambda,
                       sparse=sparse, gc=gc, chunk_size=chunk_size)
    est2 = LLregression(X0, Y2, x=est1$x, nEval=nEval, kernel.type=kernel.type, D=D, 
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
    type.est = "LL"
    
    return(listN(x, X0, X1, estimator, type.est, density, kernel.type, h, method.h, lambda))
}

LLfieldAdaptive <- function(X0, X1, x=NULL, nEval=2500, kernel.type="gauss", D=NULL,
                           method.h=NULL, h=NULL, lambda=NULL,
                           sparse=FALSE, gc=FALSE, chunk_size=nrow(x), alpha = 0.5) {
    lambda = getLocalBandwidth(X0, kernel.type=kernel.type, D=D, method.h=method.h, h=h,
                             sparse=sparse, gc=gc, chunk_size=chunk_size, alpha = alpha)
    est = LLfield(X0, X1, x=x, nEval=nEval, kernel.type=kernel.type, D=D, 
                 method.h=method.h, h=h, lambda=lambda,
                 sparse=sparse, gc=gc, chunk_size=chunk_size)
    est$lambda = lambda
    return(est)
}