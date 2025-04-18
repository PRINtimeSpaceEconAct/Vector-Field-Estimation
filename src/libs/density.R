densityEst2d <- function(X, x=NULL, nEval=2500, kernel.type="gauss", D=NULL, 
                          method.h=NULL, h=NULL, lambda=NULL, 
                          sparse=FALSE, gc=FALSE, chunk_size=nrow(x)) {
    
    type.est = "density"
    resultEst = kernelMethod(X=X,x=x,nEval=nEval,kernel.type=kernel.type,D=D,
                             method.h=method.h,h=h,lambda=lambda,
                             sparse=sparse,gc=gc,chunk_size=chunk_size,type.est=type.est)
    return(resultEst)

}
    
densityEst2dAdaptive <- function(X, x=NULL, nEval=2500, kernel.type="gauss", D=NULL,
                                method.h=NULL, h=NULL,
                                sparse=FALSE, gc=FALSE, chunk_size=nrow(X), alpha = 0.5){

    lambda = getLocalBandwidth(X, kernel.type=kernel.type, D=D, method.h=method.h, h=h,
                                sparse=sparse, gc=gc, chunk_size=chunk_size, alpha = alpha)
    est = densityEst2d(X, x=x, nEval=nEval, kernel.type=kernel.type, D=D, method.h=method.h,
                        h=h, lambda=lambda, sparse=sparse, gc=gc, chunk_size=chunk_size)
    return(est)
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



