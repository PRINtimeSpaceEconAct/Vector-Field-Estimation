#' Performs 2D kernel density estimation
#' 
#' @param X Matrix of input points (nObs x 2)
#' @param x Matrix of evaluation points (nEval x 2, if NULL, generated internally)
#' @param nEval Number of evaluation points if x is NULL (default: 2500)
#' @param kernel.type Type of kernel function to use (default: "gauss")
#' @param D Pre-computed distance components (if NULL, computed internally)
#' @param method.h Method for bandwidth selection (if NULL, specified h is used)
#' @param h Bandwidth parameter (if NULL, selected by method.h)
#' @param lambda Vector of local bandwidths for adaptive estimation (nEval)
#' @param gc Whether to force garbage collection (default: FALSE)
#' @param chunk_size Number of points to process at once (default: nrow(x))
#' 
#' @return A list containing density estimation results and parameters:
#'   \item{estimator}{Vector of density estimates at evaluation points (nEval)}
#'   \item{x}{Matrix of evaluation points (nEval x 2)}
#'   \item{h}{Bandwidth parameter used}
#'   \item{method.h}{Method used for bandwidth selection}
#'   \item{kernel.type}{Type of kernel function used}
#'   \item{lambda}{Vector of local bandwidths used, if applicable (nObs)}
densityEst2d <- function(X, x=NULL, nEval=2500, kernel.type="gauss", D=NULL, 
                          method.h=NULL, h=NULL, lambda=NULL, 
                          gc=FALSE, chunk_size=nrow(x)) {
    
    type.est = "density"
    resultEst = kernelMethod(X=X,x=x,nEval=nEval,kernel.type=kernel.type,D=D,
                             method.h=method.h,h=h,lambda=lambda,
                             gc=gc,chunk_size=chunk_size,type.est=type.est)
    return(resultEst)

}
    
#' Performs 2D kernel density estimation with adaptive bandwidth
#' 
#' @param X Matrix of input points (nObs x 2)
#' @param x Matrix of evaluation points (nEval x 2, if NULL, generated internally)
#' @param nEval Number of evaluation points if x is NULL (default: 2500)
#' @param kernel.type Type of kernel function to use (default: "gauss")
#' @param D Pre-computed distance components (if NULL, computed internally)
#' @param method.h Method for bandwidth selection (if NULL, specified h is used)
#' @param h Bandwidth parameter (if NULL, selected by method.h)
#' @param gc Whether to force garbage collection (default: FALSE)
#' @param chunk_size Number of points to process at once (default: nrow(X))
#' @param alpha Sensitivity parameter for adaptive bandwidth (default: 0.5)
#' 
#' @return A list containing density estimation results with adaptive bandwidths:
#'   \item{estimator}{Vector of density estimates at evaluation points (nEval)}
#'   \item{x}{Matrix of evaluation points (nEval x 2)}
#'   \item{h}{Bandwidth parameter used}
#'   \item{method.h}{Method used for bandwidth selection}
#'   \item{kernel.type}{Type of kernel function used} 
#'   \item{lambda}{Vector of local bandwidths used (nObs)}
densityEst2dAdaptive <- function(X, x=NULL, nEval=2500, kernel.type="gauss", D=NULL,
                                method.h=NULL, h=NULL,
                                gc=FALSE, chunk_size=nrow(X), alpha = 0.5){

    lambda = getLocalBandwidth(X, kernel.type=kernel.type, D=D, method.h=method.h, h=h,
                                gc=gc, chunk_size=chunk_size, alpha = alpha)
    est = densityEst2d(X, x=x, nEval=nEval, kernel.type=kernel.type, D=D, method.h=method.h,
                        h=h, lambda=lambda, gc=gc, chunk_size=chunk_size)
    return(est)
}


#' Calculates local bandwidths for adaptive kernel estimation based on pilot density
#' 
#' @param X Matrix of input points (nObs x 2)
#' @param kernel.type Type of kernel function to use (default: "gauss")
#' @param D Pre-computed distance components (if NULL, computed internally)
#' @param method.h Method for bandwidth selection (if NULL, specified h is used)
#' @param h Bandwidth parameter (if NULL, selected by method.h)
#' @param gc Whether to force garbage collection (default: FALSE)
#' @param chunk_size Number of points to process at once (default: 1024)
#' @param alpha Sensitivity parameter for adaptive bandwidth (default: 0.5)
#' 
#' @return Vector of local bandwidths for each input point (nObs)
getLocalBandwidth <- function(X, kernel.type="gauss", D=NULL, method.h=NULL, h=NULL,
                              gc=FALSE, chunk_size=1024, alpha = 0.5){
    nObs = nrow(X)
    nEval = nObs
    
    pilotDensity = densityEst2d(X, x=X, nEval=nEval, kernel.type=kernel.type,
                                D=D, method.h=method.h, h=h,
                                gc=gc, chunk_size=chunk_size)
    g = exp(mean(log(pilotDensity$estimator)))
    lambda = (pilotDensity$estimator/g)^(-alpha)
    
    return(lambda)
}



