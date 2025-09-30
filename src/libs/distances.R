require(Matrix)
require(spatstat)
 
#' Computes distance components between data points and evaluation points
#' 
#' @param X Matrix of data points (nObs x 2)
#' @param x Matrix of evaluation points (nEval x 2)#' 
#' @return A list containing two matrices of differences (nObs x nEval) as matrices:
#'         \item{z1}{Differences in the first dimension}
#'         \item{z2}{Differences in the second dimension}
computeDcomponents <- function(X,x){
    
    nObs = nrow(X)
    nEval = nrow(x)
    MObs1 = matrix(rep(X[,1],nEval),nObs,nEval,byrow=FALSE)
    MObs2 = matrix(rep(X[,2],nEval),nObs,nEval,byrow=FALSE)
    MEval1 = matrix(rep(x[,1],nObs),nObs,nEval,byrow=TRUE)
    MEval2 = matrix(rep(x[,2],nObs),nObs,nEval,byrow=TRUE)
  

    z1 = MObs1 - MEval1
    z2 = MObs2 - MEval2

    return(listN(z1,z2))
}


#' Computes Mahalanobis distances between data points and evaluation points
#' 
#' @param z1 Matrix of differences in first dimension (nObs x nEval)
#' @param z2 Matrix of differences in second dimension (nObs x nEval)
#' @param A Inverse of covariance matrix of observations (2 x 2)
#' @param den Vector of length nObs for denominators
#' 
#' @return Matrix of Mahalanobis distances (nObs x nEval)
mahalanobis <- function(z1,z2,A,den){
    # z1,z2 = matrices of differences on both components (nObs x nEval)
    # A inverse of covariance matrix of observations (2 x 2)
    # den vector of length nObs for denominators

    nObs = nrow(z1)
    nEval = ncol(z1)

    QF = A[1,1]*z1^2 + (A[1,2] + A[2,1])*z1*z2 + A[2,2]*z2^2

    Mmahalanobis = sweep(QF,1,den,"/")

    return(Mmahalanobis)
}
    

#' Computes distance components between data points and evaluation points for 1D data
#'
#' @param X Vector of data points (nObs)
#' @param x Vector of evaluation points (nEval)
#' @return A matrix of differences (nObs x nEval).
computeDcomponents1d <- function(X,x){

    nObs = length(X)
    nEval = length(x)
    
    MObs = matrix(rep(X, nEval), nObs, nEval, byrow=FALSE)
    MEval = matrix(rep(x, nObs), nObs, nEval, byrow=TRUE)

    return(MEval - MObs)
}
    

