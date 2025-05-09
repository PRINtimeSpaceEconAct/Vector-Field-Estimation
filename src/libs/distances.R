require(Matrix)
require(spatstat)
 
#' Computes distance components between data points and evaluation points
#' 
#' @param X Matrix of data points (nObs x 2)
#' @param x Matrix of evaluation points (nEval x 2)
#' @param A Covariance matrix of X (not used)
#' @param dMax Maximum distance to consider, 0 otherwise (not used)
#' @param sparse Whether to return sparse matrices (default: TRUE)
#' 
#' @return A list containing two matrices of differences (nObs x nEval) as regular or sparse matrices:
#'         \item{z1}{Differences in the first dimension}
#'         \item{z2}{Differences in the second dimension}
computeDcomponents <- function(X,x,A=NULL,dMax=NULL,sparse=TRUE){
    # X = matrix of data points (nObs x 2)
    # x = matrix of evaluation points (nEval x 2)
    # dMax = maximum distance to consider, 0 otherwise (not used)
    # A matrix of covariance of X (not used)
    # returns two matrices of differences (nObs x nEval) as SparseMatrix 
    

        # test  
    # X = matrix(nrow=10,runif(20))
    # x = matrix(nrow=5,runif(10))

    nObs = nrow(X)
    nEval = nrow(x)
    MObs1 = matrix(rep(X[,1],nEval),nObs,nEval,byrow=FALSE)
    MObs2 = matrix(rep(X[,2],nEval),nObs,nEval,byrow=FALSE)
    MEval1 = matrix(rep(x[,1],nObs),nObs,nEval,byrow=TRUE)
    MEval2 = matrix(rep(x[,2],nObs),nObs,nEval,byrow=TRUE)
  
    if (sparse == TRUE){
        z1 = as(MEval1 - MObs1,"sparseMatrix")
        z2 = as(MEval2 - MObs2,"sparseMatrix")
    }else{
        z1 = MEval1 - MObs1
        z2 = MEval2 - MObs2
    }

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
    

