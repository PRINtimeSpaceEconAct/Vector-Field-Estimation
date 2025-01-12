require(Matrix)
require(spatstat)

computeD <- function(X,x,dMax){
    # X = matrix of data points (nObs x 2)
    # x = matrix of evaluation points (nEval x 2)
    # dMax = maximum distance to consider, 0 otherwise
    # returns a matrix of distances (nObs x nEval) as SparseMatrix
    
        # test
    # X = matrix(nrow=10,runif(20)) 
    # x = matrix(nrow=5,runif(10)) 
    # dMax = 0.5
    D = crossdist.default(X[,1],X[,2],x[,1],x[,2],period=NULL)
    
    D[D >= dMax] = 0
    D = as(D,"sparseMatrix")
    return(D)
    
}
    
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

mahalanobis <- function(z1,z2,A,den){
    # z1,z2 = matrices of differences on both components (nObs x nEval)
    # A matrix of covariance of observations (2 x 2)
    # den vector of length nObs for denominators

    nObs = nrow(z1)
    nEval = ncol(z1)

    QF = A[1,1]*z1^2 + (A[1,2] + A[2,1])*z1*z2 + A[2,2]*z2^2

    Mmahalanobis = sweep(QF,1,den,"/")

    return(Mmahalanobis)
}
    

