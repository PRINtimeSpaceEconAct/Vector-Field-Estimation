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
    


