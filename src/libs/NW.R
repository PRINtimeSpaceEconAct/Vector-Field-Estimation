NWregression <- function(X, Y, x=NULL, nEval=2500, kernel.type="gauss", D=NULL, 
                        method.h=NULL, h=NULL, lambda = NULL, 
                        sparse=FALSE, gc=FALSE, chunk_size=nrow(x)) {

    type.est = "NW"
    resultEst = kernelMethod(X=X,x=x,nEval=nEval,kernel.type=kernel.type,D=D,
                             method.h=method.h,h=h,lambda=lambda,
                             sparse=sparse,gc=gc,chunk_size=chunk_size,type.est=type.est, Y=Y)
    
    return(resultEst)
}

# Adaptive bandwidth version of NW regression
NWregressionAdaptive <- function(X, Y, x=NULL, nEval=2500, kernel.type="gauss", D=NULL,
                                method.h=NULL, h=NULL, lambda=NULL,
                                sparse=FALSE, gc=FALSE, chunk_size=nrow(x), alpha = 0.5) {
    
    # Get adaptive bandwidths using the density-based approach
    lambda = getLocalBandwidth(X, kernel.type=kernel.type, D=D, method.h=method.h, h=h,
                             sparse=sparse, gc=gc, chunk_size=chunk_size, alpha = alpha)
    
    # Apply NW regression with adaptive bandwidths
    est = NWregression(X, Y, x=x, nEval=nEval, kernel.type=kernel.type, D=D, 
                      method.h=method.h, h=h, lambda=lambda,
                      sparse=sparse, gc=gc, chunk_size=chunk_size)
    est$lambda = lambda
    return(est)
}

NWfield <- function(X0, X1, x=NULL, nEval=2500, kernel.type="gauss", D=NULL,
                    method.h=NULL, h=NULL, lambda=NULL,
                    sparse=FALSE, gc=FALSE, chunk_size=nrow(x),
                    hOpt = FALSE, nGridh = 10) {

    Y1 = X1[,1] - X0[,1]
    Y2 = X1[,2] - X0[,2]
    
    if (hOpt == TRUE){
        list.h = define_h_method.h(X0, NULL ,"silverman", kernel.type)
        hStart = list.h$h/10
        hEnd = list.h$h*2
        hGrid = exp(log(hStart) + (log(hEnd) - log(hStart)) * (0:(nGridh-1))/(nGridh-1))
        
        RSS = array(NA,dim = nGridh)
        for (i in 1:nGridh){
            if (DEBUG) print(paste("Computing h ",i,"/",nGridh,sep=""))
            hi = hGrid[i]
            # est1 = NWregression(X0, Y1, x=x, nEval=nEval, kernel.type=kernel.type, D=D, 
            #                     method.h=method.h, h=hi, lambda=lambda,
            #                     sparse=sparse, gc=gc, chunk_size=chunk_size)
            # est2 = NWregression(X0, Y2, x=est1$x, nEval=nEval, kernel.type=kernel.type, D=D, 
            #                     method.h=method.h, h=hi, lambda=lambda,
            #                     sparse=sparse, gc=gc, chunk_size=chunk_size)
            # VFhi = list(x = est1$x, estimator = cbind(est1$estimator, est2$estimator))
            # X1Hat = forecastDiscrete(X0,VFhi,speedFactor=1,nPeriods=1)
            
            
            # stima alla cazzo fatta sulle obs con un forecast solo
            est1 = NWregression(X0, Y1, x=X0, nEval=nEval, kernel.type=kernel.type, D=D, 
                                method.h=method.h, h=hi, lambda=lambda,
                                sparse=sparse, gc=gc, chunk_size=chunk_size)
            est2 = NWregression(X0, Y2, x=est1$x, nEval=nEval, kernel.type=kernel.type, D=D, 
                                method.h=method.h, h=hi, lambda=lambda,
                                sparse=sparse, gc=gc, chunk_size=chunk_size)
            
            X1Hat = X0 + cbind(est1$estimator, est2$estimator)
            RSS[i] = sum(rowSums((X1Hat - X1)^2,na.rm=TRUE))
            
            # maxLength = max(sqrt(rowSums(VFhi$estimator)^2),na.rm=T)
            # nPeriods = ceil(10*maxLength)
            # X1Hat = forecastDiscrete(X0,VFhi,speedFactor=1/nPeriods,nPeriods=nPeriods)
            # RSS[i] = sum(rowSums((X1Hat[,,nPeriods] - X1)^2,na.rm=TRUE))
        }

        h = hGrid[which.min(RSS)]
        if (DEBUG) print(hGrid)
        if (DEBUG) print(RSS)
        if (DEBUG) print(paste("Optimal h: ",h))
    }
    
    est1 = NWregression(X0, Y1, x=x, nEval=nEval, kernel.type=kernel.type, D=D, 
                      method.h=method.h, h=h, lambda=lambda,
                      sparse=sparse, gc=gc, chunk_size=chunk_size)
    est2 = NWregression(X0, Y2, x=est1$x, nEval=nEval, kernel.type=kernel.type, D=D, 
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
    type.est = "NW"
    
    return(listN(x, X0, X1, estimator, type.est, density, kernel.type, h, method.h, lambda))
}

NWfieldAdaptive <- function(X0, X1, x=NULL, nEval=2500, kernel.type="gauss", D=NULL,
                            method.h=NULL, h=NULL, lambda=NULL,
                            sparse=FALSE, gc=FALSE, chunk_size=nrow(x), alpha = 0.5) {
    
    lambda = getLocalBandwidth(X0, kernel.type=kernel.type, D=D, method.h=method.h, h=h,
                             sparse=sparse, gc=gc, chunk_size=chunk_size, alpha = alpha)
    est = NWfield(X0, X1, x=x, nEval=nEval, kernel.type=kernel.type, D=D, 
                      method.h=method.h, h=h, lambda=lambda,
                      sparse=sparse, gc=gc, chunk_size=chunk_size)
    est$lambda = lambda
    return(est)
}

