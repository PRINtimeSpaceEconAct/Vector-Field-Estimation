normKernelBiv_vec <- function(y,x,
                              lambdas=rep(1,dim(x)[1]),
                              invS=solve(cov(x)),
                              detS=det(cov(x)),
                              weiNorm=rep(1/dim(x)[1],dim(x)[1]),
                              h.bandwidth=c(0.96*nrow(x)^(-1/6),0.96*nrow(x)^(-1/6))){
    
    #y: vector of NEval x 2
    #x: matrix with observation of dimensions NObs x 2
    # weiNorm: weights of observations
    
  
    normPDF <- function(z){
        1/(2*pi)*exp(-z/2)
    }
    
    NObs = nrow(x)
    NEval = nrow(y)
    MEval1 = matrix(rep(y[,1],NObs),NEval,NObs)
    MEval2 = matrix(rep(y[,2],NObs),NEval,NObs)
    MObs1 = matrix(rep(x[,1],NEval),NEval,NObs,byrow=TRUE)
    MObs2 = matrix(rep(x[,2],NEval),NEval,NObs,byrow=TRUE)
    MLambdas = matrix(rep(lambdas,NEval),NEval,NObs,byrow=TRUE)
    
    M1 = MEval1 - MObs1
    M2 = MEval2 - MObs2 
    
    # quadratic form 
    QF = invS[1,1]*M1^2 + invS[1,2]*M1*M2 + invS[2,1]*M1*M2 + invS[2,2]*M2^2
    
    # argument of standard gaussian kernel 
    MArg = (1/prod(h.bandwidth))*QF/MLambdas^2
    MK =  normPDF(MArg)
    MWei = matrix(rep(weiNorm,NEval),NEval,NObs,byrow=TRUE)
    fx = (detS)^(-1/2)/(prod(h.bandwidth))*rowSums(MWei*MK/MLambdas^2)
    
    return(fx)    
    
}

normKernelBivAdaptive <- function(evalPoints_Stack,x,
                                 invS=solve(cov(x)),
                                 detS=det(cov(x)),
                                 weiNorm=rep(1/nrow(x),nrow(x)),
                                 alpha=0.5, ngrid,
                                 h.bandwidth=c(0.96*nrow(x)^(-1/6),0.96*nrow(x)^(-1/6)),
                                 adaptive=TRUE){
    
    ###################################
    ## Calculation of lambdas from the density pilot
    
    ##Pilot estimate, Silverman (1986), p. 101
    # pb = txtProgressBar(min = 1, max = nrow(x), initial = 0) 
    # densityPilot <- foreach(i=1:nrow(x), .combine=c, .verbose=FALSE) %do%{  
    #     setTxtProgressBar(pb,i)
    #     normKernelBiv(y=x[i,],x=x,invS=invS,detS=detS,weiNorm=weiNorm,lambdas=rep(1,dim(x)[1]), 
    #                   h.bandwidth=h.bandwidth)
    #     # print(paste("i = ",i,"/",nrow(x)))
    # }
    # close(pb)
    
    if (adaptive==TRUE){
        print("Start pilot estimate")
        densityPilot <- normKernelBiv_vec(y=x,x=x,invS=invS,detS=detS,
                                          weiNorm=weiNorm,lambdas=rep(1,dim(x)[1]))
        print("end estimate pilot")
        
        # correct 
        lambdas <- (densityPilot/exp(sum(weiNorm*log(densityPilot))))^(-alpha)
    }else{
        lambdas <- rep(1,dim(x)[1])
    }
    
    print("Start estimate of the stochastic kernel")
    estimateDensity_Stack <- normKernelBiv_vec(y=evalPoints_Stack,x=x,
                                          invS=invS,detS=detS,
                                          weiNorm=weiNorm,lambdas=lambdas,
                                          h.bandwidth=h.bandwidth)
    print("End estimate of the stochastic kernel")
    return(estimateDensity_Stack) # nobs x (ngrid*ngrid) 
}                    