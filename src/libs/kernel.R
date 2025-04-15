kernelMethod <- function(X, x=NULL, nEval=2500, kernel.type="gauss", D=NULL, 
                             method.h=NULL, h=NULL, lambda=NULL, 
                             sparse=FALSE, gc=FALSE, chunk_size=nrow(x), type.est=NULL, Y=NULL) {
    
    nObs = nrow(X)
    covX = cov(X)
    invS = solve(covX)
    sqrtinvS = expm::sqrtm(solve(covX))
    detS = det(covX)
    
    Z = X %*% sqrtinvS
    z = x %*% sqrtinvS
    
    if (is.null(x)) { x = defineEvalPoints(X,nEval) }
    nEval =  nrow(x)
    if (DEBUG) {
        print(paste("nEval: ",nEval))
        print(paste("nObs: ",nObs)) }
    if (is.null(chunk_size)) { chunk_size = nrow(x) }
    if (is.null(lambda)) { lambda = rep(1,nObs) }
    
    if (is.null(type.est)) { stop("type.est not specified") }
    if (is.null(h) | is.null(method.h)) { 
        list.h = define_h_method.h(Z,h,method.h, kernel.type)
        h = list.h$h
        method.h = list.h$method.h }
    kernelFunction = defineKernel(kernel.type)
    
    
    density = numeric(nEval)
    kernel_sum = numeric(nEval)
    estimator = numeric(nEval)


    chunks = split(seq_len(nEval), ceiling(seq_len(nEval)/chunk_size))
    
    # Start computing the trace of the H matrix needed for AICc if LL is used
    if (type.est == "LL") {
        partialTraceHLL = numeric(length(chunks))
    }
    
    if (DEBUG) {
        print(paste("Computing ",length(chunks)," chunks"))
        print(paste("Chunk size: ",chunk_size)) }
    
    # start estimate
    for(i in 1:length(chunks)) {
        if (DEBUG) print(paste("Computing chunk ",i,"/",length(chunks),sep=""))
        
        chunk = chunks[[i]]
        z_chunk = z[chunk, ,drop=FALSE]
        
        if (is.null(D)){ D_chunk = computeDcomponents(Z, z_chunk, sparse=sparse) 
        } else { D_chunk = list(z1=D$z1[,chunk], z2=D$z2[,chunk]) }
        
        # Kernel computation
        K = kernelFunction(sweep(D_chunk$z1, 1, h * lambda, "/"),sweep(D_chunk$z2, 1, h * lambda, "/"))
        if (gc == TRUE){ gc() }
        
        K_scaled = sweep(K, 1, lambda^2, "/")
        
        # insert bootstrap here 
        
        # switch between density, NW
        switch(type.est,
               "density" = {
                   estimator[chunk] = computeTerms(D_chunk, Y, h, detS, K_scaled, type.est) 
                   density[chunk] = estimator[chunk]},
               "NW" = {
                   estimator[chunk] = computeTerms(D_chunk, Y, h, detS, K_scaled, type.est)
                   density[chunk] = computeTerms(D_chunk, Y, h, detS, K_scaled, "density")  },
               "LL" = {
                   LLoutputs = computeTerms(D_chunk, Y, h, detS, K_scaled, type.est) 
                   estimator[chunk] = LLoutputs$solutions
                   print("culo")
                   partialTraceHLL[i] = LLoutputs$partialTraceHLL
                   density[chunk] = computeTerms(D_chunk, Y, h, detS, K_scaled, "density")},
               {
                   stop(paste("Invalid type.est:", type.est, ". Must be either 'density' or 'NW'."))
               }
        )
        if (gc == TRUE){ gc() }
    }
    if (type.est == "LL") {
        return(listN(x, estimator, density, h, method.h, kernel.type, partialTraceHLL))
    } else {
        return(listN(x, estimator, density, h, method.h, kernel.type))
    }
}



computeTerms <- function(distances, Y, h, detS, K_scaled, type.est){
    nObs = dim(distances$z1)[1]
    nEval = dim(distances$z1)[2]
    d1 = distances$z1
    d2 = distances$z2

    switch(type.est,
        "density" = {
            return(1/(nObs * h^2 * sqrt(detS)) * colSums(K_scaled))
        },
        "NW" = {
            K_scaledY = sweep(K_scaled, 1, Y, "*")
            numerator = colSums(K_scaledY)
            denominator = colSums(K_scaled)
            return(numerator/denominator)
        },
        "LL" = {
            K_scaledY = sweep(K_scaled, 1, Y, "*")
            K_scaled_d1 = K_scaled * d1
            K_scaled_d2 = K_scaled * d2
            
            # to compute tr(H)
            if (nEval <= nObs){
                tau = array(c(K_scaled,K_scaled_d1,K_scaled_d2),c(nObs,nEval,3))
            }
            
            S0 = colSums(K_scaled)
            S10 = colSums(K_scaled_d1)
            S20 = colSums(K_scaled_d2)
            S12 = colSums(K_scaled_d1 * d2)
            S11 = colSums(K_scaled * d1^2)
            S22 = colSums(K_scaled * d2^2)
            T0 = colSums(K_scaledY)
            T1 = colSums(K_scaledY * d1)
            T2 = colSums(K_scaledY * d2)
        
            MNumerator = array(NA,dim=c(3,3,nEval))
            MNumerator[1,1,] = T0
            MNumerator[1,2,] = S10
            MNumerator[1,3,] = S20
            MNumerator[2,1,] = T1
            MNumerator[2,2,] = S11
            MNumerator[2,3,] = S12
            MNumerator[3,1,] = T2
            MNumerator[3,2,] = S12
            MNumerator[3,3,] = S22
            
            MDenominator = MNumerator
            MDenominator[1,1,] = S0
            MDenominator[2,1,] = S10
            MDenominator[3,1,] = S20
            
            bConstant = rbind(T0,T1,T2)

            solve_system <- function(i) {
                # if ( cond(MDenominator[,,i]) > 1/.Machine$double.eps ) {
                if ( det(MDenominator[,,i]) == 0 ) {
                    return(rep(NaN,3))
                }
                return(as.numeric(ginv(MDenominator[,,i]) %*% bConstant[,i]))
                # return(as.numeric(solve(MDenominator[,,i]) %*% bConstant[,i]))
            }
            
            # take only the first component of each element of the list 
            solutions <- matrix(unlist(purrr::map(1:nEval, ~ solve_system(.x))),nrow=3,byrow = FALSE)[1,]
            
            

            # to compute tr(H)
            if (nEval <= nObs){
                
                solve_trH <- function(i){
                    if ( det(MDenominator[,,i]) == 0 ) {
                        return(rep(NaN,3))
                    }
                    return((ginv(MDenominator[,,i]) %*% t(tau[,i,]))[1,i])
                }
                partialTraceHLL <- sum( unlist(purrr::map(1:nEval, ~ solve_trH(.x))) ,na.rm=T)
                
                # partialTraceHLL = ((S11 * S22 - S12^2)/(S0 * (S11 * S22 - S12^2)
                #                                         - S10 * (S10 * S22 - S12 * S20)
                #                                         + S20 * (S20 * S12 - S11 * S20)))

            } else { partialTraceHLL = NULL }
            
            return(listN(solutions, partialTraceHLL))
        }
    )

}

defineKernel <- function(kernel.type) {
    if (is.null(kernel.type) || kernel.type == "epa"){ return(epaKernel) }
    else if (kernel.type == "gauss") { return(gaussKernel) } 
    else { stop("kernel not recognized") }
}

epaKernel <- function(z1,z2){
    epaKernel1d(z1)*epaKernel1d(z2)
}

epaKernel1d <- function(z){
    3/4*(1-z^2)*(abs(z) <= 1)
}

gaussKernel <- function(z1,z2){
    gaussKernel1d(z1)*gaussKernel1d(z2)
}

gaussKernel1d <- function(z){
    1/sqrt(2*pi)*exp(-z^2/2)
}

