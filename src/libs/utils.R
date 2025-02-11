listN <- function(...){
    # automatically give names to list elements = var name
    anonList <- list(...)
    names(anonList) <- as.character(substitute(list(...)))[-1]
    anonList
}

chunkSize <- function(){
    freeMem = Sys.procmem()$freeram
    return(10000)
}

determinant3 <- function(a00, a01, a02, a10, a11, a12, a20, a21, a22){
    return(a00 * (a11 * a22 - a12 * a21) - a01 * (a10 * a22 - a12 * a20) + a02 * (a10 * a21 - a11 * a20))
}

defineEvalPoints <- function(X,nEval){
    xGrid = seq(from=min(X[,1]), to=max(X[,1]), length.out=round(sqrt(nEval)))
    yGrid = seq(from=min(X[,2]), to=max(X[,2]), length.out=round(sqrt(nEval)))
    x = as.matrix(expand.grid(xGrid,yGrid))
    return(x)
}

define_h_method.h <- function(X, h, method.h, kernel.type="gauss") {
    # Check for conflicting bandwidth specifications
    if (!is.null(h) && !is.null(method.h)) {
        stop("Cannot specify both h and method.h. Please provide only one.")
    }
    
    # Define constant based on kernel type
    c = if (kernel.type == "epa") 1.77 else 0.96
    
    # If h is provided directly, use it regardless of method.h
    if (!is.null(h)) {
        # Keep h as is
    } else if (is.null(method.h)) {
        # Default method when both h and method.h are null
        method.h = "silverman"
        h = c * nrow(X)^(-1/6)
    } else {
        # Select h based on specified method
        h = switch(method.h,
                   "silverman" = c * nrow(X)^(-1/6),
                   "sj" = bw.SJ(X),
                   "botev" = botev(X),
                   stop("method.h not recognized")
        )
    }
    
    if (DEBUG) print(paste("Using h = ", h, "and method = ", method.h))
    
    return(listN(h, method.h))
}

interp2d <- function(z,x,VF){
    # z = where to interpolate <- (nInterp x 2)
    # x = evaluation points where f is computed <- (nEval x 2)
    # VF = Vector Field to interpolate <- (nEval x 2)
    
    VFz1 = interp1d(z,x,VF[,1])
    VFz2 = interp1d(z,x,VF[,2])
    
    return(cbind(VFz1,VFz2))
}

interp1d <- function(z,x,f){
    # z = where to interpolate <- (nInterp x 2)
    # x = evaluation points where f is computed <- (nEval x 2)
    # f = function to interpolate <- (nEval)
    
    xS = sort(unique(x[,1]))
    yS = sort(unique(x[,2]))
    nEval = length(xS)
    interpZ = interp2(xS,yS,matrix(f,length(xS),length(yS),byrow=T),z[,1],z[,2])
    return(interpZ)
}

# fx = x[,1]^2 + 0.1*x[,2]^3
# image.plot(x = unique(x[,1]), y = unique(x[,2]), z = matrix((fx), nrow=sqrt(nEval), ncol=sqrt(nEval)),xlab="x",ylab="y",main="error norm rel")
# image.plot(x = unique(z[,1]), y = unique(z[,2]), z = matrix((interp1d(z,x,fx)), nrow=length(unique(z[,1])), ncol=length(unique(z[,2]))),xlab="x",ylab="y",main="error norm rel")


