
significanceVF <- function(est,X0,X1,alpha = 0.05){
    

    h = est$h
    
    n = nrow(X0)
    Y = X1-X0
    Yhat = interp2d(X0,est$x,est$estimator)
    eps = Y - Yhat    
    
    # homoskedastic errors
    Sigma = cov(eps, use = "complete.obs")

    if (est$kernel.type == "epa"){ k = 3/5 
    } else if (est$kernel == "gauss") { k = 1/(2*sqrt(pi)) }
    
    
    Var = array(data = NA, dim = c(2,2,nrow(x)))
    ChiSquare_stat = rep(NA, nrow(x))     
    p_values = rep(NA, nrow(x))
    signif = rep(NA, nrow(x))
    for (i in 1:nrow(x)){
        Var[,,i] = (k^2 * Sigma / (est$density[i] * n*h^2)) 
        ChiSquare_stat[i] = t(est$estimator[i,]) %*% solve(Var[,,i]) %*% est$estimator[i,]
        p_values[i] = 1-pchisq(ChiSquare_stat[i], 2)
        p_values[i] = ifelse( is.na(p_values[i]), 1 , p_values[i])
        signif[i] = ifelse(p_values[i] < alpha, 1, 0)
    }
    
    
    return(listN(Var,ChiSquare_stat,p_values,signif))
}






