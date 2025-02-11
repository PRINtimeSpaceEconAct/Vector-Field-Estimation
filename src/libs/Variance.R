
significance <- function(est,X0,X1){
    
    Y = X1-X0
    Yhat = interp2d(X0,est$x,est$estimator)
    eps = Y - Yhat    
    
}
    




