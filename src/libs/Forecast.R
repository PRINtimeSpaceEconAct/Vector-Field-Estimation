forecastDiscrete <- function(z,estimatedVF,speedFactor=1,nPeriods=1){
    # z = points to forecast <- (nPts x 2)  
    # v = estimatedVF <- output of the function that estimate the VF
    # speedfactor = scaling for the length of arrows (i.e slows time)
    # number of periods to forecast
    # zForecasted = forecasted points <- (nPts x 2 x nPeriods)
    # if nPeriods == 1, then zForecasted <- (nPts x 2)
    
    # set to zero the NA to perform the interpolation    
    estimatedVF$estimator[is.na(estimatedVF$estimator)] = 0
    
    if (nPeriods > 1) {
        zForecasted <- array(0,dim=c(nrow(z),2,nPeriods))
        
        if (DEBUG) print(paste("Forecasting period ",1, "over",nPeriods))
        zForecasted[,,1] <- z +  speedFactor*interp2d(z,estimatedVF$x,estimatedVF$estimator) 
        for (t in 2:nPeriods){
            if (DEBUG) print(paste("Forecasting period ",t, "over",nPeriods))
            zForecasted[,,t] <- zForecasted[,,t-1] +  speedFactor*interp2d(zForecasted[,,t-1],estimatedVF$x,estimatedVF$estimator)
        }
    }
    else {
        zForecasted = matrix(0,nrow=nrow(z),ncol=ncol(z))
        zForecasted <- z + speedFactor*interp2d(z,estimatedVF$x,estimatedVF$estimator)
    }
    
    return (zForecasted)
}

# nPeriods = 10
# forecast = forecastDiscrete(X0, est_field_adaptive, speedFactor = 1, nPeriods = nPeriods)
# for (i in 1:nPeriods) {
#     plot(est_field_adaptive$x, type = "n", xlab = "X", ylab = "Y", main = "Estimated Vector Field")
#         arrows(est_field_adaptive$x[,1], est_field_adaptive$x[,2],
#         est_field_adaptive$x[,1] + est_field_adaptive$estimator[,1],
#         est_field_adaptive$x[,2] + est_field_adaptive$estimator[,2],
#         length = 0.05, angle = 15, col = "blue")
#     points(forecast[,,i],cex=0.5)
#     Sys.sleep(0.1)
#     }
