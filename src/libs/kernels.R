

epaKernel <- function(z){
    # to be computed on z S^-1 z
    2/pi*(1-z)*(z <= 1)
}

gaussKernel <- function(z){
    # to be computed on z S^-1 z
    1/(2*pi)*exp(-z/2)
}



