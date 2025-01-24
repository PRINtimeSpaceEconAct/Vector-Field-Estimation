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
