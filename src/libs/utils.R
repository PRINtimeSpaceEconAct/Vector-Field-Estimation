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