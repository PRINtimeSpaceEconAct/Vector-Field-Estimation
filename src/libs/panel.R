#' Computes weights for the within transformation.
#'
#' @param X A 3D array representing the panel data of size nObs x 2 x nT.
#' @param nEval_chunk The number of evaluation points.
#' @return A list containing the weights `omega_a`, `omega_b`, and `omega_c`.
#' @export
compute_weights <- function(X, nEval_chunk) {
    # This is a placeholder for the actual weight computation.
    # The actual implementation will be provided later.
    stop("compute_weights is not implemented yet.")
}

#' Applies the FE/TE transformations to a given dataset
apply_transformation <- function(data, nT_data, omega_a, omega_b, omega_c, FE, TE, nEval_chunk) {
    data_dims <- dim(data)
    nObs_data <- data_dims[1]
    
    data_ext <- array(data, dim = c(data_dims, nEval_chunk))

    if (!FE && !TE) {
        return(data_ext)
    }

    term_a <- term_b <- term_c <- array(0, dim = dim(data_ext))

    if (FE) { # Corresponds to "0,."
        omega_a_perm <- aperm(omega_a, c(1, 3, 2)) # i, e, s
        data_slice_perm <- aperm(data, c(1, 3, 2))   # i, s, d
        term_a_sum <- array(0, dim = c(nObs_data, 2, nEval_chunk))
        for (i in 1:nObs_data) {
            for (e in 1:nEval_chunk) {
                term_a_sum[i,,e] <- omega_a_perm[i, e, ] %*% data_slice_perm[i, , ]
            }
        }
        term_a_sum_ext <- array(term_a_sum, dim = c(nObs_data, 2, nEval_chunk, nT_data))
        term_a <- aperm(term_a_sum_ext, c(1, 2, 4, 3))
    }

    if (TE) { # Corresponds to ".,0"
        omega_b_perm <- aperm(omega_b, c(1, 3, 2)) # j, e, s
        data_slice_perm <- aperm(data, c(1, 3, 2)) # j, s, d
        term_b_sum <- array(0, dim = c(nEval_chunk, nT_data, 2))
        for(e in 1:nEval_chunk){
            for(s in 1:nT_data){
                term_b_sum[e,s,] <- colSums(omega_b_perm[,e,s] * data_slice_perm[,s,])
            }
        }
        term_b_slice <- aperm(term_b_sum, c(3,2,1))
        term_b_slice_ext <- array(term_b_slice, dim = c(2, nT_data, nEval_chunk, nObs_data))
        term_b <- aperm(term_b_slice_ext, c(4, 1, 2, 3))
    }

    if (FE && TE) { # Corresponds to "0,0"
        term_c_sum <- array(0, dim = c(2, nEval_chunk))
        for(e in 1:nEval_chunk){
            for(s in 1:nT_data){
                term_c_sum[,e] <- term_c_sum[,e] + colSums(omega_c[,s,e] * data[,,s])
            }
        }
        term_c_sum_ext <- array(term_c_sum, dim = c(2, nEval_chunk, nObs_data, nT_data))
        term_c <- aperm(term_c_sum_ext, c(3, 1, 4, 2))
    }

    if (FE && !TE) { # "0,."
        transformed_data <- data_ext - term_a
    } else if (!FE && TE) { # ".,0"
        transformed_data <- data_ext - term_b
    } else if (FE && TE) { # "0,0"
        transformed_data <- data_ext - term_a - term_b + term_c
    } else {
        transformed_data <- data_ext
    }
    return(transformed_data)
}

#' Applies within-group transformations to panel data.
#'
#' This function computes transformations on a panel time series X.
#' It can perform within-individual, within-time, or two-ways transformations.
#'
#' @param X A 3D array of size nObs x 2 x nT representing the panel data.
#' @param FE A logical value indicating whether to include fixed effects (within-individual transformation).
#' @param TE A logical value indicating whether to include time effects (within-time transformation).
#' @param uniform_weights A logical value. If TRUE, uniform weights are used.
#'                        If FALSE, weights are computed using `weights_function`.
#' @param weights_function A function to compute the weights. It is used when
#'                         `uniform_weights` is FALSE. It should take `X` and `nEval_chunk`
#'                         as input and return a list of weights.
#' @param nEval_chunk The number of evaluation points for the weights. Required if `uniform_weights` is FALSE.
#'
#' @return A 4D array of size nObs x 2 x nT x nEval_chunk with the transformed data.
#'
#' @details The transformations are defined as follows:
#' - `FE=TRUE`, `TE=FALSE`:   \(X_{i,t}^{0,\cdot}(x) = X_{i,t} - \sum_{s=1}^{T}\omega_{i,s}^{a}(x)X_{i,s}\)
#' - `FE=FALSE`, `TE=TRUE`:   \(X_{i,t}^{\cdot,0}(x) = X_{i,t} - \sum_{j=1}^{N}\omega_{j,t}^{b}(x)X_{j,t}\)
#' - `FE=TRUE`, `TE=TRUE`:   \(X_{i,t}^{0,0}(x) = X_{i,t} - \sum_{s=1}^{T}\omega_{i,s}^{a}(x)X_{i,s} - \sum_{j=1}^{N}\omega_{j,t}^{b}(x)X_{j,t} + \sum_{j=1}^{N}\sum_{s=1}^{T}\omega_{j,s}^{c}(x)X_{j,s}\)
#'
#' @export
within_transform <- function(X,
                             FE = FALSE,
                             TE = FALSE,
                             uniform_weights = TRUE,
                             nEval_chunk = 1000) {

    dims <- dim(X)
    if (length(dims) != 3 || dims[2] != 2) {
        stop("Input X must be a 3D array with second dimension of size 2.")
    }
    nObs <- dims[1]
    nT <- dims[3]

    # Compute Y as the time difference of X
    Y <- X[,, 2:nT, drop = FALSE] - X[,, 1:(nT - 1), drop = FALSE]

    # Compute weights for X
    if (uniform_weights) {
        omega_a_X <- array(1/nT, dim = c(nObs, nT, nEval_chunk))
        omega_b_X <- array(1/nObs, dim = c(nObs, nT, nEval_chunk))
        omega_c_X <- array(1/(nObs * nT), dim = c(nObs, nT, nEval_chunk))
    } else {
        weights_X <- compute_weights(X, nEval_chunk)
        omega_a_X <- weights_X$omega_a
        omega_b_X <- weights_X$omega_b
        omega_c_X <- weights_X$omega_c
    }

    # Transform X and compute X0_star
    X_transformed <- apply_transformation(X, nT, omega_a_X, omega_b_X, omega_c_X, FE, TE, nEval_chunk)
    X0_star <- X_transformed[,, 1:(nT - 1), ]

    # Compute weights for Y
    nT_Y <- nT - 1
    if (uniform_weights) {
        omega_a_Y <- array(1/nT_Y, dim = c(nObs, nT_Y, nEval_chunk))
        omega_b_Y <- array(1/nObs, dim = c(nObs, nT_Y, nEval_chunk))
        omega_c_Y <- array(1/(nObs * nT_Y), dim = c(nObs, nT_Y, nEval_chunk))
    } else {
        # Assuming compute_weights can handle data with different nT
        weights_Y <- compute_weights(Y, nEval_chunk)
        omega_a_Y <- weights_Y$omega_a
        omega_b_Y <- weights_Y$omega_b
        omega_c_Y <- weights_Y$omega_c
    }
    
    # Transform Y
    Y_star <- apply_transformation(Y, nT_Y, omega_a_Y, omega_b_Y, omega_c_Y, FE, TE, nEval_chunk)

    # Unroll X0_raw
    X0_raw <- X[,,1:(nT-1), drop=FALSE]
    X0_raw_unrolled <- aperm(X0_raw, c(3, 1, 2))
    dim(X0_raw_unrolled) <- c((nT-1) * nObs, 2)

    # Unroll X0_star
    X0_star_unrolled <- aperm(X0_star, c(3, 1, 2, 4))
    dim(X0_star_unrolled) <- c((nT-1) * nObs, 2, nEval_chunk)

    # Unroll Y_star
    Y1_unrolled <- aperm(Y_star[,1,,], c(2, 1, 3))
    dim(Y1_unrolled) <- c((nT-1) * nObs, nEval_chunk)
    Y2_unrolled <- aperm(Y_star[,2,,], c(2, 1, 3))
    dim(Y2_unrolled) <- c((nT-1) * nObs, nEval_chunk)

    return(listN(X0_star_unrolled, Y1_unrolled, Y2_unrolled, X0_raw_unrolled))
}

compute_derivative_term <- function(X, X_star, x=NULL, nEval=2500, kernel.type="gauss", D=NULL, 
                             method.h=NULL, h=NULL, lambda=NULL, 
                             sparse=FALSE, gc=FALSE, chunk_size=nrow(x), Y=NULL) {
    
    nObs = nrow(X)
    covX = cov(X)
    invS = solve(covX)
    sqrtinvS = expm::sqrtm(solve(covX))
    detS = det(covX)
    
    Z = X %*% sqrtinvS
    Z_star <- simplify2array(purrr::map(1:dim(X_star)[3], ~ X_star[,,.x] %*% sqrtinvS))
    z = x %*% sqrtinvS
    
    if (is.null(x)) { x = defineEvalPoints(X,nEval) }
    nEval =  nrow(x)
    if (DEBUG) {
        print(paste("nEval: ",nEval))
        print(paste("nObs: ",nObs)) }
    if (is.null(chunk_size)) { chunk_size = nrow(x) }
    if (is.null(lambda)) { lambda = rep(1,nObs) }
    
    if (is.null(h) | is.null(method.h)) { 
        list.h = define_h_method.h(Z,h,method.h, kernel.type)
        h = list.h$h
        method.h = list.h$method.h }
    kernelFunction = defineKernel(kernel.type)
    
    
    estimator <- matrix(NA, nrow = nEval, ncol = 2)


    chunks = split(seq_len(nEval), ceiling(seq_len(nEval)/chunk_size))
    
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
        
        K_scaledY = K_scaled * Y[, chunk, drop = FALSE]
        d1 = Z_star[,1,chunk]
        d2 = Z_star[,2,chunk]
    
        # Calculate local moments

        S12 = colSums(K_scaled * d1 * d2)
        S11 = colSums(K_scaled * d1^2)
        S22 = colSums(K_scaled * d2^2)

        T1 = colSums(K_scaledY * d1)
        T2 = colSums(K_scaledY * d2)
    
        # MDenominator[,,k] is the S_k matrix for the k-th observation IN THE CHUNK
        MDenominator = array(NA,dim=c(2,2,chunk_size))
        MDenominator[1,1,] = S11
        MDenominator[1,2,] = S12
        MDenominator[2,1,] = S12
        MDenominator[2,2,] = S22
        
        # Keep calculation for solutions (estimator)
        bConstant = rbind(T1,T2)
        solve_system <- function(k) { # k is index within chunk (1 to nEval_chunk)
            if ( det(MDenominator[,,k]) == 0 ) {
                return(rep(NaN,2))
            }
            return(as.numeric(ginv(MDenominator[,,k]) %*% bConstant[,k]))
        }
        solutions <- matrix(unlist(purrr::map(1:chunk_size, ~ solve_system(.x))),
                            nrow=2,byrow = FALSE)

        estimator[chunk, ] <- t(solutions)

        if (gc == TRUE){ gc() }
    }
    return(listN(x, estimator, h, method.h, kernel.type))
}