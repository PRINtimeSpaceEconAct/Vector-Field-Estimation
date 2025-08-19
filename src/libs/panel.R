#' Computes weights for the within transformation.
#'
#' @param X A 3D array representing the panel data of size nObs x 2 x nT.
#' @param nEval_chunk The number of evaluation points.
#' @return A list containing the weights `omega_a`, `omega_b`, and `omega_c`.
#' @export
compute_weights <- function(X, x=NULL, nEval=2500, kernel.type="gauss", D=NULL, 
                             method.h=NULL, h=NULL, lambda=NULL,
                             sparse=FALSE, gc=FALSE, chunk_size=NULL){    
    dims = dim(X)
    nObs = dims[1]
    nT = dims[3]
    X_unrolled = aperm(X, c(1, 3, 2))
    dim(X_unrolled) = c(nObs * nT, 2)
    nTotalObs = nObs * nT

    covX = cov(X_unrolled)
    invS = solve(covX)
    sqrtinvS = expm::sqrtm(solve(covX))
    detS = det(covX)
    
    Z = X_unrolled %*% sqrtinvS
    
    if (is.null(x)) { x = defineEvalPoints(X_unrolled,nEval) }
    nEval =  nrow(x)
    z = x %*% sqrtinvS

    if (DEBUG) {
        print(paste("nEval: ",nEval))
        print(paste("nObs: ",nObs))
        print(paste("nT: ",nT))
    }
    if (is.null(chunk_size)) { chunk_size = nEval }
    if (is.null(lambda)) { lambda = rep(1, nTotalObs) }
    
    if (is.null(h) | is.null(method.h)) { 
        list.h = define_h_method.h(Z,h,method.h, kernel.type)
        h = list.h$h
        method.h = list.h$method.h 
    }
    kernelFunction = defineKernel(kernel.type)
    
    omega_a <- array(NA, dim = c(nObs, nT, nEval))
    omega_b <- array(NA, dim = c(nObs, nT, nEval))
    omega_c <- array(NA, dim = c(nObs, nT, nEval))

    chunks = split(seq_len(nEval), ceiling(seq_len(nEval)/chunk_size))
    
    if (DEBUG) {
        print(paste("Computing ",length(chunks)," chunks"))
        print(paste("Chunk size: ",chunk_size)) 
    }
    
    # start estimate
    for(i in 1:length(chunks)) {
        if (DEBUG) print(paste("Computing chunk ",i,"/",length(chunks),sep=""))
        
        chunk = chunks[[i]]
        current_chunk_size = length(chunk)
        z_chunk = z[chunk, ,drop=FALSE]
        
        if (is.null(D)){ 
            D_chunk = computeDcomponents(Z, z_chunk, sparse=sparse) 
        } else { 
            D_chunk = list(z1=D$z1[,chunk], z2=D$z2[,chunk]) 
        }
        
        # Kernel computation
        K = kernelFunction(sweep(D_chunk$z1, 1, h * lambda, "/"),sweep(D_chunk$z2, 1, h * lambda, "/"))
        
        dim(K) <- c(nObs, nT, current_chunk_size)

        # Compute denominators
        denom_a <- apply(K, c(1, 3), sum, na.rm = TRUE)
        denom_b <- apply(K, c(2, 3), sum, na.rm = TRUE)
        denom_c <- apply(K, 3, sum, na.rm = TRUE)
        
        # To prevent division by zero
        denom_a[denom_a == 0] <- 1
        denom_b[denom_b == 0] <- 1
        denom_c[denom_c == 0] <- 1

        # Compute weights for the chunk
        omega_a[,,chunk] <- sweep(K, c(1, 3), denom_a, "/")
        omega_b[,,chunk] <- sweep(K, c(2, 3), denom_b, "/")
        omega_c[,,chunk] <- sweep(K, 3, denom_c, "/")

        if (gc == TRUE){ gc() }
    }

    return(listN(omega_a, omega_b, omega_c))
}

#' Applies fixed effect (FE) and/or time effect (TE) transformations to panel data.
#'
#' This function takes a 3D array of panel data and applies within-individual (FE),
#' within-time (TE), or two-way transformations based on the provided weights.
#' The transformations subtract weighted averages from the original data.
#'
#' @param data A 3D array of panel data with dimensions (nObs, 2, nT).
#' @param nT_data The number of time periods in `data`.
#' @param omega_a A 3D array of weights for the individual-specific transformation.
#' @param omega_b A 3D array of weights for the time-specific transformation.
#' @param omega_c A 3D array of weights for the overall transformation.
#' @param FE A logical value. If TRUE, the fixed-effect (within-individual) transformation is applied.
#' @param TE A logical value. If TRUE, the time-effect (within-time) transformation is applied.
#' @param nEval_chunk The number of evaluation points, determining the fourth dimension of the output array.
#' @return A 4D array containing the transformed data, with dimensions (nObs, 2, nT, nEval_chunk).
#' @export
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
                             nEval_chunk = 1000,
                             x = NULL,
                             kernel.type = "gauss",
                             D = NULL,
                             method.h = NULL,
                             h = NULL,
                             lambda = NULL,
                             sparse = FALSE,
                             gc = FALSE,
                             chunk_size = NULL) {

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
        weights_X <- compute_weights(X = X,
                                     x = x,
                                     nEval = nEval_chunk,
                                     kernel.type = kernel.type,
                                     D = D,
                                     method.h = method.h,
                                     h = h,
                                     lambda = lambda,
                                     sparse = sparse,
                                     gc = gc,
                                     chunk_size = chunk_size)
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
        # For Y, we reset data-dependent parameters
        weights_Y <- compute_weights(X = Y,
                                     x = NULL,
                                     nEval = nEval_chunk,
                                     kernel.type = kernel.type,
                                     D = NULL,
                                     method.h = method.h,
                                     h = NULL,
                                     lambda = NULL,
                                     sparse = sparse,
                                     gc = gc,
                                     chunk_size = chunk_size)
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
    Y1_star_unrolled <- aperm(Y_star[,1,,], c(2, 1, 3))
    dim(Y1_star_unrolled) <- c((nT-1) * nObs, nEval_chunk)
    Y2_star_unrolled <- aperm(Y_star[,2,,], c(2, 1, 3))
    dim(Y2_star_unrolled) <- c((nT-1) * nObs, nEval_chunk)
    
    # Unroll also Y
    Y_unrolled <- aperm(Y, c(3, 1, 2))
    dim(Y_unrolled) <- c((nT-1) * nObs, 2)
    Y1_unrolled = Y_unrolled[,1]
    Y2_unrolled = Y_unrolled[,2]
    
    return(listN(X0_star_unrolled, Y1_star_unrolled, Y2_star_unrolled, X0_raw_unrolled,Y1_unrolled,Y2_unrolled))
}

#' Computes the derivative term using a local linear estimator.
#'
#' This function estimates the derivative of a dependent variable Y with respect to
#' independent variables X at specified evaluation points x. It uses a local
#' linear kernel regression approach, handling panel data transformations.
#'
#' @param X A numeric matrix of observations (nObs x 2).
#' @param X_star A 3D array representing the transformed independent variables (nObs x 2 x nEval).
#' @param x A matrix of evaluation points (nEval x 2).
#' @param nEval The number of evaluation points, used if x is not provided.
#' @param kernel.type The type of kernel to use (e.g., "gauss").
#' @param D A list containing precomputed distance components.
#' @param method.h The method for bandwidth selection.
#' @param h The bandwidth value.
#' @param lambda A numeric vector of observation-specific weights.
#' @param sparse A logical indicating whether to use sparse matrices.
#' @param gc A logical indicating whether to perform garbage collection.
#' @param chunk_size The number of evaluation points to process in each chunk.
#' @param Y A matrix representing the dependent variable (nObs x nEval).
#' @return A list containing the evaluation points `x`, the estimated derivatives `estimator`,
#'         the bandwidth `h`, the bandwidth selection method `method.h`, and the kernel type `kernel.type`.
compute_derivative_term <- function(X, X_star, x=NULL, nEval=2500, kernel.type="gauss", D=NULL, 
                             method.h=NULL, h=NULL, lambda=NULL, 
                             sparse=FALSE, gc=FALSE, chunk_size=nrow(x), Y=NULL) {
    
    nObs = nrow(X)
    covX = cov(X)
    invS = solve(covX)
    sqrtinvS = expm::sqrtm(solve(covX))
    detS = det(covX)
    nEval = nrow(x)
    
    Z = X %*% sqrtinvS
    Z_star <- simplify2array(purrr::map(1:nEval, ~ X_star[,,.x] %*% sqrtinvS))
    z = x %*% sqrtinvS
    
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
        current_chunk_size = length(chunk)
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
        MDenominator = array(NA,dim=c(2,2,current_chunk_size))
        MDenominator[1,1,] = S11
        MDenominator[1,2,] = S12
        MDenominator[2,1,] = S12
        MDenominator[2,2,] = S22
        
        # Keep calculation for solutions (estimator)
        bConstant = rbind(T1,T2)
        solve_system <- function(k) { # k is index within chunk (1 to current_chunk_size)
            if ( det(MDenominator[,,k]) == 0 ) {
                return(rep(NaN,2))
            }
            
            # GINV
            # sol = as.numeric(ginv(MDenominator[,,k]) %*% bConstant[,k])
            
            # RIDGE
            lambda = 1e-5 # small ridge parameter to avoid singularity
            A = MDenominator[,,k]
            sol = as.numeric(solve(A + lambda * diag(2), bConstant[,k]))
            
            return(sol)
        }
        solutions <- matrix(unlist(purrr::map(1:current_chunk_size, ~ solve_system(.x))),
                            nrow=2,byrow = FALSE)

        estimator[chunk, ] <- t(solutions)

        if (gc == TRUE){ gc() }
    }
    return(listN(x, estimator, h, method.h, kernel.type))
}

#' Reconstructs the scalar potential m(x) from gradient estimates
#'
#' Given estimates of the gradient field \(\beta(x) = \nabla m(x)\) and an
#' anchoring condition \(m(x_0) = m_0\), reconstructs the potential at
#' evaluation points using
#'\preformatted{
#' m(x) = m_0 + (1/NT) * sum_{i,t} [ beta(x_0)^T (X_{i,t} - x_0)
#'                                   - beta(x)^T   (X_{i,t} - x)   ]
#'}
#' where `NT` is the number of unrolled observed points used in estimation.
#'
#' Input relationships and guidance:
#' - X_obs_unrolled (NT x 2): unrolled design points used elsewhere (e.g. `X0_raw`).
#' - x_eval (nEval x 2): grid of evaluation points where m(x) is returned.
#' - beta (nEval x 2): gradient estimates at each row of `x_eval`.
#' - x0 (length-2), beta_0 (length-2): anchor point and its gradient; m_0 = m(x0)
#'   is computed by `compute_m0` using a mean-matching moment.
#'
#' @param X_obs_unrolled Numeric matrix (NT x 2) of observed regressor points.
#' @param x_eval Numeric matrix (nEval x 2) of evaluation points.
#' @param beta Numeric matrix (nEval x 2): gradient estimates at `x_eval`.
#' @param m_0 Numeric scalar: anchoring value `m(x0)`.
#' @param x0 Numeric vector length 2: anchor point `x0`.
#' @param beta_0 Numeric vector length 2: gradient at `x0`.
#' @return Numeric vector length nEval with reconstructed m(x).
#' @export
compute_m <- function(X_obs_unrolled, x_eval, beta, m_0, x0, beta_0) {

    NT = dim(X_obs_unrolled)[1]

    # Compute the first term of the sum: sum_{i,t} beta(x0)^T * (X_{i,t} - x0)
    # This does not depend on the evaluation points x
    diff_x0 <- sweep(X_obs_unrolled, 2, x0, "-")
    term1 <- diff_x0 %*% beta_0
    S1 <- sum(term1)

    # Compute the second term of the sum: sum_{i,t} beta(x)^T * (X_{i,t} - x)
    # This depends on the evaluation points x
    
    # Use computeDcomponents to get (x - X)
    # The result z1, z2 are matrices of size (N*T) x nEval
    # z1[it, e] = x[e,1] - X_unrolled[it, 1]
    D_components <- computeDcomponents(X_obs_unrolled, x_eval, sparse = FALSE)
    
    # We need (X - x), so we negate the results
    Diff1 <- -D_components$z1
    Diff2 <- -D_components$z2
    
    # S2 will be a vector of size nEval
    # For each evaluation point e, we compute S2[e] = sum_{i,t} (beta_e^T * (X_{i,t} - x_e))
    # beta_e^T * (X_{it} - x_e) = beta[e,1]*(X_{it,1}-x_{e,1}) + beta[e,2]*(X_{it,2}-x_{e,2})
    # This is beta[e,1]*Diff1[it,e] + beta[e,2]*Diff2[it,e]
    
    # Vectorized computation of the second term
    term2_1 <- sweep(Diff1, 2, beta[, 1], "*")
    term2_2 <- sweep(Diff2, 2, beta[, 2], "*")
    term2 <- term2_1 + term2_2
    
    S2 <- colSums(term2)
    
    # Final computation of m(x)
    # m(x) = m_0 + (1/N/T) * (S1 - S2)
    m_x <- m_0 + (1 / NT) * (S1 - S2)
    
    return(m_x)
}

#' Computes the anchoring constant m_0 via mean-matching
#'
#' Fixes the additive indeterminacy of m(·) by imposing the sample mean
#' condition
#'\preformatted{
#'   (1/NT) * sum_{i,t} Y_{i,t}  =  (1/NT) * sum_{i,t} m(X_{i,t}).
#'}
#' Using the same representation as in `compute_m` and solving for `m_0` gives
#'\preformatted{
#' m_0 = (1/NT) * [ sum_{i,t} Y_{i,t}
#'                  - sum_{i,t} beta(x_0)^T (X_{i,t} - x_0)
#'                  + (1/NT) * sum_{j,s} beta(X_{j,s})^T sum_{i,t} (X_{i,t} - X_{j,s}) ].
#'}
#' The code implements this as three terms:
#' - Sy = sum of Y over all unrolled observations
#' - S1 = sum of the reference linear term beta_0^T (X_{i,t} - x0)
#' - S2 = (1/NT) * sum_{j,s} beta(X_{j,s})^T sum_{i,t} (X_{i,t} - X_{j,s})
#'   (the outer 1/NT in the return line yields the overall 1/NT^2 on the
#'   pairwise double sum).
#'
#' Inputs should be aligned with the rest of the pipeline:
#' - X_obs_unrolled (NT x 2): same unrolled design matrix used for derivatives (e.g. `X0_raw`).
#' - Y_obs_unrolled (NT): component of observed changes aligned with rows of X (e.g. `Y1` or `Y2`).
#' - beta (NT x 2): gradient evaluated at each observed row (e.g. `derivative_obs_*$estimator`).
#' - x0, beta_0: chosen anchor and its gradient; using x0 at the sample mean
#'   typically makes S1 ≈ 0 and improves stability.
#'
#' @param X_obs_unrolled Numeric matrix (NT x 2) of observed regressor points.
#' @param Y_obs_unrolled Numeric vector length NT with observed dependent values.
#' @param beta Numeric matrix (NT x 2): gradients at the observed points.
#' @param x0 Numeric vector length 2: anchor point.
#' @param beta_0 Numeric vector length 2: gradient at `x0`.
#' @return Numeric scalar with the integration constant m_0.
compute_m0 <- function(X_obs_unrolled, Y_obs_unrolled, beta, x0, beta_0) {
    NT = dim(X_obs_unrolled)[1]
    
    # S1 = sum_{i,t} beta(x0)^T (X_{i,t} - x0)
    diff_x0 <- sweep(X_obs_unrolled, 2, x0, "-")
    term1 <- diff_x0 %*% beta_0
    S1 <- sum(term1)
    
    # S2 = (1/NT) * sum_{j,s} beta(X_{j,s})^T sum_{i,t} (X_{i,t} - X_{j,s})
    D_components <- computeDcomponents(X_obs_unrolled, X_obs_unrolled, sparse = FALSE)
    
    Diff1 <- -D_components$z1
    Diff2 <- -D_components$z2
    
    term2_1 <- sweep(Diff1, 2, beta[, 1], "*")
    term2_2 <- sweep(Diff2, 2, beta[, 2], "*")
    term2 <- term2_1 + term2_2
    
    S2 <- 1/NT * sum(term2)
    
    # Sy = sum_{i,t} Y_{i,t}
    Sy = sum(Y_obs_unrolled) 
    
    m_0 <- (1 / NT) * (Sy - S1 + S2)
    
    return(m_0)
}

#' Estimates a vector field from panel data with fixed and time effects.
#'
#' This function implements a procedure to estimate a vector field from panel data,
#' allowing for the removal of fixed and time effects. It wraps several steps:
#' data transformation, derivative estimation, and vector field reconstruction.
#'
#' @param X A 3D array of size nObs x 2 x nT representing the panel data.
#' @param x A matrix of evaluation points of size nEval x 2.
#' @param nEval The number of evaluation points for the weights.
#' @param FE A logical value indicating whether to include fixed effects.
#' @param TE A logical value indicating whether to include time effects.
#' @param uniform_weights A logical value for using uniform weights in transformation.
#' @param kernel.type The kernel function to be used.
#' @param method.h The method for bandwidth selection.
#' @param chunk_size The size of chunks for processing.
#' @param sparse Logical, for sparse matrix calculations.
#' @param gc Logical, for garbage collection.
#' @return A list containing the estimated vector field `VF_hat`, the evaluation points `x`, and other intermediate results.
#' @export
estimate_panel_vf <- function(X,
                              x,
                              nEval,
                              FE = TRUE,
                              TE = TRUE,
                              uniform_weights = TRUE,
                              kernel.type = "gauss",
                              method.h = "silverman",
                              h = NULL,
                              chunk_size = 512,
                              sparse = FALSE,
                              gc = FALSE) {

    # 1. Within transformation
    Filtered <- within_transform(X,
                                 FE = FE, TE = TE,
                                 uniform_weights = uniform_weights, nEval_chunk = nEval,
                                 x = x, kernel.type = kernel.type,
                                 method.h = method.h, h = h, chunk_size = chunk_size)

    X0_raw <- Filtered$X0_raw_unrolled
    X0_star <- Filtered$X0_star_unrolled
    Y1_star <- Filtered$Y1_star_unrolled
    Y2_star <- Filtered$Y2_star_unrolled
    Y1 <- Filtered$Y1_unrolled
    Y2 <- Filtered$Y2_unrolled
    
    nObs = dim(X)[1]
    nT = dim(X)[3]

    # 2. Estimate derivatives at evaluation points
    derivative_estimator_1 <- compute_derivative_term(X0_raw, X_star = X0_star, x=x,
                                                      kernel.type=kernel.type, D=NULL,
                                                      method.h=method.h, h=h, lambda=NULL,
                                                      sparse=sparse, gc=gc, chunk_size=chunk_size, Y=Y1_star)

    derivative_estimator_2 <- compute_derivative_term(X0_raw, X_star = X0_star, x=x,
                                                      kernel.type=kernel.type, D=NULL,
                                                      method.h=method.h, h=h, lambda=NULL,
                                                      sparse=sparse, gc=gc, chunk_size=chunk_size, Y=Y2_star)

    # 3. Estimate derivatives at observed points
  
    X_star_obs <- X0_star[,, rep(1, nrow(X0_raw))]
    Y1_star_obs <- Y1_star[, rep(1, nrow(Y1_star))]
    Y2_star_obs <- Y2_star[, rep(1, nrow(Y2_star))]
    
    derivative_obs_1 <- compute_derivative_term(X0_raw, X_star=X_star_obs, x=X0_raw,
                                                kernel.type=kernel.type, D=NULL,
                                                method.h=method.h, h=h, lambda=NULL,
                                                sparse=sparse, gc=gc, chunk_size=chunk_size, Y=Y1_star_obs)

    derivative_obs_2 <- compute_derivative_term(X0_raw, X_star=X_star_obs, x=X0_raw,
                                                kernel.type=kernel.type, D=NULL,
                                                method.h=method.h, h=h, lambda=NULL,
                                                sparse=sparse, gc=gc, chunk_size=chunk_size, Y=Y2_star_obs)

    # 4. Compute m0
    meanPoint <- colSums(X0_raw) / nrow(X0_raw)
    iBest <- which.min(rowSums(sweep(x, 2, meanPoint)^2))

    m10 <- compute_m0(X_obs_unrolled=X0_raw, Y_obs_unrolled=Y1, beta=derivative_obs_1$estimator, x0=x[iBest,], beta_0=derivative_estimator_1$estimator[iBest,])
    m20 <- compute_m0(X_obs_unrolled=X0_raw, Y_obs_unrolled=Y2, beta=derivative_obs_2$estimator, x0=x[iBest,], beta_0=derivative_estimator_2$estimator[iBest,])
    
    # 5. Compute vector field
    VF_hat1 <- compute_m(X_obs_unrolled = X0_raw, x_eval = x, beta=derivative_estimator_1$estimator, m_0=m10, x0=x[iBest,], beta_0=derivative_estimator_1$estimator[iBest,])
    VF_hat2 <- compute_m(X_obs_unrolled = X0_raw, x_eval = x, beta=derivative_estimator_2$estimator, m_0=m20, x0=x[iBest,], beta_0=derivative_estimator_2$estimator[iBest,])

    estimator <- cbind(VF_hat1, VF_hat2)
    h <- derivative_estimator_1$h

    return(listN(estimator, x, X0_raw, X, nEval, FE, TE, uniform_weights, kernel.type, method.h, h, chunk_size, sparse, gc,
                 derivative_estimator_1, derivative_estimator_2, derivative_obs_1, derivative_obs_2, 
                 iBest, m10, m20, VF_hat1, VF_hat2, Y1, Y2))
}

#' Reconstructs fixed and time effects from an estimated panel vector field model.
#'
#' After estimating the vector field, this function computes the individual fixed effects (alpha_i)
#' and time effects (gamma_t) based on the residuals between the observed changes and the
#' estimated vector field. The function assumes that `X_obs` has `(nT-1)*nObs` rows.
#'
#' @param panel_estimator A list object returned by `estimate_panel_vf`, containing all estimation results.
#' @param X_obs A matrix of observation points where the effects are to be computed, typically `panel_estimator$X0_raw`.
#' @param FE A logical value. If TRUE, computes fixed effects.
#' @param TE A logical value. If TRUE, computes time effects.
#' @return A list containing the estimated fixed effects `alpha_i` (a matrix of size nObs x 2)
#'         and time effects `gamma_t` (a matrix of size (nT-1) x 2).
get_effects <- function(panel_estimator, X_obs, FE=TRUE, TE=TRUE) {
    x <- panel_estimator$x
    X <- panel_estimator$X
    nObs <- dim(X)[1]
    nT <- dim(X)[3]

    X_obs_estimator <- panel_estimator$X0_raw
    derivative_obs_1 <- panel_estimator$derivative_obs_1
    derivative_obs_2 <- panel_estimator$derivative_obs_2
    derivative_estimator_1 <- panel_estimator$derivative_estimator_1
    derivative_estimator_2 <- panel_estimator$derivative_estimator_2
    iBest <- panel_estimator$iBest
    m10 <- panel_estimator$m10
    m20 <- panel_estimator$m20
    Y1 <- panel_estimator$Y1
    Y2 <- panel_estimator$Y2

    # reconstruct FE
    VF_hat1 = compute_m(X_obs_unrolled = X_obs_estimator, x_eval = X_obs, beta=derivative_obs_1$estimator, m_0=m10, x0=x[iBest,], beta_0=derivative_estimator_1$estimator[iBest,])
    VF_hat2 = compute_m(X_obs_unrolled = X_obs_estimator, x_eval = X_obs, beta=derivative_obs_2$estimator, m_0=m20, x0=x[iBest,], beta_0=derivative_estimator_2$estimator[iBest,])

    YObs = aperm(array(cbind(Y1,Y2),dim = c(nT-1,nObs,2)),c(2,3,1))
    VFObs = aperm(array(cbind(VF_hat1,VF_hat2),dim = c(nT-1,nObs,2)),c(2,3,1))

    alpha_i <- NULL
    gamma_t <- NULL

    if (FE) {
        alpha_i = apply(YObs - VFObs, MARGIN = c(1, 2), FUN = sum)/(nT-1)
    }
    if (TE) {
        gamma_t = t(apply(YObs - VFObs, MARGIN = c(2, 3), FUN = sum)/nObs)
    }
    return(listN(alpha_i, gamma_t))
}