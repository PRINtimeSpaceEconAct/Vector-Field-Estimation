# Let's create a classic ill-conditioned matrix: the Hilbert matrix
hilbert <- function(n) {
    i <- 1:n
    j <- 1:n
    1 / (outer(i - 1, j, "+"))
}

A <- hilbert(10)

# The condition number will be very large
# kappa() is a good estimate
kappa(A)

set.seed(42)
n <- 10
A <- hilbert(n)
x_true <- 1:n
b <- A %*% x_true

# This will likely throw a warning
x_naive <- solve(A, b)
#> Warning in solve.default(A, b): system is computationally singular:
#> reciprocal condition number = 6.2402e-14

# Compare with the true solution
print(x_naive)
#>  [1]  1.000000  1.999999  3.000006  3.999978  5.000050  5.999912  7.000096
#>  [8]  7.999914  9.000030 10.000010

# Calculate the error
sum((x_naive - x_true)^2)
#> [1] 2.193136e-08

svd_A <- svd(A)
U <- svd_A$u
D <- diag(svd_A$d)
V <- svd_A$v

# The singular values show the problem - they get incredibly small
print(svd_A$d)
#>  [1] 1.664468e+00 2.871578e-01 2.305413e-02 1.134737e-03 3.655890e-05
#>  [6] 7.978254e-07 1.171866e-08 1.116035e-10 6.641618e-13 2.302307e-15

# Compute the pseudoinverse of A by inverting non-tiny singular values
# Set a tolerance to avoid dividing by nearly zero
tol <- .Machine$double.eps
d_inv <- 1 / svd_A$d
d_inv[svd_A$d < tol] <- 0 # Crucial step: treat tiny values as zero
D_inv <- diag(d_inv)

# A+ = V D+ U'
A_pinv <- V %*% D_inv %*% t(U)

# The solution is x = A+ b
x_svd <- A_pinv %*% b

# Compare with the true solution
print(x_svd)
#>            [,1]
#>  [1,]  1.000000
#>  [2,]  2.000000
#>  [3,]  3.000000
#>  [4,]  4.000000
#>  [5,]  5.000000
#>  [6,]  6.000000
#>  [7,]  7.000000
#>  [8,]  8.000000
#>  [9,]  9.000000
#> [10,] 10.000000

# Calculate the error
sum((x_svd - x_true)^2)

x_MASS = MASS::ginv(A) %*% b
sum((x_MASS - x_true)^2)

# x_qr <- qr.solve(A, b)
x_qr <- solve(qr(A,LAPACK=TRUE,tol=1e-10),b)

# Compare with the true solution
print(x_qr)
#>  [1]  1.000000e+00  2.000000e+00  3.000000e+00  4.000000e+00  5.000000e+00
#>  [6]  6.000000e+00  7.000000e+00  8.000000e+00  9.000000e+00 1.000000e+01

# Calculate the error
sum((x_qr - x_true)^2)
#> [1] 1.053805e-24

