library(slam)
library(gurobi)
library(Matrix)

# --- Kernel functions ---
rbf_kernel <- function(x, y, sigma = 1) {
  exp(-sum((x - y)^2) / (2 * sigma^2))
}

linear_kernel <- function(x, y) {
  sum(x * y)
}

polynomial_kernel <- function(x, y, degree = 3, coef0 = 1) {
  (sum(x * y) + coef0)^degree
}

gaussian_kernel <- function(x, y, sigma = 1, ...) {
  d2 <- sum((x - y)^2)
  exp(-d2 / (2 * sigma^2))
}

laplacian_kernel <- function(x, y, sigma = 1) {
  exp(-sum(abs(x - y)) / sigma)
}

#   Compute Kernel Gram Matrix (Ka)   #

compute_kernel_matrix <- function(X, kernel_func, ...) {
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.numeric(X)) stop("Matrix X must be numeric")

  n <- nrow(X)
  K <- matrix(0, n, n)

  for (i in 1:n) {
    for (j in i:n) {
      if (identical(kernel_func, linear_kernel)) {
        val <- kernel_func(X[i, ], X[j, ])
      } else {
        val <- kernel_func(X[i, ], X[j, ], ...)
      }
      K[i, j] <- val
      K[j, i] <- val
    }
  }

  return(K)
}


#' Estimate CATE weights via RKHS kernel balancing
#'
#' @param A Treatment vector (0/1)
#' @param V Subgroup label vector
#' @param S Source population indicator (1 = source)
#' @param X Covariates
#' @param a Treatment of interest (0/1)
#' @param v Subgroup value of interest
#' @param kernel_func Kernel function (default = rbf)
#' @param lambda Regularization strength
#' @param sigma2 Small constant for numerical stability
#' @return Weight vector gamma of length n
estimate_cate_weights <- function(A, V, S, X, a, v,
                                  kernel_func,
                                  lambda = 1,
                                  sigma2 = 1e-6,
                                  ...) {
  n <- nrow(X)
  I_av <- as.numeric(A == a & V == v)
  e_vs <- as.numeric(V == v & S == 1)

  K <- compute_kernel_matrix(X, kernel_func, ...)
  I_diag <- diag(I_av)
  Q <- 2 * (I_diag %*% K %*% I_diag + lambda * sigma2 * diag(n))
  Q <- (Q + t(Q)) / 2  # ensure symmetry

  b <- 2 * as.vector(t(e_vs) %*% K %*% I_diag)

  model <- list(
    Q = slam::as.simple_triplet_matrix(Q / 2),
    obj = -b,
    modelsense = "min",
    A = matrix(I_av, nrow = 1),
    rhs = 1,
    sense = "=",
    lb = rep(0, n)  # non-negative weights
  )

  result <- gurobi::gurobi(model, params = list(OutputFlag = 0))

  if (is.null(result$x)) {
    warning("Gurobi failed — using fallback uniform weights.")
    gamma <- I_av / sum(I_av)
  } else {
    gamma <- result$x
    gamma[I_av == 0] <- 0
    gamma <- gamma / sum(gamma)
  }

  return(gamma)
}







#' Estimate Conditional Average Treatment Effect (CATE)
#'
#' @param Y Outcome vector
#' @param A Treatment vector (0/1)
#' @param V Subgroup label vector
#' @param S Source indicator (1 = source)
#' @param X Covariate matrix
#' @param a Treatment group of interest (0/1)
#' @param v Subgroup value
#' @param kernel_func Kernel function
#' @param lambda Regularization parameter
#' @param sigma2 Small constant for numerical stability
#' @return Estimated CATE
estimate_cate <- function(Y, A, S, V, v, m1_hat, m0_hat, gamma1, gamma0) {
  # Proportion of source units in subgroup
  a <- mean((V == v) & (S == 1))
  if (a == 0) {
    warning("No source units with V == ", v)
    return(NA_real_)
  }
  alpha <- 1 / a

  # weights are numeric vectors
  gamma1 <- as.numeric(gamma1)
  gamma0 <- as.numeric(gamma0)

  # Influence-function terms (vectorized)
  T1 <- ((A == 1) & (V == v)) * gamma1 * (Y - m1_hat)
  T2 <- ((A == 0) & (V == v)) * gamma0 * (Y - m0_hat)
  T3 <- ((V == v) & (S == 1)) * (m1_hat - m0_hat)

  # CATE estimate
  cate <- mean(alpha * (T1 - T2 + T3))
  return(cate)
}



