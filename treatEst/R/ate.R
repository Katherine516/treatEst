library(slam)
library(gurobi)
library(Matrix)

#' RBF Kernel Function
#' @param x,y Vectors
#' @param sigma Kernel bandwidth
rbf_kernel <- function(x, y, sigma = 1, ...) {
  exp(-sum((x - y)^2) / (2 * sigma^2))
}


#' Linear Kernel Function
linear_kernel <- function(x, y, ...) {
  sum(x * y)
}

#' Polynomial Kernel Function
#' @param degree Polynomial degree
#' @param coef0 Constant term
polynomial_kernel <- function(x, y, degree = 3, coef0 = 1, ...) {
  (sum(x * y) + coef0)^degree
}

#' Gaussian Kernel Function
gaussian_kernel <- function(x, y, sigma = 1, ...) {
  d2 <- sum((x - y)^2)
  exp(-d2 / (2 * sigma^2))
}

#' Laplacian Kernel Function
laplacian_kernel <- function(x, y, sigma = 1, ...) {
  d1 <- sum(abs(x - y))
  exp(-d1 / sigma)
}

#' Compute Kernel Gram Matrix
#' @param X Covariate matrix
#' @param kernel_func Kernel function
#' @param ... Additional kernel arguments
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



#' Estimate ATE Weights using Kernel Balancing
#'
#' This function estimates balancing weights to minimize the imbalance in RKHS.
#'
#' @param A Treatment assignment vector (0/1)
#' @param X Covariate matrix (n x p)
#' @param a Treatment group of interest (0 or 1)
#' @param kernel_func Kernel function (default = rbf_kernel)
#' @param lambda Regularization parameter
#' @param sigma2 Noise variance for stability (default = 1e-6)
#' @param ... Additional kernel arguments
#' @return Weight vector (length = n)
estimate_ate_weights <- function(A, X, a, kernel_func, lambda = 1, sigma2 = 1e-6, ...) {
  n <- nrow(X)
  idx <- which(A == a)
  if (length(idx) == 0) {
    warning("No samples with A == a found.")
    return(rep(0, n))
  }

  K <- compute_kernel_matrix(X, kernel_func, ...)

  K_sub <- K[idx, idx, drop = FALSE]
  p <- -rowSums(K[idx, ]) / (n^2)

  Q <- (K_sub / n^2) + diag(lambda * sigma2 / n^2, length(idx))
  Q <- (Q + t(Q)) / 2
  Q <- Q / 2

  Q_triplet <- slam::as.simple_triplet_matrix(Q)

  model <- list(
    Q = Q_triplet,
    obj = p,
    modelsense = "min",
    A = matrix(1, nrow = 1, ncol = length(idx)),
    rhs = 1,
    sense = "=",
    lb = rep(0, length(idx))
  )

  params <- list(OutputFlag = 0)
  result <- gurobi::gurobi(model, params)

  gamma <- numeric(n)
  gamma[idx] <- result$x
  return(gamma)
}



# AIPW estimator for ATE at treatment a
estimate_ate_aipw <- function(Y, A, a, m_a_hat, gamma_a) {
  gamma_weighted <- gamma_a * (A == a)
  return(mean(m_a_hat) - mean(gamma_weighted * (m_a_hat - Y)))
}



#' Estimate Average Treatment Effect
#'
#' @param Y Outcome vector
#' @param A Treatment vector (0/1)
#' @param X Covariate matrix
#' @param kernel_func Kernel function
#' @param lambda Regularization parameter
#' @param sigma2 Small noise term
#' @param ... Additional kernel arguments
#' @return Estimated ATE value
estimate_ate <- function(Y, A, m1_hat, m0_hat, gamma1, gamma0) {
  psi1 <- estimate_ate_aipw(Y, A, 1, m1_hat, gamma1)
  psi0 <- estimate_ate_aipw(Y, A, 0, m0_hat, gamma0)
  return(psi1 - psi0)
}

