library(slam)
library(gurobi)
#library(reticulate)
#source_python("gp_simu_gate.py")

# Kernel functions
rbf_kernel <- function(x, y, sigma = 1) {
  exp(-sum((x - y)^2) / (2 * sigma^2))
}

linear_kernel <- function(x, y) {
  sum(x * y)
}

polynomial_kernel <- function(x, y, degree = 3, coef0 = 1) {
  (sum(x * y) + coef0)^degree
}

gaussian_kernel <- function(x, y, sigma = 1) {
  d2 <- sum((x - y)^2)
  exp(-d2 / (2 * sigma^2))
}

laplacian_kernel <- function(x, y, sigma = 1) {
  d1 <- sum(abs(x - y))
  exp(-d1 / sigma)
}


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


#' Estimate kernel balancing weights for QTE using Gurobi
#'
#' This function minimizes the worst-case imbalance in a RKHS subject to a regularization penalty.
#'
#' @param A Treatment assignment vector (0/1)
#' @param X Covariate matrix (rows = samples)
#' @param a Treatment group of interest (0 or 1)
#' @param tau Quantile level (0 < tau < 1)
#' @param kernel_func Kernel function (default is RBF)
#' @param lambda Regularization parameter controlling weight complexity (default = 1)
#' @param sigma2 Noise variance parameter (default = 1)
#'
#' @return A vector of balancing weights gamma
estimate_qte_weights <- function(A, X, a, K, lambda = 1, sigma2 = 1e-6) {
  if (!is.matrix(X)) X <- as.matrix(X)
  n <- nrow(X)
  idx <- which(A == a)

  if (length(idx) == 0) {
    warning("No samples with A == ", a)
    return(rep(0, n))
  }

  # Submatrix of K for treated units
  K_sub <- K[idx, idx, drop = FALSE]

  # Linear term
  e_n <- rep(1, n)
  p <- -as.vector(K[idx, ] %*% e_n) / (n^2)

  # Quadratic term
  Q <- K_sub / n^2 + diag(lambda * sigma2 / n^2, length(idx))
  Q <- (Q + t(Q)) / 2
  Q <- Q / 2  # Gurobi expects 1/2 x'Qx

  Q_triplet <- slam::as.simple_triplet_matrix(Q)

  model <- list(
    Q = Q_triplet,
    obj = p,
    modelsense = "min",
    A = matrix(1, nrow = 1, ncol = length(idx)),
    rhs = 1,
    sense = "=",
    lb = rep(0, length(idx))  # Lower bound only
  )

  result <- gurobi::gurobi(model, params = list(OutputFlag = 0))

  # Construct full gamma vector
  gamma <- numeric(n)
  gamma[idx] <- result$x

  return(gamma)
}


#' AIPW Estimator for QTE at quantile tau
#'
#' @param Y Outcome vector
#' @param A Treatment vector (0/1)
#' @param a Treatment arm of interest
#' @param tau Quantile level (0 < tau < 1)
#' @param m_a_hat Estimated conditional quantile function at tau
#' @param gamma_a Balancing weights for A == a group
#'
#' @return Quantile estimate for treatment group a at level tau
estimate_qte_aipw <- function(Y, A, a, tau, m_a_hat, gamma_a) {
  gamma_weighted <- gamma_a * as.numeric(A == a)
  aug_term <- gamma_weighted * (as.numeric(Y <= m_a_hat) - tau)
  est <- mean(aug_term) + mean(m_a_hat)  # This should be scalar
  return(est)
}




#' Estimate Quantile Treatment Effect (QTE)
#'
#' @param Y Outcome vector
#' @param A Treatment vector (0/1)
#' @param tau Quantile level (0 < tau < 1)
#' @param m1_hat Estimated quantiles for A=1
#' @param m0_hat Estimated quantiles for A=0
#' @param gamma1 Weights for A=1
#' @param gamma0 Weights for A=0
#'
#' @return Estimated QTE (numeric scalar)
estimate_qte <- function(Y, A, tau, m1_hat, m0_hat, gamma1, gamma0) {
  q1 <- estimate_qte_aipw(Y, A, 1, tau, m1_hat, gamma1)
  q0 <- estimate_qte_aipw(Y, A, 0, tau, m0_hat, gamma0)
  return(as.numeric(q1) - as.numeric(q0))  # Ensure scalar subtraction
}



