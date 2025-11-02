library(testthat)

test_that("QTE weights return correct length and are numeric", {
  set.seed(123)
  n <- 100
  p <- 5

  X <- matrix(rnorm(n * p), nrow = n)
  A <- rbinom(n, 1, 0.5)
  a <- 1

  # Compute kernel matrix K
  K <- compute_kernel_matrix(X, rbf_kernel, sigma = 1)

  gamma <- estimate_qte_weights(A, X, a = a, K = K)

  expect_type(gamma, "double")
  expect_length(gamma, n)
})
