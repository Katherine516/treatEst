library(testthat)

test_that("ATE weights return correct length, are numeric, and sum to 1 within groups", {
  set.seed(123)
  n <- 100
  p <- 5

  X <- matrix(rnorm(n * p), nrow = n)
  A <- rbinom(n, 1, 0.5)

  # Simulate outcome model for m_a_hat
  m0_true <- X %*% c(1, -1, 0.5, 0, 0)
  m1_true <- X %*% c(2, 0, -0.5, 0, 0)
  sigma <- 1
  Y0 <- m0_true + rnorm(n, 0, sigma)
  Y1 <- m1_true + rnorm(n, 0, sigma)
  Y <- ifelse(A == 1, Y1, Y0)

  m0_hat <- m0_true + rnorm(n, 0, 0.1)
  m1_hat <- m1_true + rnorm(n, 0, 0.1)

  lambda <- 0.1
  sigma2 <- sigma^2

  gamma1 <- estimate_ate_weights(A, X, a = 1, kernel_func = rbf_kernel, lambda = lambda, sigma2 = sigma2, sigma = sigma)
  gamma0 <- estimate_ate_weights(A, X, a = 0, kernel_func = rbf_kernel, lambda = lambda, sigma2 = sigma2, sigma = sigma)

  # Check types and lengths
  expect_type(gamma1, "double")
  expect_length(gamma1, n)
  expect_type(gamma0, "double")
  expect_length(gamma0, n)

  # Check weights sum to 1 within groups (allow small tolerance)
  expect_equal(sum(gamma1[A == 1]), 1, tolerance = 1e-6)
  expect_equal(sum(gamma0[A == 0]), 1, tolerance = 1e-6)

  # Check ATE estimate is numeric and reasonable
  ate_estimate <- estimate_ate(Y, A, m1_hat, m0_hat, gamma1, gamma0)
  expect_type(ate_estimate, "double")
  expect_length(ate_estimate, 1)
  expect_true(abs(ate_estimate) < 10)  # sanity check
})

