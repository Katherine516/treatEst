library(testthat)

test_that("CATE weights return correct length and are numeric", {
  set.seed(123)
  n <- 100
  p <- 5

  X <- matrix(rnorm(n * p), nrow = n)
  A <- rbinom(n, 1, 0.5)
  V <- sample(1:3, n, replace = TRUE)
  S <- rbinom(n, 1, 0.7)

  a <- 1
  v <- 1

  # Compute a data-driven sigma: median pairwise distance
  pairwise_dist <- dist(X)  # Euclidean distance by default
  sigma_default <- median(pairwise_dist)

  gamma <- estimate_cate_weights(
    A, V, S, X,
    a = a, v = v,
    kernel_func = rbf_kernel,
    sigma = sigma_default
  )

  expect_type(gamma, "double")
  expect_length(gamma, n)
})

