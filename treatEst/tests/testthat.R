# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

library(testthat)
library(treatEst)
test_check("treatEst")

devtools::document()     # Generates NAMESPACE & man/
devtools::load_all()     # Load package without installing
devtools::check()        # Run full CRAN check
devtools::test()         # Run all testthat tests

test_that("ATE weights return correct length", {
  n <- 100
  dta <- data_generation(n, correct = 1, te_function = function(x, y) 1, sigma = 1)
  X <- as.matrix(dta[, c("Z1", "Z2")])
  gamma <- estimate_ate_weights(A = dta$A, X = X)
  expect_length(gamma, n)
})

test_that("QTE estimator returns numeric", {
  n <- 100
  dta <- data_generation(n, correct = 1, te_function = function(x, y) 0.5, sigma = 1)
  X <- as.matrix(dta[, c("Z1", "Z2")])
  gamma1 <- estimate_qte_weights(dta$A, X, a = 1)
  gamma0 <- estimate_qte_weights(dta$A, X, a = 0)
  est <- estimate_qte(Y = dta$Y, A = dta$A, tau = 0.5,
                      m1_hat = dta$Y1, m0_hat = dta$Y0,
                      gamma1 = gamma1, gamma0 = gamma0)
  expect_type(est, "double")
})


