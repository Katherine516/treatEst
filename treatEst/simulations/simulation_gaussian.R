library(quantreg)
library(gurobi)

source("R/synthetic_data.R")

# --- Kernels ---
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

# --- Adaptive bandwidth for kernels ---
adaptive_sigma <- function(X) {
  as.numeric(median(dist(X)))
}

# --- Conditional mean estimates using linear regression ---
estimate_conditional_means <- function(Y, A, X) {
  m1 <- lm(Y ~ ., data = data.frame(Y = Y[A == 1], X[A == 1, , drop = FALSE]))
  m0 <- lm(Y ~ ., data = data.frame(Y = Y[A == 0], X[A == 0, , drop = FALSE]))

  pred_df <- as.data.frame(X)
  list(
    m1_hat = predict(m1, newdata = pred_df),
    m0_hat = predict(m0, newdata = pred_df)
  )
}

# --- Conditional quantile estimates using quantile regression ---
estimate_conditional_quantiles <- function(Y, A, X, tau) {
  m1 <- rq(Y ~ ., tau = tau, data = data.frame(Y = Y[A == 1], X[A == 1, , drop = FALSE]))
  m0 <- rq(Y ~ ., tau = tau, data = data.frame(Y = Y[A == 0], X[A == 0, , drop = FALSE]))

  pred_df <- as.data.frame(X)
  list(
    m1_hat = predict(m1, newdata = pred_df),
    m0_hat = predict(m0, newdata = pred_df)
  )
}

# --- Evaluation metrics ---
calc_metrics <- function(estimates, truths) {
  bias <- mean(estimates - truths)
  rmse <- sqrt(mean((estimates - truths)^2))
  list(bias = bias, rmse = rmse)
}

# --- Main simulation parameters ---
n_sim <- 100
n <- 1000
te_true <- 2
tau_vec <- c(0.25, 0.5, 0.75)

# Storage for results
ate_estimates <- numeric(n_sim)
cate_estimates <- numeric(n_sim)
qte_estimates <- matrix(NA_real_, n_sim, length(tau_vec))
colnames(qte_estimates) <- paste0("QTE_tau_", tau_vec)

# Select kernel and parameters here
kernel_to_use <- gaussian_kernel
kernel_params <- list(sigma = NULL)  # sigma will be adaptive

for (i in seq_len(n_sim)) {
  set.seed(123 + i)

  # Generate data for ATE and QTE
  dta <- data_generation(n = n, correct = 1, te_function = NULL, sigma = 1)
  dta$te <- te_true
  dta$Y1 <- dta$Y0 + dta$te
  dta$Y <- ifelse(dta$A == 1, dta$Y1, dta$Y0)

  X <- as.matrix(dta[, c("X1", "X2")])
  sigma <- adaptive_sigma(X)
  if (is.null(kernel_params$sigma)) {
    kernel_params$sigma <- sigma
  }

  # Estimate conditional means
  cond_means <- estimate_conditional_means(dta$Y, dta$A, X)

  # Estimate weights for A=1 and A=0 groups
  gamma1 <- estimate_ate_weights(dta$A, X, a = 1,
                                 kernel = kernel_to_use,
                                 sigma2 = sigma,
                                 lambda = 1,
                                 sigma = kernel_params$sigma)
  gamma0 <- estimate_ate_weights(dta$A, X, a = 0,
                                 kernel = kernel_to_use,
                                 sigma2 = sigma,
                                 lambda = 1,
                                 sigma = kernel_params$sigma)

  # Estimate ATE
  ate_estimates[i] <- estimate_ate(dta$Y, dta$A,
                                   cond_means$m1_hat,
                                   cond_means$m0_hat,
                                   gamma1,
                                   gamma0)

  # Generate data for CATE
  dta2 <- data_genera(n = n,
                      treffect_A = 2,
                      beta_W_Y = 1, beta_W_S = 1, beta_V_S = 1,
                      beta_S_Y = 1, beta_W_A = 1, beta_V_A = 1,
                      beta_V_Y = 1, heteroTE = TRUE, probS = 0.5)

  X2 <- as.matrix(dta2[, c("W", "V")])
  sigma_cate <- adaptive_sigma(X2)
  if (is.null(kernel_params$sigma)) {
    kernel_params$sigma <- sigma_cate
  }

  # Estimate conditional means for CATE
  cond_means_cate <- estimate_conditional_means(dta2$Y, dta2$A, X2)

  # Estimate weights for CATE
  gamma1_cate <- estimate_cate_weights(dta2$A, dta2$V, dta2$S, X2,
                                       a = 1, v = 1,
                                       kernel = kernel_to_use,
                                       sigma2 = sigma_cate,
                                       lambda = 1,
                                       sigma = kernel_params$sigma)
  gamma0_cate <- estimate_cate_weights(dta2$A, dta2$V, dta2$S, X2,
                                       a = 0, v = 1,
                                       kernel = kernel_to_use,
                                       sigma2 = sigma_cate,
                                       lambda = 1,
                                       sigma = kernel_params$sigma)

  cate_estimates[i] <- estimate_cate(Y=dta2$Y, A=dta2$A, S=dta2$S, V=dta2$V, v = 1,
                                     m1_hat = cond_means_cate$m1_hat,
                                     m0_hat = cond_means_cate$m0_hat,
                                     gamma1 = gamma1_cate,
                                     gamma0 = gamma0_cate)

  # QTE estimation loop
  for (j in seq_along(tau_vec)) {
    tau <- tau_vec[j]

    cond_quants <- estimate_conditional_quantiles(dta$Y, dta$A, X, tau)

    K <- compute_kernel_matrix(X, kernel_to_use, sigma = kernel_params$sigma)

    gamma1_qte <- estimate_qte_weights(dta$A, X, a = 1,
                                       K = K,
                                       sigma2 = sigma,
                                       lambda = 1)
    gamma0_qte <- estimate_qte_weights(dta$A, X, a = 0,
                                       K = K,
                                       sigma2 = sigma,
                                       lambda = 1)

    qte_estimates[i, j] <- estimate_qte(dta$Y, dta$A, tau,
                                        cond_quants$m1_hat,
                                        cond_quants$m0_hat,
                                        gamma1_qte,
                                        gamma0_qte)
  }
}

# --- Evaluation summaries ---
cat("===== ATE Results =====\n")
ate_metrics <- calc_metrics(ate_estimates, rep(te_true, n_sim))
cat(sprintf("Bias: %.4f, RMSE: %.4f\n\n", ate_metrics$bias, ate_metrics$rmse))

cat("===== CATE Results =====\n")
heteroTE <- TRUE   # set to FALSE if you want homogeneous TE
cate_true <- if (heteroTE) te_true + 1 else te_true
cate_metrics <- calc_metrics(cate_estimates, rep(cate_true, n_sim))
cat(sprintf("Bias: %.4f, RMSE: %.4f\n\n", cate_metrics$bias, cate_metrics$rmse))

cat("===== QTE Results =====\n")
for (j in seq_along(tau_vec)) {
  qte_metrics <- calc_metrics(qte_estimates[, j], rep(te_true, n_sim))
  cat(sprintf("Tau %.2f — Bias: %.4f, RMSE: %.4f\n",
              tau_vec[j], qte_metrics$bias, qte_metrics$rmse))
}



