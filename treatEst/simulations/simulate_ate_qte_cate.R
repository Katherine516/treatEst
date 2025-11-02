library(gbm)        # for CATE modeling
library(MASS)       # for quantile regression
library(gurobi)     # kernel balancing
library(quantreg)   # for QTE
library(dplyr)

# Simulation settings
set.seed(123)
n_sim <- 100      # number of repetitions
n <- 500          # sample size
sigma <- 1
tau <- 0.5        # quantile for QTE
correct <- 1      # use correct model
tr_effectA <- 2   # treatment effect (for ATE/CATE/QTE)

# Store results
ate_estimates <- numeric(n_sim)
ate_truths    <- numeric(n_sim)

qte_estimates <- numeric(n_sim)
qte_truths    <- numeric(n_sim)

cate_biases   <- numeric(n_sim)

### Simulation loop
for (i in 1:n_sim) {
  ## ---- ATE / QTE ----
  dta <- data_generation(n, correct, tr_effectA, sigma)

  # ATE
  ate_est <- mean(dta$Y[dta$A == 1]) - mean(dta$Y[dta$A == 0])
  true_ate <- mean(dta$Y1) - mean(dta$Y0)

  ate_estimates[i] <- ate_est
  ate_truths[i] <- true_ate

  # QTE (simple plugin)
  qte_est <- quantile(dta$Y[dta$A == 1], tau) - quantile(dta$Y[dta$A == 0], tau)
  true_qte <- quantile(dta$Y1, tau) - quantile(dta$Y0, tau)

  qte_estimates[i] <- qte_est
  qte_truths[i] <- true_qte

  ## ---- CATE ----
  dta_cate <- data_genera(n, tr_effectA,
                          beta_W_Y = 1, beta_W_S = 0.5,
                          beta_V_S = 1, beta_S_Y = 1,
                          beta_W_A = 1, beta_V_A = 1, beta_V_Y = 1,
                          heteroTE = TRUE, probS = 0.5)

  # Fit model (A, W, S, V → Y)
  fit <- gbm(Y ~ A * (W + S + V), data = dta_cate, distribution = "gaussian", n.trees = 100)

  # Predicted CATE
  pred1 <- predict(fit, newdata = data.frame(A = 1, W = dta_cate$W, S = dta_cate$S, V = dta_cate$V), n.trees = 100)
  pred0 <- predict(fit, newdata = data.frame(A = 0, W = dta_cate$W, S = dta_cate$S, V = dta_cate$V), n.trees = 100)
  pred_cate <- pred1 - pred0
  true_cate <- dta_cate$Y1 - dta_cate$Y0

  cate_biases[i] <- mean(pred_cate - true_cate)
}

### Evaluation
# ATE
ate_bias <- mean(ate_estimates - ate_truths)
ate_mse <- mean((ate_estimates - ate_truths)^2)
ate_ci_coverage <- mean((ate_estimates > ate_truths - 1.96 * sd(ate_estimates)) &
                          (ate_estimates < ate_truths + 1.96 * sd(ate_estimates)))

# QTE
qte_bias <- mean(qte_estimates - qte_truths)
qte_mse <- mean((qte_estimates - qte_truths)^2)
qte_ci_coverage <- mean((qte_estimates > qte_truths - 1.96 * sd(qte_estimates)) &
                          (qte_estimates < qte_truths + 1.96 * sd(qte_estimates)))

# CATE
cate_bias_mean <- mean(cate_biases)
cate_bias_sd <- sd(cate_biases)

### Print results
cat("==== ATE Evaluation ====\n")
cat("Bias:", ate_bias, "\n")
cat("MSE:", ate_mse, "\n")
cat("Coverage (approx):", ate_ci_coverage, "\n\n")

cat("==== QTE Evaluation ====\n")
cat("Bias:", qte_bias, "\n")
cat("MSE:", qte_mse, "\n")
cat("Coverage (approx):", qte_ci_coverage, "\n\n")

cat("==== CATE Evaluation ====\n")
cat("Avg Bias in CATE:", cate_bias_mean, "\n")
cat("SD of CATE Bias:", cate_bias_sd, "\n")
