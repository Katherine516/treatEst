# for ATE and QTE (ATE and QTE are just the mean and quantile difference between Y1 and Y0)
# use these R functions for generating synthetic data:

# Data generation function

data_generation <- function(n, correct, te_function, sigma) {
  # Data
  X1 <- rnorm(n, 0, sigma)
  X2 <- rnorm(n, 0, sigma)

  prt <- 1 / (1 + exp(-(-0.5 + 0.5 * (X1 + X2))))

  A <- rbinom(n, 1, prt)
  Y0 <- 0.5 * (X1 + X2) + rnorm(n, 0, sigma)
  # Apply the quantile-specific treatment effect with predefined quantiles

  Y1 <- Y0 + te
  Y <- Y0 * (1 - A) + Y1 * A
  dta <- data.frame(Y, X1, X2, A)

  wrong2 <- log(abs(dta$X2))
  wrong1 <- (dta$X2) / (exp(dta$X1))

  Z1 <- correct * dta$X1 + (1 - correct) * (wrong1)
  Z2 <- correct * dta$X2 + (1 - correct) * (wrong2)
  Z <- cbind(Z1, Z2)

  colnames(Z) <- c("Z1", "Z2")

  dta <- data.frame(dta, Z, Y1, Y0)
  return(dta)

}



# For CATE:

data_genera <- function(n,
                        beta_W_A = 1, beta_V_A = 1,
                        beta_W_Y = 1, beta_V_Y = 1,
                        beta_W_S = 0, beta_V_S = 0, beta_S_Y = 0,
                        treffect_A = 2, heteroTE = TRUE, probS = NULL){


  #W = runif(n,0,1)
  W = rnorm(n, 0,1)
  V = rbinom(n, 1, 0.5)
  Z <- (W)/(1 + W) + 2


  if (is.null(probS)) {
    linear_probS <- beta_W_S * W + beta_V_S * V + beta_S_Y
    probS <- 1 / (1 + exp(-linear_probS))
  }
  S <- rbinom(n, 1, probS)

  nS1 = sum(S)
  nS0 = n - nS1
  AS1 = rbinom(nS1, 1, 0.5)

  WS0 = W[S==0]
  VS0 = V[S==0]
  #prAS0 = (1+exp( -(-0.1*beta_W_A + beta_W_A*WS0 + beta_V_A*VS0) ))^-1
  prAS0 <- 1 / (1 + exp(-(beta_W_A * WS0 + beta_V_A * VS0)))
  AS0 <- rbinom(nS0, 1, prAS0)

  A = rep(NA,n)
  A[S==1] = AS1
  A[S==0] = AS0

  Y0 =  beta_W_Y * W + beta_V_Y * V + rnorm(n)
  if (heteroTE) {
    Y1 <- Y0 + V + treffect_A
  } else {
    Y1 <- Y0 + treffect_A
  }

  # Observed outcome
  Y <- ifelse(A == 1, Y1, Y0)

  data <- data.frame(Y, Y1, Y0, W, S, A, V, Z)

  return(data)

}



