make_ssl_fixture_data <- function(n = 100, p = 20, seed = 20260902) {
  set.seed(seed)
  x <- matrix(rnorm(n * p), n, p)
  colnames(x) <- paste0("V", seq_len(p))
  beta_true <- c(1, -0.5, 0.8, -0.3, 0.5, rep(0, p - 5))
  lp <- as.numeric(x %*% beta_true)
  time <- round(rexp(n, rate = exp(lp - mean(lp))) * 5, 4) + 0.01
  status <- sample(c(0L, 1L, 2L), n, replace = TRUE, prob = c(0.25, 0.55, 0.20))
  y <- cbind(time = time, status = status)
  list(x = x, y = y)
}
