# Verifies the RcppArmadillo CV-metric kernels in src/cv_fastcrrp.cpp
# reproduce their pure-R references. These tests require the compiled package
# (run devtools::document() then devtools::load_all() / install first).

make_cr_data <- function(n = 120, p = 4, seed = 1) {
  set.seed(seed)
  x      <- matrix(rnorm(n * p), n, p)
  beta   <- c(1, -0.5, 0, 0.25)[seq_len(p)]
  lp     <- as.numeric(x %*% beta)
  time   <- round(rexp(n, rate = exp(lp - mean(lp))) * 10, 4) + 0.01
  # status: 0 censored, 1 cause-1, 2 competing; keep times distinct-ish
  status <- sample(c(0L, 1L, 2L), n, replace = TRUE, prob = c(0.3, 0.5, 0.2))
  list(x = x, time = time, status = status, lp = lp)
}

test_that("cpp_km_censoring matches survival::survfit censoring KM", {
  skip_if_not_installed("survival")
  d <- make_cr_data()
  cpp <- cpp_km_censoring(d$time, d$status)

  cens_status <- ifelse(d$status == 0, 1, 0)
  km <- survival::survfit(survival::Surv(d$time, cens_status) ~ 1)

  get_G <- function(t) {
    idx <- findInterval(t, km$time)
    probs <- rep(1, length(t)); v <- idx > 0
    probs[v] <- km$surv[idx[v]]; probs
  }
  idx_m <- findInterval(d$time - 1e-08, km$time)
  G_minus <- rep(1, length(d$time)); vm <- idx_m > 0
  G_minus[vm] <- km$surv[idx_m[vm]]

  expect_equal(cpp$G_Ti,       get_G(d$time), tolerance = 1e-10)
  expect_equal(cpp$G_Ti_minus, G_minus,       tolerance = 1e-10)
})

test_that("cpp_cv_cindex 'wolbers' matches measure_ssl_psdh", {
  d <- make_cr_data()
  coef <- matrix(c(1, -0.5, 0, 0.25,
                   0.2, 0.2, 0.2, 0.2), nrow = 4)
  eval_time <- as.numeric(quantile(d$time, 0.5))

  got <- cpp_cv_cindex(d$x, coef, d$time, d$status,
                       tuning = "wolbers",
                       evaluation_time = eval_time)

  ref <- apply(d$x %*% coef, 2, function(rs) {
    measure_ssl_psdh(y_true = cbind(d$time, d$status),
                     risk_score = rs, evaluation_time = eval_time)
  })

  expect_equal(got, ref, tolerance = 1e-10)
})

test_that("cpp_cv_cindex 'normal' matches survival::concordance (distinct times)", {
  skip_if_not_installed("survival")
  d <- make_cr_data()
  # de-duplicate times so Harrell tie handling cannot diverge from survival
  d$time <- d$time + seq_along(d$time) * 1e-6
  coef <- matrix(c(1, -0.5, 0, 0.25), nrow = 4)

  got <- cpp_cv_cindex(d$x, coef, d$time, d$status,
                       tuning = "normal", evaluation_time = 0,
                       reverse = FALSE)

  rs  <- as.numeric(d$x %*% coef)
  ref <- survival::concordance(
    survival::Surv(d$time, d$status == 1) ~ rs
  )$concordance

  expect_equal(got[1], ref, tolerance = 1e-8)
})
