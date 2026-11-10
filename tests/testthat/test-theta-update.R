test_that("update_mixture_prob_count matches (a + nnz) / (a + b + p)", {
  betas <- c(1.2, 0, -0.5, 0, 0.3)  # 3 nonzero, p = 5
  p <- length(betas)

  # (a, b) = (1, 1): (1 + 3) / (1 + 1 + 5) = 4/7
  expect_equal(update_mixture_prob_count(betas, a = 1, b = 1),
               (1 + 3) / (1 + 1 + 5), tolerance = 1e-14)

  # (a, b) = (1, p): (1 + 3) / (1 + 5 + 5) = 4/11
  expect_equal(update_mixture_prob_count(betas, a = 1, b = p),
               (1 + 3) / (1 + p + p), tolerance = 1e-14)
})

test_that("expected_inclusion_probs with theta matches Bayes formula directly", {
  betas <- c(0.5, 0.0, -0.3, 0.0, 0.1)
  s0 <- 0.04; s1 <- 0.5; theta <- 0.3

  pips_got <- expected_inclusion_probs(s0, s1, pips = rep(0.1, 5), betas, theta = theta)

  d1 <- dlaplace(betas, mu = 0, b = s1)
  d0 <- dlaplace(betas, mu = 0, b = s0)
  pips_expected <- theta * d1 / (theta * d1 + (1 - theta) * d0)

  expect_equal(pips_got, pips_expected, tolerance = 1e-14)
})

test_that("expected_inclusion_probs errors when theta is outside (0, 1)", {
  betas <- c(0.1, -0.2, 0.0)
  pips  <- rep(0.1, 3)
  expect_error(
    expected_inclusion_probs(0.04, 0.5, pips, betas, theta = 0),
    regexp = "strictly between 0 and 1"
  )
  expect_error(
    expected_inclusion_probs(0.04, 0.5, pips, betas, theta = 1),
    regexp = "strictly between 0 and 1"
  )
  expect_error(
    expected_inclusion_probs(0.04, 0.5, pips, betas, theta = 1.5),
    regexp = "strictly between 0 and 1"
  )
})

test_that("count-rule fit returns correct theta fields", {
  dat <- make_ssl_fixture_data()
  x <- dat$x; y <- dat$y
  theta_a <- 1; theta_b <- ncol(x)
  initial_sparsity <- 0.05

  fit <- fit_ssl_psdh(x, y,
    ss = c(0.04, 0.5), initial_sparsity = initial_sparsity,
    theta_a = theta_a, theta_b = theta_b,
    theta_update = "count",
    init_method = "LASSO_cv")

  expect_equal(fit$theta_update, "count")
  expect_equal(fit$theta_path[1], initial_sparsity)
  expect_equal(fit$theta,
               update_mixture_prob_count(fit$coefficients$Estimate,
                                         theta_a, theta_b),
               tolerance = 1e-12)
})

test_that("default initial_sparsity is 0.5", {
  expect_equal(formals(fit_ssl_psdh)$initial_sparsity, 0.5)
})

test_that("theta_update = 'expected' returns correct theta_update field", {
  dat <- make_ssl_fixture_data()
  fit <- fit_ssl_psdh(dat$x, dat$y,
    ss = c(0.04, 0.5), initial_sparsity = 0.05,
    theta_update = "expected",
    init_method = "LASSO_cv")

  expect_equal(fit$theta_update, "expected")
})
