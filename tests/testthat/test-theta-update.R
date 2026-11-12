test_that("update_mixture_prob('count') matches (a + nnz) / (a + b + p)", {
  betas <- c(1.2, 0, -0.5, 0, 0.3)  # 3 nonzero, p = 5
  p <- length(betas)

  # (a, b) = (1, 1): (1 + 3) / (1 + 1 + 5) = 4/7
  expect_equal(update_mixture_prob(betas = betas, a = 1, b = 1,
                                   theta_update = "count"),
               (1 + 3) / (1 + 1 + 5), tolerance = 1e-14)

  # (a, b) = (1, p): (1 + 3) / (1 + 5 + 5) = 4/11
  expect_equal(update_mixture_prob(betas = betas, a = 1, b = p,
                                   theta_update = "count"),
               (1 + 3) / (1 + p + p), tolerance = 1e-14)
})

test_that("update_mixture_prob('expected') matches (sum(p) + a - 1) / (p + a + b - 2)", {
  pips <- c(0.9, 0.05, 0.6, 0.02, 0.3)
  p    <- length(pips)
  expect_equal(update_mixture_prob(pips = pips, a = 1, b = 1,
                                   theta_update = "expected"),
               (sum(pips) + 1 - 1) / (p + 1 + 1 - 2),
               tolerance = 1e-14)
  expect_equal(update_mixture_prob(pips = pips, a = 2, b = p,
                                   theta_update = "expected"),
               (sum(pips) + 2 - 1) / (p + 2 + p - 2),
               tolerance = 1e-14)
})

test_that("update_mixture_prob('expected_clamp') floors at 0.001", {
  pips_tiny <- rep(1e-8, 100)
  expect_equal(
    update_mixture_prob(pips = pips_tiny, a = 1, b = 1,
                        theta_update = "expected_clamp"),
    0.001, tolerance = 1e-14
  )
})

test_that("update_mixture_prob('expected_clamp') equals 'expected' above the floor", {
  pips <- c(0.9, 0.05, 0.6, 0.02, 0.3)
  unclamped <- update_mixture_prob(pips = pips, a = 1, b = 1,
                                   theta_update = "expected")
  expect_gt(unclamped, 0.001)
  expect_equal(
    update_mixture_prob(pips = pips, a = 1, b = 1,
                        theta_update = "expected_clamp"),
    unclamped,
    tolerance = 1e-14
  )
})

test_that("update_mixture_prob errors when the required input is NULL", {
  expect_error(
    update_mixture_prob(pips = c(0.1, 0.2), theta_update = "count"),
    regexp = "requires 'betas'"
  )
  expect_error(
    update_mixture_prob(betas = c(1, 0), theta_update = "expected"),
    regexp = "requires 'pips'"
  )
  expect_error(
    update_mixture_prob(betas = c(1, 0), theta_update = "expected_clamp"),
    regexp = "requires 'pips'"
  )
})

test_that("expected_inclusion_probs with scalar theta matches Bayes formula", {
  betas <- c(0.5, 0.0, -0.3, 0.0, 0.1)
  s0 <- 0.04; s1 <- 0.5; theta <- 0.3

  pips_got <- expected_inclusion_probs(s0, s1, betas, theta = theta)

  d1 <- dlaplace(betas, mu = 0, b = s1)
  d0 <- dlaplace(betas, mu = 0, b = s0)
  pips_expected <- theta * d1 / (theta * d1 + (1 - theta) * d0)

  expect_equal(pips_got, pips_expected, tolerance = 1e-14)
})

test_that("expected_inclusion_probs accepts a length-p theta vector", {
  betas <- c(0.5, 0.0, -0.3, 0.0, 0.1)
  s0 <- 0.04; s1 <- 0.5
  theta_vec <- c(0.1, 0.2, 0.3, 0.4, 0.5)

  pips_got <- expected_inclusion_probs(s0, s1, betas, theta = theta_vec)

  d1 <- dlaplace(betas, mu = 0, b = s1)
  d0 <- dlaplace(betas, mu = 0, b = s0)
  pips_expected <- theta_vec * d1 /
                   (theta_vec * d1 + (1 - theta_vec) * d0)

  expect_equal(pips_got, pips_expected, tolerance = 1e-14)
})

test_that("expected_inclusion_probs errors when theta length is neither 1 nor p", {
  betas <- c(0.1, -0.2, 0.0)
  expect_error(
    expected_inclusion_probs(0.04, 0.5, betas, theta = c(0.1, 0.2)),
    regexp = "length 1 or length\\(betas\\)"
  )
  expect_error(
    expected_inclusion_probs(0.04, 0.5, betas, theta = c(0.1, 0.2, 0.3, 0.4)),
    regexp = "length 1 or length\\(betas\\)"
  )
})

test_that("expected_inclusion_probs errors when theta is outside (0, 1)", {
  betas <- c(0.1, -0.2, 0.0)
  expect_error(
    expected_inclusion_probs(0.04, 0.5, betas, theta = 0),
    regexp = "strictly between 0 and 1"
  )
  expect_error(
    expected_inclusion_probs(0.04, 0.5, betas, theta = 1),
    regexp = "strictly between 0 and 1"
  )
  expect_error(
    expected_inclusion_probs(0.04, 0.5, betas, theta = 1.5),
    regexp = "strictly between 0 and 1"
  )
})

test_that("count-rule fit: theta_path[1] = initial_sparsity and theta follows count rule", {
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
               update_mixture_prob(betas = fit$coefficients$Estimate,
                                   a = theta_a, b = theta_b,
                                   theta_update = "count"),
               tolerance = 1e-12)
})

test_that("expected-rule fit: theta_path[1] = initial_sparsity and theta follows expected rule", {
  dat <- make_ssl_fixture_data()
  x <- dat$x; y <- dat$y
  theta_a <- 1; theta_b <- ncol(x)
  initial_sparsity <- 0.05

  fit <- fit_ssl_psdh(x, y,
    ss = c(0.04, 0.5), initial_sparsity = initial_sparsity,
    theta_a = theta_a, theta_b = theta_b,
    theta_update = "expected",
    init_method = "LASSO_cv")

  expect_equal(fit$theta_update, "expected")
  expect_equal(fit$theta_path[1], initial_sparsity)
  expect_equal(fit$theta,
               update_mixture_prob(pips = fit$pips,
                                   a = theta_a, b = theta_b,
                                   theta_update = "expected"),
               tolerance = 1e-12)
})

test_that("expected_clamp fit: theta_path stays >= 0.001 and theta follows clamped rule", {
  dat <- make_ssl_fixture_data()
  x <- dat$x; y <- dat$y
  theta_a <- 1; theta_b <- ncol(x)
  initial_sparsity <- 0.05

  fit <- fit_ssl_psdh(x, y,
    ss = c(0.04, 0.5), initial_sparsity = initial_sparsity,
    theta_a = theta_a, theta_b = theta_b,
    theta_update = "expected_clamp",
    init_method = "LASSO_cv")

  expect_equal(fit$theta_update, "expected_clamp")
  expect_true(all(fit$theta_path[-1] >= 0.001))
  expect_equal(fit$theta,
               update_mixture_prob(pips = fit$pips,
                                   a = theta_a, b = theta_b,
                                   theta_update = "expected_clamp"),
               tolerance = 1e-12)
})

test_that("default initial_sparsity is 0.5", {
  expect_equal(formals(fit_ssl_psdh)$initial_sparsity, 0.5)
})
