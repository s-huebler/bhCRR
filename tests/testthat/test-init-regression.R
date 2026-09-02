# This test pins the current behaviour of fit_ssl_psdh()'s built-in CV-LASSO
# initialization as a regression fixture, so the upcoming init refactor can be
# proven behaviour-preserving.
# The fit_ssl_psdh() call now passes init_method = "LASSO_cv" (the refactored
# form) instead of the old inline cv_fastCrrp_cpp() block.

test_that("fit_ssl_psdh CV-LASSO init reproduces stored init_ref", {
  fixture_path <- test_path("fixtures/init-regression.rds")
  skip_if(
    !file.exists(fixture_path),
    "Regression fixture not found — run tests/testthat/fixtures/make-init-regression-fixture.R to generate it"
  )

  fix <- readRDS(fixture_path)
  x   <- fix$x
  y   <- fix$y

  init_lam_path <- 10^seq(log10(0.1), log10(0.001), length = 10)
  cv_out <- cv_fastCrrp_cpp(x, y[, 1], y[, 2], k = 5, penalty = "LASSO",
                            lambda_path = init_lam_path,
                            tuning = "wolbers", eval_quantile = 0.5)
  init_got <- cv_out$full_model$coef[, cv_out$lambda == cv_out$lambda_min]

  expect_equal(init_got, fix$init_ref, tolerance = 1e-10)
})

test_that("fit_ssl_psdh with init_method = 'LASSO_cv' reproduces stored fit_ref", {
  fixture_path <- test_path("fixtures/init-regression.rds")
  skip_if(
    !file.exists(fixture_path),
    "Regression fixture not found — run tests/testthat/fixtures/make-init-regression-fixture.R to generate it"
  )

  fix <- readRDS(fixture_path)
  x   <- fix$x
  y   <- fix$y

  fit <- fit_ssl_psdh(x, y,
    ss = c(0.04, 0.5), initial_sparsity = 0.05,
    theta_a = 1, theta_b = ncol(x),
    maxit = 50, epsilon = 1e-04,
    init_method = "LASSO_cv")

  ref <- fix$fit_ref
  expect_equal(fit$coefficients$Estimate, ref$Estimate,       tolerance = 1e-10)
  expect_equal(fit$pips,                  ref$pips,           tolerance = 1e-10)
  expect_equal(fit$penalty.factor,        ref$penalty.factor, tolerance = 1e-10)
  expect_equal(fit$iterations,            ref$iterations)
  expect_equal(fit$conv,                  ref$conv)
  expect_equal(fit$ss,                    ref$ss,             tolerance = 1e-10)
})

test_that("fit_ssl_psdh with init = <numeric> works and sets $init_method to 'supplied'", {
  fixture_path <- test_path("fixtures/init-regression.rds")
  skip_if(
    !file.exists(fixture_path),
    "Regression fixture not found — run tests/testthat/fixtures/make-init-regression-fixture.R to generate it"
  )

  fix <- readRDS(fixture_path)
  x   <- fix$x
  y   <- fix$y

  fit <- fit_ssl_psdh(x, y,
    ss = c(0.04, 0.5), initial_sparsity = 0.05,
    theta_a = 1, theta_b = ncol(x),
    maxit = 50, epsilon = 1e-04,
    init = fix$init_ref)

  expect_equal(fit$init_method, "supplied")
  expect_equal(fit$init, fix$init_ref, tolerance = 1e-14)
  expect_null(fit$init_meta)
})

test_that("fit_ssl_psdh errors when neither init nor init_method is supplied", {
  fixture_path <- test_path("fixtures/init-regression.rds")
  skip_if(
    !file.exists(fixture_path),
    "Regression fixture not found — run tests/testthat/fixtures/make-init-regression-fixture.R to generate it"
  )

  fix <- readRDS(fixture_path)
  expect_error(
    fit_ssl_psdh(fix$x, fix$y, ss = c(0.04, 0.5), initial_sparsity = 0.05,
                 theta_a = 1, theta_b = ncol(fix$x), maxit = 2),
    regexp = "now requires an initialization"
  )
})

test_that("fit_ssl_psdh errors when both init and init_method are supplied", {
  fixture_path <- test_path("fixtures/init-regression.rds")
  skip_if(
    !file.exists(fixture_path),
    "Regression fixture not found — run tests/testthat/fixtures/make-init-regression-fixture.R to generate it"
  )

  fix <- readRDS(fixture_path)
  expect_error(
    fit_ssl_psdh(fix$x, fix$y, ss = c(0.04, 0.5), initial_sparsity = 0.05,
                 theta_a = 1, theta_b = ncol(fix$x), maxit = 2,
                 init = rep(0, ncol(fix$x)), init_method = "zero"),
    regexp = "not both"
  )
})

test_that("fit_ssl_psdh errors with defunct message when init_lam_path is supplied", {
  fixture_path <- test_path("fixtures/init-regression.rds")
  skip_if(
    !file.exists(fixture_path),
    "Regression fixture not found — run tests/testthat/fixtures/make-init-regression-fixture.R to generate it"
  )

  fix <- readRDS(fixture_path)
  expect_error(
    fit_ssl_psdh(fix$x, fix$y, ss = c(0.04, 0.5), initial_sparsity = 0.05,
                 theta_a = 1, theta_b = ncol(fix$x), maxit = 2,
                 init_method = "zero",
                 init_lam_path = c(0.1, 0.01)),
    regexp = "defunct"
  )
})
