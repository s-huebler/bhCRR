# Tests for the deferred Breslow jump computation introduced in fit_ssl_psdh.
# getBreslowJumps = FALSE must not alter update_betas() coefficients; the
# final refit must populate breslowJump / uftime correctly.

.load_breslow_fixture <- function() {
  path <- test_path("fixtures/init-regression.rds")
  skip_if(
    !file.exists(path),
    "Regression fixture absent — run tests/testthat/fixtures/make-init-regression-fixture.R"
  )
  readRDS(path)
}

# ---------------------------------------------------------------------------
# 1. Correctness assumption: getBreslowJumps does not affect $coef
# ---------------------------------------------------------------------------

test_that("update_betas $coef is identical whether getBreslowJumps is TRUE or FALSE", {
  fix <- .load_breslow_fixture()
  pw  <- rep(1, ncol(fix$x))
  lam <- 1 / nrow(fix$x)

  m_true  <- update_betas(pw, fix$y[, 1], fix$y[, 2], fix$x,
                          lambda = lam, getBreslowJumps = TRUE,  max.iter = 1000)
  m_false <- update_betas(pw, fix$y[, 1], fix$y[, 2], fix$x,
                          lambda = lam, getBreslowJumps = FALSE, max.iter = 1000)

  expect_equal(m_true$coef, m_false$coef, tolerance = 1e-12,
    info = paste(
      "getBreslowJumps must not alter the coefficient estimates.",
      "If this fails, fastCrrp is not deterministic given identical arguments,",
      "which invalidates the deferred-refit approach."
    )
  )
})

test_that("update_betas with getBreslowJumps=TRUE returns a matrix breslowJump", {
  fix <- .load_breslow_fixture()
  pw  <- rep(1, ncol(fix$x))
  lam <- 1 / nrow(fix$x)
  m   <- update_betas(pw, fix$y[, 1], fix$y[, 2], fix$x,
                      lambda = lam, getBreslowJumps = TRUE, max.iter = 1000)
  expect_true(is.matrix(m$breslowJump) || is.data.frame(m$breslowJump))
  expect_true(!is.null(m$uftime))
})

# ---------------------------------------------------------------------------
# 2. fit$final_model_object has breslowJump and uftime populated
# ---------------------------------------------------------------------------

test_that("fit_ssl_psdh final_model_object has a non-degenerate breslowJump and uftime", {
  fix <- .load_breslow_fixture()
  fit <- fit_ssl_psdh(fix$x, fix$y,
    ss = c(0.04, 0.5), initial_sparsity = 0.05,
    theta_a = 1, theta_b = ncol(fix$x),
    maxit = 50, epsilon = 1e-4,
    init_method = "LASSO_cv")

  bj <- fit$final_model_object$breslowJump
  expect_true(is.matrix(bj) || is.data.frame(bj),
    label = "breslowJump should be a matrix/data.frame, not scalar FALSE or NULL")
  expect_gt(nrow(bj), 0L)
  expect_true(!is.null(fit$final_model_object$uftime))
  expect_gt(length(fit$final_model_object$uftime), 0L)
})

# ---------------------------------------------------------------------------
# 3. predict_from_ssl_psdh() runs and returns finite risk scores
# ---------------------------------------------------------------------------

test_that("predict_from_ssl_psdh returns finite risk scores on fixture data", {
  fix <- .load_breslow_fixture()
  fit <- fit_ssl_psdh(fix$x, fix$y,
    ss = c(0.04, 0.5), initial_sparsity = 0.05,
    theta_a = 1, theta_b = ncol(fix$x),
    maxit = 50, epsilon = 1e-4,
    init_method = "LASSO_cv")

  pred_time <- median(fix$y[, 1])
  risk <- predict_from_ssl_psdh(fit, newx = fix$x, prediction_time = pred_time)

  expect_equal(length(risk), nrow(fix$x))
  expect_true(all(is.finite(risk)))
  expect_true(all(risk >= 0 & risk <= 1))
})

# ---------------------------------------------------------------------------
# 4. $conv is the EM convergence flag, not the final refit's converged field
# ---------------------------------------------------------------------------

test_that("fit$conv reflects EM convergence, not the final refit's converged slot", {
  fix <- .load_breslow_fixture()
  fit <- fit_ssl_psdh(fix$x, fix$y,
    ss = c(0.04, 0.5), initial_sparsity = 0.05,
    theta_a = 1, theta_b = ncol(fix$x),
    maxit = 50, epsilon = 1e-4,
    init_method = "LASSO_cv")

  expect_true(is.logical(fit$conv))
  # Fixture converges, so this should be TRUE; the important thing is that
  # it is a scalar logical, not the final_mod$converged integer.
  expect_equal(length(fit$conv), 1L)
})
