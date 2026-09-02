# Tests for R/init_methods.r
# Internal dot-prefixed functions are visible here because testthat runs inside
# the package namespace under devtools::load_all() / test_check().

# ---------------------------------------------------------------------------
# Shared fixture loader (used by tests that need real data)
# ---------------------------------------------------------------------------

.load_init_fixture <- function() {
  path <- test_path("fixtures/init-regression.rds")
  skip_if(
    !file.exists(path),
    "Regression fixture absent — run tests/testthat/fixtures/make-init-regression-fixture.R"
  )
  readRDS(path)
}

# ---------------------------------------------------------------------------
# 1. .LASSO_cv reproduces the stored init_ref
# ---------------------------------------------------------------------------

test_that(".LASSO_cv$init matches stored init_ref (proves which.max is behaviour-preserving)", {
  fix <- .load_init_fixture()
  got <- .LASSO_cv(fix$x, fix$y)
  expect_equal(
    got$init, fix$init_ref,
    tolerance = 1e-12,
    label     = "which.max index selection",
    expected.label = "stored init_ref",
    info = paste(
      "which.max(cv_c_index) must select the same column as the old",
      "lambda == lambda_min floating-point equality.",
      "If this fails, the index-based refactor is NOT behaviour-preserving."
    )
  )
})

test_that(".LASSO_cv returns a list with numeric init of length ncol(x) and non-null meta", {
  fix <- .load_init_fixture()
  res <- .LASSO_cv(fix$x, fix$y)
  expect_true(is.list(res))
  expect_true(is.numeric(res$init))
  expect_equal(length(res$init), ncol(fix$x))
  expect_true(is.list(res$meta))
  expect_true(!is.null(res$meta$lambda_min))
  expect_true(res$meta$elapsed_sec >= 0)
})

# ---------------------------------------------------------------------------
# 2. .LASSO_bic on fixture data
# ---------------------------------------------------------------------------

test_that(".LASSO_bic returns a finite length-p vector with df > 0", {
  fix <- .load_init_fixture()
  res <- .LASSO_bic(fix$x, fix$y)
  expect_true(is.list(res))
  expect_true(is.numeric(res$init))
  expect_equal(length(res$init), ncol(fix$x))
  expect_true(all(is.finite(res$init)))
  expect_gt(sum(res$init != 0), 0)
  expect_gt(res$meta$df, 0L)
})

# ---------------------------------------------------------------------------
# 3. .LASSO_bic errors informatively when every path solution is null
# ---------------------------------------------------------------------------

test_that(".LASSO_bic errors when all path solutions are null (absurdly large lambda)", {
  set.seed(1)
  x_tiny <- matrix(rnorm(50 * 5), 50, 5)
  colnames(x_tiny) <- paste0("V", 1:5)
  y_tiny <- cbind(
    time   = rexp(50) + 0.01,
    status = sample(c(0L, 1L, 2L), 50, replace = TRUE, prob = c(0.3, 0.5, 0.2))
  )
  expect_error(
    .LASSO_bic(x_tiny, y_tiny, lambda_path = c(1e4, 1e5, 1e6)),
    regexp = "non-null|No finite"
  )
})

# ---------------------------------------------------------------------------
# 4. .zero_init
# ---------------------------------------------------------------------------

test_that(".zero_init returns rep(0, p) with method meta", {
  fix <- .load_init_fixture()
  res <- .zero_init(fix$x, fix$y)
  expect_equal(res$init, rep(0, ncol(fix$x)))
  expect_equal(res$meta$method, "zero")
})

test_that(".zero_init output passes .validate_init", {
  fix <- .load_init_fixture()
  p   <- ncol(fix$x)
  res <- .zero_init(fix$x, fix$y)
  val <- .validate_init(res, p, "zero")
  expect_equal(val$init, rep(0, p))
})

# ---------------------------------------------------------------------------
# 5. .resolve_init_method
# ---------------------------------------------------------------------------

test_that(".resolve_init_method resolves built-in names case-insensitively", {
  r1 <- .resolve_init_method("LASSO_cv")
  expect_equal(r1$label, "LASSO_cv")
  expect_identical(r1$fn, .LASSO_cv)

  r2 <- .resolve_init_method("lasso_cv")
  expect_equal(r2$label, "LASSO_cv")
  expect_identical(r2$fn, .LASSO_cv)

  r3 <- .resolve_init_method("LASSO_BIC")
  expect_equal(r3$label, "LASSO_bic")
  expect_identical(r3$fn, .LASSO_bic)
})

test_that(".resolve_init_method passes a function through with label 'custom'", {
  my_fn <- function(x, y) list(init = rep(0, ncol(x)), meta = NULL)
  res   <- .resolve_init_method(my_fn)
  expect_equal(res$label, "custom")
  expect_identical(res$fn, my_fn)
})

test_that(".resolve_init_method resolves a name from an explicitly passed envir", {
  resolve_from_local <- function() {
    local_init_fn <- function(x, y) list(init = rep(0.1, ncol(x)), meta = NULL)
    .resolve_init_method("local_init_fn", envir = environment())
  }
  res <- resolve_from_local()
  expect_equal(res$label, "custom")
  expect_true(is.function(res$fn))
})

test_that(".resolve_init_method errors on unknown string and lists built-in names", {
  expect_error(
    .resolve_init_method("not_a_method"),
    regexp = "LASSO_cv"
  )
})

test_that(".resolve_init_method errors on non-character/non-function input", {
  expect_error(
    .resolve_init_method(42L),
    regexp = "integer"
  )
})

# ---------------------------------------------------------------------------
# 6. .validate_init
# ---------------------------------------------------------------------------

test_that(".validate_init accepts a bare numeric vector", {
  val <- .validate_init(c(0.1, 0.2, 0.3), p = 3L, method_label = "test_method")
  expect_equal(val$init, c(0.1, 0.2, 0.3))
  expect_null(val$meta)
})

test_that(".validate_init accepts list(init = , meta = )", {
  val <- .validate_init(
    list(init = c(1, 2), meta = list(foo = "bar")),
    p = 2L, method_label = "test_method"
  )
  expect_equal(val$init, c(1, 2))
  expect_equal(val$meta$foo, "bar")
})

test_that(".validate_init errors on wrong length and names the method label", {
  expect_error(
    .validate_init(c(1, 2, 3), p = 5L, method_label = "my_custom_fn"),
    regexp = "my_custom_fn"
  )
})

test_that(".validate_init errors on non-numeric and names the method label", {
  expect_error(
    .validate_init(list(init = c("a", "b"), meta = NULL), p = 2L, method_label = "bad_fn"),
    regexp = "bad_fn"
  )
})

test_that(".validate_init errors on NA", {
  expect_error(
    .validate_init(c(1, NA, 3), p = 3L, method_label = "na_fn"),
    regexp = "na_fn"
  )
})

test_that(".validate_init errors on NaN", {
  expect_error(
    .validate_init(c(1, NaN, 3), p = 3L, method_label = "nan_fn"),
    regexp = "nan_fn"
  )
})

test_that(".validate_init errors on Inf", {
  expect_error(
    .validate_init(c(1, Inf, 3), p = 3L, method_label = "inf_fn"),
    regexp = "inf_fn"
  )
})

# ---------------------------------------------------------------------------
# 7. Custom initializer end-to-end: resolve -> call -> validate (bare-vector return)
# ---------------------------------------------------------------------------

test_that("custom bare-vector initializer round-trips through resolve/call/validate", {
  set.seed(1)
  x_small <- matrix(rnorm(30 * 4), 30, 4)
  colnames(x_small) <- paste0("V", 1:4)
  y_small <- cbind(
    time   = rexp(30) + 0.01,
    status = sample(c(0L, 1L, 2L), 30, replace = TRUE, prob = c(0.3, 0.5, 0.2))
  )

  bare_fn <- function(x, y, ...) rep(0.1, ncol(x))

  resolved <- .resolve_init_method(bare_fn)
  expect_equal(resolved$label, "custom")

  raw <- resolved$fn(x_small, y_small)
  val <- .validate_init(raw, p = ncol(x_small), method_label = "custom")
  expect_equal(val$init, rep(0.1, ncol(x_small)))
  expect_null(val$meta)
})
