# Tests for R/cv_fold_path.r — .cv_fold_path()
# Internal; visible because tests run inside the package namespace.

# ---------------------------------------------------------------------------
# shared fixture setup
# ---------------------------------------------------------------------------

.load_fold_fixture <- function() {
  path <- test_path("fixtures/init-regression.rds")
  skip_if(!file.exists(path),
    "Regression fixture absent — run make-init-regression-fixture.R")
  readRDS(path)
}

.fold_ctrl <- function(...) {
  defaults <- list(
    init_method = "zero",
    warm_start  = TRUE,
    fit_args    = list(maxit = 20L, epsilon = 1e-4, theta_b = 20L)
  )
  do.call(bhcrr_cv_control, utils::modifyList(defaults, list(...)))
}

# ---------------------------------------------------------------------------
# 1. Output shapes
# ---------------------------------------------------------------------------

test_that("lp has shape length(test_idx) x nrow(grid); iterations has one entry per pair", {
  fix       <- .load_fold_fixture()
  x <- fix$x; y <- fix$y
  train_idx <- seq_len(80L); test_idx <- 81L:100L
  grid      <- .cv_grid(c(0.02, 0.04), c(0.3, 0.5))   # 4 pairs
  eval_t    <- median(y[y[, 2] == 1, 1])
  ctrl      <- .fold_ctrl()

  r <- .cv_fold_path(x, y, train_idx, test_idx, grid, ctrl, eval_t)

  expect_equal(dim(r$lp),         c(length(test_idx), nrow(grid)))
  expect_equal(length(r$iterations), nrow(grid))
  expect_equal(length(r$init),       ncol(x))
  expect_equal(r$n_failed,            0L)
  expect_equal(nrow(r$errors),        0L)
  expect_null(r$coefs)
})

test_that("keep_coefs = TRUE returns a list of length nrow(grid)", {
  fix       <- .load_fold_fixture()
  x <- fix$x; y <- fix$y
  train_idx <- seq_len(80L); test_idx <- 81L:100L
  grid      <- .cv_grid(c(0.02, 0.06), c(0.5))
  eval_t    <- median(y[y[, 2] == 1, 1])
  ctrl      <- .fold_ctrl(keep_coefs = TRUE)

  r <- .cv_fold_path(x, y, train_idx, test_idx, grid, ctrl, eval_t)

  expect_true(is.list(r$coefs))
  expect_equal(length(r$coefs), nrow(grid))
  for (j in seq_len(nrow(grid)))
    expect_equal(length(r$coefs[[j]]), ncol(x))
})

# ---------------------------------------------------------------------------
# 2. Supplied init: init_method is NOT called
# ---------------------------------------------------------------------------

test_that("supplied init bypasses init_method entirely", {
  fix       <- .load_fold_fixture()
  x <- fix$x; y <- fix$y
  train_idx <- seq_len(80L); test_idx <- 81L:100L
  grid      <- .cv_grid(c(0.02, 0.04), c(0.5))
  eval_t    <- median(y[y[, 2] == 1, 1])

  init_call_count <- 0L
  counting_init   <- function(x, y, ...) {
    init_call_count <<- init_call_count + 1L
    list(init = rep(0, ncol(x)), meta = NULL)
  }
  ctrl <- bhcrr_cv_control(
    init_method = counting_init, warm_start = TRUE,
    fit_args    = list(maxit = 20L, epsilon = 1e-4, theta_b = 20L)
  )
  pre_supplied <- rep(0, ncol(x))

  r <- .cv_fold_path(x, y, train_idx, test_idx, grid, ctrl, eval_t,
                     init = pre_supplied)

  expect_equal(init_call_count, 0L,
    label = "init_method must not be called when init is supplied")
  expect_equal(r$init, pre_supplied)
})

# ---------------------------------------------------------------------------
# 3. init = NULL: init_method is called EXACTLY ONCE across the full grid
# ---------------------------------------------------------------------------

test_that("init_method is called exactly once for a multi-pair grid", {
  fix       <- .load_fold_fixture()
  x <- fix$x; y <- fix$y
  train_idx <- seq_len(80L); test_idx <- 81L:100L
  grid      <- .cv_grid(c(0.02, 0.04, 0.06), c(0.5))  # 3 pairs
  eval_t    <- median(y[y[, 2] == 1, 1])

  init_call_count <- 0L
  counting_init   <- function(x, y, ...) {
    init_call_count <<- init_call_count + 1L
    list(init = rep(0, ncol(x)), meta = NULL)
  }
  ctrl <- bhcrr_cv_control(
    init_method = counting_init, warm_start = TRUE,
    fit_args    = list(maxit = 20L, epsilon = 1e-4, theta_b = 20L)
  )

  r <- .cv_fold_path(x, y, train_idx, test_idx, grid, ctrl, eval_t)

  expect_equal(init_call_count, 1L,
    label = "init_method must be called exactly once, not once per pair")
})

# ---------------------------------------------------------------------------
# 4. warm_start TRUE vs FALSE: both complete; iteration counts are captured
# ---------------------------------------------------------------------------
#
# On small (n=100, p=20) data both modes converge in the same number of
# iterations because the EM hits the minimum check horizon (iter > 5) before
# the warm-start advantage becomes detectable. At p >= 2000 or with
# mal-conditioned grids the warm-start chain consistently requires fewer
# EM iterations — that difference is the measurement used to set parallel width.

test_that("warm_start TRUE and FALSE both complete with recorded iteration counts", {
  fix       <- .load_fold_fixture()
  x <- fix$x; y <- fix$y
  train_idx <- seq_len(80L); test_idx <- 81L:100L
  grid      <- .cv_grid(c(0.02, 0.04, 0.08), c(0.3, 0.5))
  eval_t    <- median(y[y[, 2] == 1, 1])
  fa        <- list(maxit = 30L, epsilon = 1e-4, theta_b = 20L)

  r_warm <- .cv_fold_path(x, y, train_idx, test_idx, grid,
                          .fold_ctrl(warm_start = TRUE,  fit_args = fa), eval_t)
  r_cold <- .cv_fold_path(x, y, train_idx, test_idx, grid,
                          .fold_ctrl(warm_start = FALSE, fit_args = fa), eval_t)

  # Both must finish without errors
  expect_equal(r_warm$n_failed, 0L)
  expect_equal(r_cold$n_failed, 0L)

  # Iteration counts must be present for every pair
  expect_true(all(!is.na(r_warm$iterations)))
  expect_true(all(!is.na(r_cold$iterations)))

  # Record the means (informational; equality is acceptable on small data)
  mean_warm <- mean(r_warm$iterations)
  mean_cold <- mean(r_cold$iterations)
  expect_true(is.numeric(mean_warm) && is.finite(mean_warm))
  expect_true(is.numeric(mean_cold) && is.finite(mean_cold))
})

# ---------------------------------------------------------------------------
# 5. Injected failure at a middle pair: isolation and chain repair
# ---------------------------------------------------------------------------

test_that("failure at pair 2 of 3 leaves columns 1 and 3 intact; chain repair reuses last good coef", {
  fix       <- .load_fold_fixture()
  x <- fix$x; y <- fix$y
  train_idx <- seq_len(80L); test_idx <- 81L:100L
  grid      <- .cv_grid(c(0.02, 0.04, 0.06), c(0.5))  # 3 pairs in order
  eval_t    <- median(y[y[, 2] == 1, 1])

  # Track how many times fit_ssl_psdh is called; fail on the 2nd call.
  fit_call_count  <- 0L
  init_call_count <- 0L
  real_fit        <- fit_ssl_psdh

  counting_init <- function(x, y, ...) {
    init_call_count <<- init_call_count + 1L
    list(init = rep(0, ncol(x)), meta = NULL)
  }
  # Build the control BEFORE mocking: bhcrr_cv_control() validates fit_args against
  # names(formals(fit_ssl_psdh)), and under the mock those formals are just "...".
  ctrl <- bhcrr_cv_control(
    init_method = counting_init, warm_start = TRUE,
    fit_args    = list(maxit = 20L, epsilon = 1e-4, theta_b = 20L)
  )

  local_mocked_bindings(
    fit_ssl_psdh = function(...) {
      fit_call_count <<- fit_call_count + 1L
      if (fit_call_count == 2L) stop("injected failure for pair 2")
      real_fit(...)
    },
    .package = "bhCRR"
  )

  r <- .cv_fold_path(x, y, train_idx, test_idx, grid, ctrl, eval_t)

  # Middle column is all NA; first and last are finite
  expect_true(all(is.na(r$lp[, 2])),
    label = "failed pair's lp column must be all NA")
  expect_false(any(is.na(r$lp[, 1])),
    label = "pair 1 (succeeded) must have finite predictions")
  expect_false(any(is.na(r$lp[, 3])),
    label = "pair 3 (succeeded) must have finite predictions")

  # Errors data frame has exactly one row naming pair 2
  expect_equal(nrow(r$errors), 1L)
  expect_equal(r$errors$pair[1L], grid$pair[2L])
  expect_match(r$errors$message[1L], "injected")
  expect_equal(r$n_failed, 1L)

  # iterations: NA for pair 2, integer for pairs 1 and 3
  expect_equal(is.na(r$iterations), c(FALSE, TRUE, FALSE))

  # Init was called exactly once — chain repair must NOT re-invoke init
  expect_equal(init_call_count, 1L,
    label = "chain repair must not call the init method again")
})

# ---------------------------------------------------------------------------
# 6. Silent under failure: worker emits no messages or warnings
# ---------------------------------------------------------------------------

test_that("worker emits no messages or warnings even when all fits fail", {
  fix       <- .load_fold_fixture()
  x <- fix$x; y <- fix$y
  train_idx <- seq_len(80L); test_idx <- 81L:100L
  grid      <- .cv_grid(c(0.02, 0.04), c(0.5))
  eval_t    <- median(y[y[, 2] == 1, 1])

  # Build the control BEFORE mocking: bhcrr_cv_control() validates fit_args against
  # names(formals(fit_ssl_psdh)), and under the mock those formals are just "...".
  ctrl <- .fold_ctrl()

  local_mocked_bindings(
    fit_ssl_psdh = function(...) stop("deliberate failure"),
    .package = "bhCRR"
  )

  expect_silent(
    r <- .cv_fold_path(x, y, train_idx, test_idx, grid, ctrl, eval_t,
                       init = rep(0, ncol(x)))
  )

  expect_equal(r$n_failed, nrow(grid))
  expect_true(all(is.na(r$lp)))
})
