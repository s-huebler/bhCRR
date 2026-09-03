# Tests for R/cv.r — bhcrr_cv()

# ---------------------------------------------------------------------------
# shared helpers
# ---------------------------------------------------------------------------

.load_cv_fixture <- function() {
  path <- test_path("fixtures/init-regression.rds")
  skip_if(!file.exists(path),
    "Regression fixture absent — run make-init-regression-fixture.R")
  readRDS(path)
}

.cv_ctrl <- function(...) {
  defaults <- list(
    init_method = "zero",
    nfolds      = 3L,
    ncv         = 2L,
    seed        = 42L,
    fit_args    = list(maxit = 20L, epsilon = 1e-4, theta_b = 20L)
  )
  do.call(bhcrr_cv_control, utils::modifyList(defaults, list(...)))
}

s0_small <- c(0.02, 0.04, 0.06)
s1_small <- c(0.3,  0.5)

# ---------------------------------------------------------------------------
# 1. Basic contract: object class, tuning shape, column names, print, validate
# ---------------------------------------------------------------------------

test_that("returns bhcrr_cv; tuning has correct shape and column names; print works", {
  fix <- .load_cv_fixture()
  x <- fix$x; y <- fix$y
  ctrl <- .cv_ctrl()

  r <- bhcrr_cv(x, y, s0_small, s1_small, ctrl)

  expect_s3_class(r, "bhcrr_cv")

  expected_cols <- c("s0", "s1", "score_mean", "score_sd",
                     "score_fold_sd", "n_failed_folds", "n_not_converged")
  expect_named(r$tuning, expected_cols)
  expect_equal(nrow(r$tuning), nrow(.cv_grid(s0_small, s1_small)))

  # print() must not error and must mention "bhcrr_cv"
  expect_output(print(r), "bhcrr_cv")
})

test_that("bhcrr_tune_validate() accepts tuning without modification", {
  fix <- .load_cv_fixture()
  x <- fix$x; y <- fix$y
  ctrl <- .cv_ctrl()

  r <- bhcrr_cv(x, y, s0_small, s1_small, ctrl)

  val <- bhcrr_tune_validate(
    r$tuning, n = nrow(x), theta = 20,
    s0_seq = s0_small, s1_seq = s1_small
  )
  expect_true(is.list(val))
  expect_true("checks" %in% names(val))
})

# ---------------------------------------------------------------------------
# 2. One-pair grid
# ---------------------------------------------------------------------------

test_that("a grid with one valid pair returns tuning with one row", {
  fix <- .load_cv_fixture()
  x <- fix$x; y <- fix$y
  ctrl <- .cv_ctrl()

  r <- bhcrr_cv(x, y, s0_seq = 0.02, s1_seq = 0.5, ctrl)

  expect_equal(nrow(r$tuning), 1L)
  expect_equal(dim(r$pooled),       c(1L, ctrl$ncv))
  expect_equal(dim(r$fold_scores),  c(1L, ctrl$nfolds, ctrl$ncv))
})

# ---------------------------------------------------------------------------
# 3. eval_time source tracking
# ---------------------------------------------------------------------------

test_that("control$eval_time is used verbatim; $eval_time$source is 'supplied'", {
  fix <- .load_cv_fixture()
  x <- fix$x; y <- fix$y
  ctrl <- .cv_ctrl(eval_time = 5.0)

  r <- bhcrr_cv(x, y, s0_small, s1_small, ctrl)

  expect_equal(r$eval_time$value,  5.0)
  expect_equal(r$eval_time$source, "supplied")
})

test_that("NULL eval_time derives tau from quantile; source is 'quantile'", {
  fix <- .load_cv_fixture()
  x <- fix$x; y <- fix$y
  ctrl <- .cv_ctrl(eval_time = NULL, eval_quantile = 0.5)

  r <- bhcrr_cv(x, y, s0_small, s1_small, ctrl)

  expected_tau <- median(y[y[, 2] == 1, 1])
  expect_equal(r$eval_time$value,  expected_tau)
  expect_equal(r$eval_time$source, "quantile")
})

# ---------------------------------------------------------------------------
# 4. fold_inits round-trip: identical tuning; init never called on second run
# ---------------------------------------------------------------------------

test_that("supplying fold_inits from a first run reproduces tuning exactly without calling init", {
  fix <- .load_cv_fixture()
  x <- fix$x; y <- fix$y

  call_count_1 <- 0L
  init_fn_1 <- function(x, y, ...) {
    call_count_1 <<- call_count_1 + 1L
    list(init = rep(0, ncol(x)), meta = NULL)
  }
  ctrl1 <- .cv_ctrl(init_method = init_fn_1)

  # Pre-generate folds so both runs use identical fold assignments.
  fold_obj <- bhcrr_make_folds(y, ctrl1)

  r1 <- bhcrr_cv(x, y, s0_small, s1_small, ctrl1, folds = fold_obj)

  expect_equal(call_count_1, fold_obj$nfolds * fold_obj$ncv,
    label = "init called once per fold on first run")

  # Second run: different counter, same fold structure, supply init vectors.
  call_count_2 <- 0L
  init_fn_2 <- function(x, y, ...) {
    call_count_2 <<- call_count_2 + 1L
    list(init = rep(0, ncol(x)), meta = NULL)
  }
  ctrl2 <- .cv_ctrl(init_method = init_fn_2)

  r2 <- bhcrr_cv(x, y, s0_small, s1_small, ctrl2,
                 folds = fold_obj, fold_inits = r1$fold_inits)

  expect_equal(call_count_2, 0L,
    label = "init method must not be called when fold_inits is supplied")
  expect_equal(r1$tuning, r2$tuning,
    label = "tuning must be identical when fold_inits is replayed")
})

# ---------------------------------------------------------------------------
# 5. Partial init failure: call completes, one warning, errors logged
# ---------------------------------------------------------------------------

test_that("init failure in one fold: one warning, $errors records that fold, other folds unaffected", {
  fix <- .load_cv_fixture()
  x <- fix$x; y <- fix$y

  fail_count <- 0L
  partial_init <- function(x, y, ...) {
    fail_count <<- fail_count + 1L
    if (fail_count == 2L) stop("deliberate init failure for fold 2")
    list(init = rep(0, ncol(x)), meta = NULL)
  }
  # Build control and folds BEFORE any mock so validation uses real formals.
  ctrl     <- .cv_ctrl(init_method = partial_init)
  fold_obj <- bhcrr_make_folds(y, ctrl)

  expect_warning(
    r <- bhcrr_cv(x, y, s0_small, s1_small, ctrl, folds = fold_obj),
    regexp = "fit failure"
  )

  # Exactly one fold-level error row (pair = NA).
  expect_equal(nrow(r$errors), 1L)
  expect_true(is.na(r$errors$pair[1L]))
  expect_equal(r$errors$rep[1L],  1L)
  expect_equal(r$errors$fold[1L], 2L)
  expect_match(r$errors$message[1L], "deliberate")

  # Rep 2 (no failure) must have all-finite fold scores.
  expect_true(all(is.finite(r$fold_scores[, , 2])),
    label = "rep 2 fold scores must be finite (unaffected by rep-1 failure)")

  # $best must be non-NULL because rep 2's scores are finite.
  expect_false(is.null(r$best))
})

# ---------------------------------------------------------------------------
# 6. Total init failure: stop() with first message quoted
# ---------------------------------------------------------------------------

test_that("init failure in every fold triggers stop() quoting the first message", {
  fix <- .load_cv_fixture()
  x <- fix$x; y <- fix$y

  ctrl <- .cv_ctrl(init_method = function(x, y, ...) stop("always fails"))

  expect_error(
    bhcrr_cv(x, y, s0_small, s1_small, ctrl),
    regexp = "always fails"
  )
})

# ---------------------------------------------------------------------------
# 8. Non-convergence tracking: conv=FALSE captured; errors stays empty
# ---------------------------------------------------------------------------
#
# fit_ssl_psdh's convergence check fires only when iter > 5.  With maxit=5
# that condition is never satisfied, so every fit returns conv=FALSE
# regardless of data, init, or hyperparameters.  The fits still succeed
# (no error), so $errors must stay empty.
#
# bhcrr_cv() must emit exactly one summary non-convergence warning — distinct
# from the per-fit "fit_ssl_psdh did not converge" warnings from the fitter.

test_that("maxit=5 bypasses convergence check; all conv=FALSE; $errors empty; one summary warning", {
  fix <- .load_cv_fixture()
  x <- fix$x; y <- fix$y

  # maxit=5: convergence check (iter > 5) never fires -> conv=FALSE for all fits.
  ctrl <- .cv_ctrl(fit_args = list(maxit = 5L, epsilon = 1e-4, theta_b = 20L))

  # Count only bhcrr_cv's summary warning (pattern "fit(s) did not converge");
  # muffle all warnings to keep test output clean.
  summary_nc_count <- 0L
  r <- withCallingHandlers(
    bhcrr_cv(x, y, s0_small, s1_small, ctrl),
    warning = function(w) {
      if (grepl("fit\\(s\\) did not converge", conditionMessage(w)))
        summary_nc_count <<- summary_nc_count + 1L
      invokeRestart("muffleWarning")
    }
  )

  expect_equal(summary_nc_count, 1L,
    label = "bhcrr_cv must emit exactly one non-convergence summary warning")

  # All convergence flags must be FALSE.
  expect_true(all(!r$convergence),
    label = "every (pair, fold, rep) entry in $convergence must be FALSE")

  # n_not_converged must equal nfolds * ncv for every pair.
  n_cells <- r$folds$nfolds * r$folds$ncv
  expect_true(all(r$tuning$n_not_converged == n_cells),
    label = "n_not_converged must equal nfolds * ncv for every pair")

  # Fits completed (no errors) — non-convergence is not failure.
  expect_equal(nrow(r$errors), 0L)
})

# ---------------------------------------------------------------------------
# 10. Dimensions of $pooled, $fold_scores, $iterations, $convergence
# ---------------------------------------------------------------------------

test_that("$pooled is n_pairs x ncv and $fold_scores is n_pairs x nfolds x ncv", {
  fix <- .load_cv_fixture()
  x <- fix$x; y <- fix$y
  ctrl     <- .cv_ctrl()
  fold_obj <- bhcrr_make_folds(y, ctrl)
  n_pairs  <- nrow(.cv_grid(s0_small, s1_small))

  r <- bhcrr_cv(x, y, s0_small, s1_small, ctrl, folds = fold_obj)

  expect_equal(dim(r$pooled),      c(n_pairs, fold_obj$ncv))
  expect_equal(dim(r$fold_scores), c(n_pairs, fold_obj$nfolds, fold_obj$ncv))
  expect_equal(dim(r$iterations),  c(n_pairs, fold_obj$nfolds, fold_obj$ncv))
  expect_equal(dim(r$convergence), c(n_pairs, fold_obj$nfolds, fold_obj$ncv))
})

# ---------------------------------------------------------------------------
# Parallel correctness
# ---------------------------------------------------------------------------

skip_parallel <- function() {
  skip_on_os("windows")
  skip_if(parallel::detectCores() < 2L, "fewer than 2 cores")
}

test_that("parallel warm_start=TRUE gives identical results to sequential", {
  skip_parallel()
  fix <- .load_cv_fixture()
  x <- fix$x; y <- fix$y
  ctrl_seq <- .cv_ctrl(parallel = FALSE, workers = NULL, seed = 1L)
  ctrl_par <- .cv_ctrl(parallel = TRUE,  workers = 2L,   seed = 1L)
  fold_obj <- bhcrr_make_folds(y, ctrl_seq)

  r_seq <- suppressWarnings(bhcrr_cv(x, y, s0_small, s1_small, ctrl_seq, folds = fold_obj))
  r_par <- suppressWarnings(bhcrr_cv(x, y, s0_small, s1_small, ctrl_par, folds = fold_obj))

  expect_equal(r_seq$tuning,      r_par$tuning)
  expect_equal(r_seq$pooled,      r_par$pooled)
  expect_equal(r_seq$fold_scores, r_par$fold_scores)
  expect_equal(r_seq$iterations,  r_par$iterations)
  expect_equal(r_seq$convergence, r_par$convergence)
})

test_that("parallel warm_start=FALSE gives identical results to sequential", {
  skip_parallel()
  fix <- .load_cv_fixture()
  x <- fix$x; y <- fix$y
  ctrl_seq <- .cv_ctrl(parallel = FALSE, workers = NULL, seed = 1L, warm_start = FALSE)
  ctrl_par <- .cv_ctrl(parallel = TRUE,  workers = 2L,   seed = 1L, warm_start = FALSE)
  fold_obj <- bhcrr_make_folds(y, ctrl_seq)

  r_seq <- suppressWarnings(bhcrr_cv(x, y, s0_small, s1_small, ctrl_seq, folds = fold_obj))
  r_par <- suppressWarnings(bhcrr_cv(x, y, s0_small, s1_small, ctrl_par, folds = fold_obj))

  expect_equal(r_seq$tuning,      r_par$tuning)
  expect_equal(r_seq$pooled,      r_par$pooled)
  expect_equal(r_seq$fold_scores, r_par$fold_scores)
  expect_equal(r_seq$iterations,  r_par$iterations)
  expect_equal(r_seq$convergence, r_par$convergence)
})

test_that("sequential warm_start=FALSE matches sequential warm_start=FALSE wide dispatch invariant", {
  # The wide path (warm_start=FALSE) must be identical to sequential warm_start=FALSE.
  # This is step 3's correctness invariant.
  skip_parallel()
  fix <- .load_cv_fixture()
  x <- fix$x; y <- fix$y
  ctrl_seq  <- .cv_ctrl(parallel = FALSE, seed = 1L, warm_start = FALSE)
  ctrl_wide <- .cv_ctrl(parallel = TRUE,  seed = 1L, warm_start = FALSE, workers = 2L)
  fold_obj  <- bhcrr_make_folds(y, ctrl_seq)

  r_seq  <- suppressWarnings(bhcrr_cv(x, y, s0_small, s1_small, ctrl_seq,  folds = fold_obj))
  r_wide <- suppressWarnings(bhcrr_cv(x, y, s0_small, s1_small, ctrl_wide, folds = fold_obj))

  expect_equal(r_seq$tuning,      r_wide$tuning)
  expect_equal(r_seq$pooled,      r_wide$pooled)
  expect_equal(r_seq$fold_scores, r_wide$fold_scores)
})

test_that("$timing reports actual workers, parallel flag, and dispatch label", {
  fix <- .load_cv_fixture()
  x <- fix$x; y <- fix$y
  ctrl <- .cv_ctrl(parallel = FALSE)
  r <- suppressWarnings(bhcrr_cv(x, y, s0_small, s1_small, ctrl))

  expect_equal(r$timing$workers,  1L)
  expect_false(r$timing$parallel)
  expect_equal(r$timing$dispatch, "fold")

  ctrl_cold <- .cv_ctrl(parallel = FALSE, warm_start = FALSE)
  r2 <- suppressWarnings(bhcrr_cv(x, y, s0_small, s1_small, ctrl_cold))
  expect_equal(r2$timing$dispatch, "fold-pair")
})

test_that("fold_inits supplied + parallel: phase 1 skipped, init method not called", {
  skip_parallel()
  fix <- .load_cv_fixture()
  x <- fix$x; y <- fix$y

  call_count <- 0L
  counting_init <- function(x, y, ...) {
    call_count <<- call_count + 1L
    list(init = rep(0, ncol(x)), meta = NULL)
  }
  ctrl_seq <- .cv_ctrl(init_method = counting_init, parallel = FALSE, seed = 1L)
  fold_obj <- bhcrr_make_folds(y, ctrl_seq)
  r1 <- suppressWarnings(bhcrr_cv(x, y, s0_small, s1_small, ctrl_seq, folds = fold_obj))

  # Second run with parallel and supplied fold_inits
  call_count <- 0L
  counting_init2 <- function(x, y, ...) {
    call_count <<- call_count + 1L
    list(init = rep(0, ncol(x)), meta = NULL)
  }
  ctrl_par <- .cv_ctrl(init_method = counting_init2, parallel = TRUE, workers = 2L, seed = 1L)
  r2 <- suppressWarnings(
    bhcrr_cv(x, y, s0_small, s1_small, ctrl_par,
             folds = fold_obj, fold_inits = r1$fold_inits)
  )
  expect_equal(call_count, 0L,
    label = "init method must not be called when fold_inits is supplied (parallel path)")
})
