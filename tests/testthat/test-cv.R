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
                     "score_fold_sd", "n_failed_folds")
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
# 7. Dimensions of $pooled and $fold_scores
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
})
