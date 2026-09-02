# Tests for R/make_folds.r — bhcrr_make_folds()
#
# Rare-event configuration: n=200, 7.5% cause-1 rate (seed=7).
# At pool=0.01 rsample keeps the stratum; at pool=0.1 (>=7.5%) it pools the
# stratum away and the assignment degrades to an unstratified shuffle.

.make_rare_event_y <- function(n = 200, seed = 7) {
  set.seed(seed)
  cbind(
    time   = rexp(n) + 0.01,
    status = sample(c(0L, 1L, 2L), n, replace = TRUE,
                    prob = c(0.72, 0.08, 0.20))
  )
}

# ---------------------------------------------------------------------------
# 1. Dimensions and coverage — every observation in exactly one fold per rep
# ---------------------------------------------------------------------------

test_that("unstratified folds have correct dimensions and full coverage", {
  y    <- .make_rare_event_y()
  ctrl <- bhcrr_cv_control(nfolds = 5L, ncv = 3L, strata = "none")
  f    <- bhcrr_make_folds(y, ctrl)

  expect_equal(dim(f$foldid), c(nrow(y), 3L))
  expect_equal(f$nfolds, 5L)
  expect_equal(f$ncv,    3L)
  expect_equal(f$strata, "none")

  for (k in seq_len(f$ncv)) {
    fold_col <- f$foldid[, k]
    expect_true(all(!is.na(fold_col)))
    expect_setequal(unique(fold_col), seq_len(f$nfolds))
    # every obs in exactly one fold
    for (j in seq_len(f$nfolds)) {
      expect_true(sum(fold_col == j) > 0L)
    }
    expect_equal(length(fold_col), nrow(y))
  }
})

test_that("fold_event_counts matrix has correct dimensions and positive sums", {
  y    <- .make_rare_event_y()
  ctrl <- bhcrr_cv_control(nfolds = 5L, ncv = 2L, strata = "none")
  f    <- bhcrr_make_folds(y, ctrl)
  expect_equal(dim(f$fold_event_counts), c(2L, 5L))
  expect_equal(sum(f$fold_event_counts), sum(y[, 2] == 1L) * 2L)
})

# ---------------------------------------------------------------------------
# 2. Stratified folds: pool=0.01 keeps the stratum (spread <= 1)
# ---------------------------------------------------------------------------

test_that("strata='cause1' pool=0.01 gives spread <= 1 across all reps", {
  y    <- .make_rare_event_y()
  ctrl <- bhcrr_cv_control(nfolds = 5L, ncv = 4L,
                            strata = "cause1", pool = 0.01)
  # rsample may emit a once-per-session advisory about small strata; that is
  # expected and not the diagnostic warning we are guarding against.
  f <- suppressWarnings(bhcrr_make_folds(y, ctrl))

  counts <- f$fold_event_counts
  spread <- max(counts) - min(counts)
  expect_lte(spread, 1L,
    label = "per-fold cause-1 spread at pool=0.01")
})

# ---------------------------------------------------------------------------
# 3. pool=0.1 with a <10% event rate degrades stratification — documented trap
# ---------------------------------------------------------------------------
# rsample pools any stratum below `pool` fraction into the remainder, so a
# 7.5%-event dataset with pool=0.1 silently produces an unstratified shuffle.
# This test documents that trap: spread must exceed 1 and the diagnostic
# warning must fire.

test_that("pool=0.1 on a <10%-event dataset degrades stratification (documents rsample trap)", {
  y <- .make_rare_event_y()
  # pool=0.1 triggers a warning at bhcrr_cv_control() construction
  ctrl <- suppressWarnings(
    bhcrr_cv_control(nfolds = 5L, ncv = 4L, strata = "cause1", pool = 0.1)
  )
  # bhcrr_make_folds should emit the diagnostic spread warning
  expect_warning(
    { f <- bhcrr_make_folds(y, ctrl) },
    regexp = "uneven|spread|pool"
  )
  spread <- max(f$fold_event_counts) - min(f$fold_event_counts)
  expect_gt(spread, 1L,
    label = "spread at pool=0.1 should exceed 1 (poor stratification)")
})

# ---------------------------------------------------------------------------
# 4. Supplied foldid passes through unchanged
# ---------------------------------------------------------------------------

test_that("user-supplied foldid is returned as-is", {
  y   <- .make_rare_event_y()
  n   <- nrow(y)
  fid <- matrix(rep(1:5, length.out = n), ncol = 1L)
  ctrl <- bhcrr_cv_control(foldid = fid)
  f   <- bhcrr_make_folds(y, ctrl)
  expect_equal(f$foldid, fid)
})

test_that("supplied foldid with wrong nrow errors", {
  y    <- .make_rare_event_y()
  fid  <- matrix(1:10, ncol = 1L)   # nrow=10, not 200
  ctrl <- bhcrr_cv_control(foldid = fid)
  expect_error(bhcrr_make_folds(y, ctrl), regexp = "nrow")
})

# ---------------------------------------------------------------------------
# 5. nfolds > n is capped; nfolds == n forces ncv = 1
# ---------------------------------------------------------------------------

test_that("nfolds > n is capped at n", {
  set.seed(1)
  y_small <- cbind(time = rexp(10) + 0.01,
                   status = c(rep(1L, 3), rep(0L, 7)))
  ctrl <- bhcrr_cv_control(nfolds = 50L, ncv = 2L, strata = "none")
  f    <- suppressWarnings(bhcrr_make_folds(y_small, ctrl))
  # capped at nrow(y_small) = 10; LOO triggers ncv=1 and zero-event warning
  expect_equal(f$nfolds, nrow(y_small))
})

test_that("nfolds == n forces ncv = 1 (LOO)", {
  set.seed(1)
  n <- 15L
  y_loo <- cbind(time   = rexp(n) + 0.01,
                 status = c(rep(1L, 3), rep(0L, n - 3L)))
  ctrl <- bhcrr_cv_control(nfolds = n, ncv = 5L, strata = "none")
  # LOO has n folds with 1 obs each; most folds have 0 cause-1 events —
  # suppress that expected diagnostic to keep the test focused on dimensions.
  f <- suppressWarnings(bhcrr_make_folds(y_loo, ctrl))
  expect_equal(f$ncv,    1L)
  expect_equal(f$nfolds, n)
})

# ---------------------------------------------------------------------------
# 6. Zero cause-1 events in a fold triggers a warning
# ---------------------------------------------------------------------------

test_that("zero cause-1 events in a fold triggers an informative warning", {
  # 2 cause-1 events, 10 folds => at least 8 folds must have 0 events
  set.seed(3)
  n_s <- 40L
  y_sparse <- cbind(
    time   = rexp(n_s) + 0.01,
    status = c(1L, 1L, rep(0L, n_s - 2L))
  )
  ctrl <- bhcrr_cv_control(nfolds = 10L, ncv = 1L, strata = "none")
  expect_warning(
    bhcrr_make_folds(y_sparse, ctrl),
    regexp = "zero cause-1|zero.*fold|fold.*zero"
  )
})
