# Tests for R/autotune.r — bhcrr_autotune(), bhcrr_tune_preflight(),
# bhcrr_tune_validate()
#
# Fixture: small n and p, zero-init, very low maxit, so each test finishes in
# a few seconds even when the retry loop produces a second CV attempt.

# ---------------------------------------------------------------------------
# shared helpers
# ---------------------------------------------------------------------------

.at_data <- function(n = 60, p = 8, seed = 2026L) {
  set.seed(seed)
  x <- matrix(rnorm(n * p), n, p)
  colnames(x) <- paste0("V", seq_len(p))
  time   <- rexp(n, rate = 1) + 0.1
  status <- sample(c(0L, 1L, 2L), n, replace = TRUE, prob = c(0.3, 0.5, 0.2))
  y      <- cbind(time = time, status = status)
  list(x = x, y = y)
}

.at_mod <- function(dat) {
  fit_ssl_psdh(dat$x, dat$y, ss = c(0.03, 0.5),
               init_method = "zero", maxit = 5L)
}

.at_ctrl <- function() {
  bhcrr_cv_control(
    init_method = "zero",
    nfolds      = 3L,
    ncv         = 1L,
    seed        = 42L,
    fit_args    = list(maxit = 5L, epsilon = 1e-3, theta_b = 20L)
  )
}

# ---------------------------------------------------------------------------
# 1. Returns "bhcrr_autotune"; stable-surface fields present and well-typed
# ---------------------------------------------------------------------------

test_that("returns bhcrr_autotune class with all stable-surface fields", {
  dat  <- .at_data()
  mod  <- .at_mod(dat)
  ctrl <- .at_ctrl()

  result <- bhcrr_autotune(mod,
                           s0_seq   = c(0.02, 0.04),
                           s1_seq   = c(0.4,  0.6),
                           control  = ctrl,
                           max_widen = 0L)

  # class
  expect_s3_class(result, "bhcrr_autotune")

  # stable-surface fields present
  expect_true("all_pass"    %in% names(result))
  expect_true("flags"       %in% names(result))
  expect_true("n_attempts"  %in% names(result))
  expect_true("best"        %in% names(result))

  # well-typed
  expect_true(is.logical(result$all_pass))
  expect_length(result$all_pass, 1L)

  expect_true(is.character(result$flags))

  expect_true(is.numeric(result$n_attempts) || is.integer(result$n_attempts))
  expect_length(result$n_attempts, 1L)

  if (!is.null(result$best)) {
    expect_true("s0" %in% names(result$best))
    expect_true("s1" %in% names(result$best))
    expect_true(is.numeric(result$best$s0))
    expect_true(is.numeric(result$best$s1))
    expect_length(result$best$s0, 1L)
    expect_length(result$best$s1, 1L)
  }
})

# ---------------------------------------------------------------------------
# 2. bhcrr_tune_validate() accepts cv$tuning unchanged
# ---------------------------------------------------------------------------

test_that("bhcrr_tune_validate() accepts cv$tuning without modification", {
  dat  <- .at_data()
  mod  <- .at_mod(dat)
  ctrl <- .at_ctrl()

  s0_seq <- c(0.02, 0.04)
  s1_seq <- c(0.4,  0.6)

  cv  <- bhcrr_cv(dat$x, dat$y, s0_seq, s1_seq, control = ctrl)
  est <- mod$coefficients$Estimate
  n   <- nrow(dat$x)
  theta <- min(max(sum(est != 0) / length(est), 1e-4), 1 - 1e-4)

  val <- bhcrr_tune_validate(cv$tuning, n = n, theta = theta,
                              s0_seq = s0_seq, s1_seq = s1_seq)

  expect_true(is.list(val))
  expect_true("checks"       %in% names(val))
  expect_true("all_pass"     %in% names(val))
  expect_true("on_boundary"  %in% names(val))
  expect_true("boundary_dims" %in% names(val))
  expect_true("gap_at_opt"   %in% names(val))
  expect_true(is.data.frame(val$checks))
  expect_true(is.logical(val$all_pass))
})

# ---------------------------------------------------------------------------
# 3. max_widen = 0 produces exactly one attempt
# ---------------------------------------------------------------------------

test_that("max_widen = 0 always produces exactly one attempt", {
  dat  <- .at_data()
  mod  <- .at_mod(dat)
  ctrl <- .at_ctrl()

  result <- bhcrr_autotune(mod,
                           s0_seq    = c(0.02, 0.04),
                           s1_seq    = c(0.4,  0.6),
                           control   = ctrl,
                           max_widen = 0L)

  expect_equal(result$n_attempts, 1L)
  expect_length(result$history, 1L)
})

# ---------------------------------------------------------------------------
# 4. Supplied folds are used in every history element
#
# We use a single-point grid (s0=0.03, s1=0.5) so the optimum is trivially
# on both boundaries of s0 and s1. With max_widen=1 the loop will attempt to
# widen (unless it stalls at the feasibility cap). Regardless of the number of
# attempts, every history element must carry the folds we supplied.
# ---------------------------------------------------------------------------

test_that("supplied folds appear in every history$cv$folds", {
  dat   <- .at_data()
  mod   <- .at_mod(dat)
  ctrl  <- .at_ctrl()
  folds <- bhcrr_make_folds(dat$y, ctrl)

  result <- bhcrr_autotune(mod,
                           s0_seq    = 0.03,
                           s1_seq    = 0.5,
                           control   = ctrl,
                           folds     = folds,
                           max_widen = 1L)

  # Must have at least one attempt.
  expect_gte(result$n_attempts, 1L)

  # Every history element's bhcrr_cv object uses the folds we passed.
  for (h in result$history) {
    expect_identical(h$cv$folds$foldid, folds$foldid)
  }
})

# ---------------------------------------------------------------------------
# 5. Input validation: wrong class stops with a message
# ---------------------------------------------------------------------------

test_that("non-ssl_psdh mod and non-bhcrr_cv_control control produce clear errors", {
  dat  <- .at_data()
  mod  <- .at_mod(dat)
  ctrl <- .at_ctrl()

  # wrong mod class
  bad_mod <- list(x = dat$x, y = dat$y,
                  coefficients = data.frame(Variable = "V1", Estimate = 0))
  expect_error(
    bhcrr_autotune(bad_mod, control = ctrl),
    class = "simpleError"
  )

  # wrong control class
  expect_error(
    bhcrr_autotune(mod, control = list(nfolds = 3L)),
    class = "simpleError"
  )
})
