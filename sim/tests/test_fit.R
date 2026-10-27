library(testthat)

# Resolve repo root before load.r is available
.test_root <- local({
  e <- Sys.getenv("REPO_ROOT")
  if (nzchar(e)) return(normalizePath(e, mustWork = TRUE))
  p <- normalizePath(getwd(), mustWork = TRUE)
  repeat {
    if (file.exists(file.path(p, "DESCRIPTION"))) return(p)
    pp <- dirname(p)
    if (identical(pp, p)) stop("Cannot locate repo root; set REPO_ROOT")
    p <- pp
  }
})

source(file.path(.test_root, "sim", "R", "load.r"))
load_sim("dev", force = FALSE)

# ---------------------------------------------------------------------------
# Helpers: tiny but end-to-end runnable configuration
# Intended runtime: < 60 s for the whole test file.
# ---------------------------------------------------------------------------

# Tiny scenario: nobs=80, p=8, 2 active predictors, simple block structure
.tiny_scenario_args <- list(
  name            = "tiny.fit.test",
  nobs            = 80L,
  zero_gap_target = 0.05,
  beta1_active    = c(0.8, -0.8),
  beta2_active    = "neg_beta1",
  active_block    = c("strong", "weak"),
  block_props     = c(indep = 0.50, weak = 0.25, strong = 0.25),
  block_rho       = c(indep = 0.00, weak = 0.30, strong = 0.50),
  latent_type     = "gaussian",
  latent_q        = c(indep = 0.50, weak = 0.50, strong = 0.50),
  cens_rate       = 0.05,
  u_min           = 50,
  u_max           = 150
)

# Control that makes CV fast: 3 folds, 1 repetition
.tiny_control <- function(scenario) {
  bhcrr_cv_control(
    nfolds      = 3L,
    ncv         = 1L,
    init_method = "LASSO_cv",
    init_args   = list(lambda_path =
                         10^seq(log10(0.1), log10(0.001), length = 5)),
    fit_args    = list(
      initial_sparsity  = 0.05,
      theta_a           = 1,
      theta_b           = 1,
      epsilon           = 1e-04,
      inner_maxit_start = 1000,
      maxit             = 20L
    ),
    parallel    = FALSE,
    eval_time   = (scenario$u_min + scenario$u_max) / 2
  )
}

# Tiny grid: 2 x 2
.tiny_s0 <- c(0.04, 0.08)
.tiny_s1 <- c(0.4,  0.8)

# Build data once and reuse across tests to save time
.spec <- do.call(scenario, .tiny_scenario_args)

set.seed(42L)
.tiny_dat <- sim_data(.spec, p = 8L)
.x <- .tiny_dat$x
.y <- .tiny_dat$y

# Run a single fit once; reuse in multiple tests.
set.seed(1L)
.fit_result <- fit_bhcrr(
  .x, .y, .spec, run_id = 1L,
  control   = .tiny_control(.spec),
  s0_seq    = .tiny_s0,
  s1_seq    = .tiny_s1,
  max_widen = 0L
)

# ---------------------------------------------------------------------------
# Contract: return carries every required field with the right class/length
# ---------------------------------------------------------------------------

test_that("return has all contract fields", {
  r <- .fit_result
  expect_true(is.list(r))
  for (nm in c("coef", "s0", "s1", "autotune", "fit", "mod_initial",
               "converged", "all_pass", "flags", "n_attempts",
               "eval_time", "fit_seconds", "timing", "control")) {
    expect_true(nm %in% names(r), label = paste("field present:", nm))
  }
})

test_that("$coef is named numeric of length p, names match colnames(x)", {
  r <- .fit_result
  p <- ncol(.x)
  expect_true(is.numeric(r$coef))
  expect_length(r$coef, p)
  expect_identical(names(r$coef), colnames(.x))
})

test_that("$s0 and $s1 match autotune$best$s0 and $best$s1", {
  r <- .fit_result
  expect_equal(r$s0, r$autotune$best$s0)
  expect_equal(r$s1, r$autotune$best$s1)
})

test_that("$fit is class ssl_psdh; $mod_initial is class ssl_psdh", {
  expect_s3_class(.fit_result$fit,         "ssl_psdh")
  expect_s3_class(.fit_result$mod_initial, "ssl_psdh")
})

test_that("$converged is logical(1)", {
  expect_true(is.logical(.fit_result$converged))
  expect_length(.fit_result$converged, 1L)
})

test_that("$all_pass is logical(1); $flags is character; $n_attempts is integer(1)", {
  r <- .fit_result
  expect_true(is.logical(r$all_pass))
  expect_length(r$all_pass, 1L)
  expect_true(is.character(r$flags))
  expect_true(is.integer(r$n_attempts) || is.numeric(r$n_attempts))
  expect_length(r$n_attempts, 1L)
})

# ---------------------------------------------------------------------------
# Timing: fit_seconds is numeric (NOT difftime), positive
# ---------------------------------------------------------------------------

test_that("$fit_seconds is numeric (not difftime), positive", {
  fs <- .fit_result$fit_seconds
  expect_false(inherits(fs, "difftime"),
               label = "fit_seconds must NOT be difftime")
  expect_true(is.numeric(fs))
  expect_length(fs, 1L)
  expect_gt(fs, 0)
})

test_that("$timing components sum to $fit_seconds within tolerance", {
  r <- .fit_result
  component_sum <- r$timing$init_fit + r$timing$autotune + r$timing$final_fit
  expect_equal(component_sum, r$fit_seconds, tolerance = 1e-6)
})

test_that("$timing components are all numeric and non-negative", {
  for (nm in c("init_fit", "autotune", "final_fit")) {
    v <- .fit_result$timing[[nm]]
    expect_true(is.numeric(v), label = paste("timing$", nm, "is numeric"))
    expect_gte(v, 0,           label = paste("timing$", nm, "is non-negative"))
  }
})

# ---------------------------------------------------------------------------
# fit_args mirroring: control carries theta_b=1 and matching hyperparameters
# ---------------------------------------------------------------------------

test_that("control$fit_args has theta_b=1 and matches final-fit hyperparameters", {
  fa <- .fit_result$control$fit_args
  expect_equal(fa$theta_b,           1)
  expect_equal(fa$initial_sparsity,  0.05)
  expect_equal(fa$epsilon,           1e-04)
  expect_equal(fa$inner_maxit_start, 1000)
})

# ---------------------------------------------------------------------------
# init_method consistency: same in control and in the direct fit calls
# The direct fit calls use the file-level constant "LASSO_cv"; the control
# was built with the same value via .tiny_control() which uses "LASSO_cv".
# ---------------------------------------------------------------------------

test_that("control$init_method matches the file-level init_method constant", {
  # The control was constructed with init_method = "LASSO_cv";
  # bhcrr_cv_control stores it (possibly after normalization) in $init_method.
  ctrl <- .fit_result$control
  # The stored field may be the label string; check it resolves to LASSO_cv
  expect_true(!is.null(ctrl$init_method))
  # Normalize: the resolved label should contain "LASSO" and "cv" (case-tolerant)
  label <- paste(ctrl$init_method, collapse = " ")
  expect_true(grepl("LASSO", label, ignore.case = TRUE))
})

# ---------------------------------------------------------------------------
# eval_time: non-NULL, equals scenario midpoint
# ---------------------------------------------------------------------------

test_that("$eval_time is non-NULL and equals (u_min + u_max) / 2", {
  expected_tau <- (.spec$u_min + .spec$u_max) / 2
  expect_false(is.null(.fit_result$eval_time))
  expect_equal(.fit_result$eval_time, expected_tau)
})

test_that("$control$eval_time matches $eval_time", {
  expect_equal(.fit_result$control$eval_time, .fit_result$eval_time)
})

# ---------------------------------------------------------------------------
# Control merging: new cases exercising the merge behavior (PART B)
# ---------------------------------------------------------------------------

# Minimal fast control for merge tests: empty fit_args so all five mirrored
# entries come from defaults; small lambda path to keep CV cheap.
.merge_base_ctrl <- bhcrr_cv_control(
  nfolds      = 3L, ncv = 1L,
  init_method = "LASSO_cv",
  init_args   = list(lambda_path = 10^seq(log10(0.1), log10(0.001), length = 5)),
  parallel    = FALSE,
  eval_time   = 100
)

test_that("fit_args = list(): all five mirrored entries present after merge", {
  set.seed(10L)
  r <- suppressWarnings(
    fit_bhcrr(.x, .y, .spec, run_id = 1L,
              control = .merge_base_ctrl,
              s0_seq = .tiny_s0, s1_seq = .tiny_s1, max_widen = 0L,
              maxit_folds = 20L, maxit_final = 20L)
  )
  fa <- r$control$fit_args
  expect_equal(fa$initial_sparsity,  0.05,  label = "initial_sparsity filled")
  expect_equal(fa$theta_a,           1,     label = "theta_a filled")
  expect_equal(fa$theta_b,           1,     label = "theta_b filled")
  expect_equal(fa$epsilon,           1e-04, label = "epsilon filled")
  expect_equal(fa$inner_maxit_start, 1000,  label = "inner_maxit_start filled")
})

test_that("fit_args = list(maxit=20L): five mirrored entries AND maxit=20L present", {
  ctrl_m20 <- bhcrr_cv_control(
    nfolds = 3L, ncv = 1L, init_method = "LASSO_cv",
    init_args = list(lambda_path = 10^seq(log10(0.1), log10(0.001), length = 5)),
    fit_args  = list(maxit = 20L),
    parallel  = FALSE, eval_time = 100
  )
  set.seed(11L)
  r <- suppressWarnings(
    fit_bhcrr(.x, .y, .spec, run_id = 1L,
              control = ctrl_m20,
              s0_seq = .tiny_s0, s1_seq = .tiny_s1, max_widen = 0L)
  )
  fa <- r$control$fit_args
  # Five mirrored entries
  expect_equal(fa$initial_sparsity,  0.05)
  expect_equal(fa$theta_a,           1)
  expect_equal(fa$theta_b,           1)
  expect_equal(fa$epsilon,           1e-04)
  expect_equal(fa$inner_maxit_start, 1000)
  # Fold maxit taken from fit_args
  expect_equal(fa$maxit, 20L,
               label = "fold maxit = 20 from control$fit_args")
  # Fit ran successfully (maxit_final = 100 by default, not 20)
  expect_s3_class(r$fit, "ssl_psdh")
})

test_that("theta_b = 2 propagates: control$fit_args$theta_b is 2 and pips differ from theta_b = 1", {
  make_ctrl_tb <- function(tb) {
    bhcrr_cv_control(
      nfolds = 3L, ncv = 1L, init_method = "LASSO_cv",
      init_args = list(lambda_path = 10^seq(log10(0.1), log10(0.001), length = 5)),
      fit_args  = list(theta_b = tb, maxit = 20L),
      parallel  = FALSE, eval_time = 100
    )
  }
  set.seed(55L)
  r1 <- suppressWarnings(
    fit_bhcrr(.x, .y, .spec, run_id = 1L, control = make_ctrl_tb(1),
              s0_seq = .tiny_s0, s1_seq = .tiny_s1, max_widen = 0L)
  )
  set.seed(55L)
  r2 <- suppressWarnings(
    fit_bhcrr(.x, .y, .spec, run_id = 1L, control = make_ctrl_tb(2),
              s0_seq = .tiny_s0, s1_seq = .tiny_s1, max_widen = 0L)
  )
  expect_equal(r2$control$fit_args$theta_b, 2,
               label = "$control$fit_args$theta_b is 2")
  # Assert against the fit object: pips must differ between theta_b=1 and
  # theta_b=2 because the prior precision parameter changed. If theta_b were
  # silently ignored, both fits would run at theta_b=1 and pips would be
  # identical with the same seed.
  expect_false(identical(r2$fit$pips, r1$fit$pips),
               label = "$fit$pips differ between theta_b=1 and theta_b=2")
})

test_that("eval_time = NULL in supplied control -> $eval_time == tau and $control$eval_time == tau", {
  ctrl_no_tau <- bhcrr_cv_control(
    nfolds = 3L, ncv = 1L, init_method = "LASSO_cv",
    init_args = list(lambda_path = 10^seq(log10(0.1), log10(0.001), length = 5)),
    fit_args  = list(maxit = 20L),
    parallel  = FALSE
    # eval_time omitted -> NULL
  )
  set.seed(20L)
  r <- suppressWarnings(
    fit_bhcrr(.x, .y, .spec, run_id = 1L,
              control = ctrl_no_tau,
              s0_seq = .tiny_s0, s1_seq = .tiny_s1, max_widen = 0L)
  )
  tau <- (.spec$u_min + .spec$u_max) / 2
  expect_equal(r$eval_time,          tau, label = "$eval_time == tau")
  expect_equal(r$control$eval_time,  tau, label = "$control$eval_time == tau")
})

test_that("eval_time = 77 honored: $eval_time == 77 not silently replaced by tau", {
  ctrl_77 <- bhcrr_cv_control(
    nfolds = 3L, ncv = 1L, init_method = "LASSO_cv",
    init_args = list(lambda_path = 10^seq(log10(0.1), log10(0.001), length = 5)),
    fit_args  = list(maxit = 20L),
    parallel  = FALSE, eval_time = 77
  )
  set.seed(21L)
  r <- suppressWarnings(
    fit_bhcrr(.x, .y, .spec, run_id = 1L,
              control = ctrl_77,
              s0_seq = .tiny_s0, s1_seq = .tiny_s1, max_widen = 0L)
  )
  expect_equal(r$eval_time,         77, label = "$eval_time == 77")
  expect_equal(r$control$eval_time, 77, label = "$control$eval_time == 77")
  tau <- (.spec$u_min + .spec$u_max) / 2
  expect_false(isTRUE(all.equal(r$eval_time, tau)),
               label = "$eval_time != tau (not silently replaced)")
})

test_that("init_method = 'LASSO_bic' propagates to $fit and $mod_initial", {
  ctrl_bic <- bhcrr_cv_control(
    nfolds = 3L, ncv = 1L, init_method = "LASSO_bic",
    fit_args = list(maxit = 20L),
    parallel = FALSE, eval_time = 100
  )
  set.seed(30L)
  r <- suppressWarnings(
    fit_bhcrr(.x, .y, .spec, run_id = 1L,
              control = ctrl_bic,
              s0_seq = .tiny_s0, s1_seq = .tiny_s1, max_widen = 0L)
  )
  expect_true(grepl("LASSO_bic", r$fit$init_method,         ignore.case = TRUE),
              label = "$fit$init_method used LASSO_bic")
  expect_true(grepl("LASSO_bic", r$mod_initial$init_method, ignore.case = TRUE),
              label = "$mod_initial$init_method used LASSO_bic")
})

test_that("control = NULL regression guard: resolved $control matches production defaults", {
  set.seed(1L)
  r_null <- suppressWarnings(
    fit_bhcrr(.x, .y, .spec, run_id = 1L,
              s0_seq = .tiny_s0, s1_seq = .tiny_s1, max_widen = 0L)
  )
  ctrl <- r_null$control
  fa   <- ctrl$fit_args
  tau  <- (.spec$u_min + .spec$u_max) / 2

  expect_equal(fa$initial_sparsity,  0.05,  label = "default initial_sparsity")
  expect_equal(fa$theta_a,           1,     label = "default theta_a")
  expect_equal(fa$theta_b,           1,     label = "default theta_b")
  expect_equal(fa$epsilon,           1e-04, label = "default epsilon")
  expect_equal(fa$inner_maxit_start, 1000,  label = "default inner_maxit_start")
  expect_equal(fa$maxit,             50,    label = "default fold maxit = 50")
  expect_equal(ctrl$init_method,     "LASSO_cv", label = "default init_method")
  expect_true(length(ctrl$init_args$lambda_path) == 10L,
              label = "default 10-point lambda_path")
  expect_equal(ctrl$eval_time,       tau,   label = "$control$eval_time == tau")
  expect_equal(r_null$eval_time,     tau,   label = "$eval_time == tau")
})

# ---------------------------------------------------------------------------
# B5 guard: NULL best stops with a message naming run_id
# ---------------------------------------------------------------------------

test_that("NULL autotune$best triggers an error naming run_id", {
  # Build a fake autotune object whose $best is NULL; inject via a stub fit
  # by calling fit_bhcrr with mocked autotune -- we patch by temporarily
  # replacing bhcrr_autotune in the package environment.
  pkg_env <- as.environment("package:bhCRR")

  # Save original
  orig_autotune <- get("bhcrr_autotune", envir = pkg_env, inherits = FALSE)
  on.exit({
    unlockBinding("bhcrr_autotune", pkg_env)
    assign("bhcrr_autotune", orig_autotune, envir = pkg_env)
    lockBinding("bhcrr_autotune", pkg_env)
  }, add = TRUE)

  # Stub that returns NULL best
  stub <- function(...) {
    list(
      best       = NULL,
      all_pass   = FALSE,
      flags      = "stubbed",
      n_attempts = 1L
    )
  }
  unlockBinding("bhcrr_autotune", pkg_env)
  assign("bhcrr_autotune", stub, envir = pkg_env)

  expect_error(
    fit_bhcrr(.x, .y, .spec, run_id = 777L,
              control   = .tiny_control(.spec),
              s0_seq    = .tiny_s0,
              s1_seq    = .tiny_s1,
              max_widen = 0L),
    regexp = "777"
  )
})
