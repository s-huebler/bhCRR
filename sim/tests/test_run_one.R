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
# Fixtures -- built once. Target: < 2 minutes total for this file.
# n=80, p=10, 2 active predictors; 2x2 grid; max_widen=0; 3 folds, 1 rep.
# ---------------------------------------------------------------------------

.sargs <- list(
  name            = "run1.test",
  nobs            = 80L,
  zero_gap_target = 0.05,
  beta1_active    = c(1.0, -1.0),
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
.spec <- do.call(scenario, .sargs)

.tiny_ctrl <- bhcrr_cv_control(
  nfolds      = 3L, ncv = 1L,
  init_method = "LASSO_cv",
  init_args   = list(lambda_path = 10^seq(log10(0.1), log10(0.001), length = 5)),
  fit_args    = list(initial_sparsity = 0.05, theta_a = 1, theta_b = 1,
                     epsilon = 1e-04, inner_maxit_start = 1000, maxit = 20L),
  parallel    = FALSE,
  eval_time   = (.spec$u_min + .spec$u_max) / 2
)

.r1 <- suppressWarnings(
  run_one(.spec, p = 10L, run_id = 1L,
          fit_control = .tiny_ctrl,
          s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8), max_widen = 0L)
)

# ---------------------------------------------------------------------------
# Top-level structure
# ---------------------------------------------------------------------------

test_that("result has every documented top-level field", {
  expected <- c("run_id","p","scenario","data","fit","comparators",
                "fg_oracle","scores","meta")
  for (nm in expected)
    expect_true(nm %in% names(.r1), label = paste("top-level field", nm))
})

test_that("$run_id and $p are integers with correct values", {
  expect_identical(.r1$run_id, 1L)
  expect_identical(.r1$p,      10L)
})

test_that("$scenario is a validated bhcrr_scenario", {
  expect_true(inherits(.r1$scenario, "bhcrr_scenario"))
})

# ---------------------------------------------------------------------------
# $data
# ---------------------------------------------------------------------------

test_that("$data has exactly x, y, beta1, beta2 and nothing else", {
  expect_identical(sort(names(.r1$data)), sort(c("x","y","beta1","beta2")))
})

test_that("$data$x is n x p matrix with colnames X_1..X_p", {
  n <- .spec$nobs; p <- 10L
  expect_equal(dim(.r1$data$x), c(n, p))
  expect_identical(colnames(.r1$data$x), paste0("X_", seq_len(p)))
})

test_that("block_map and block_signature are NOT in $data", {
  expect_false("block_map"       %in% names(.r1$data))
  expect_false("block_signature" %in% names(.r1$data))
  expect_false("design_settings" %in% names(.r1$data))
})

# ---------------------------------------------------------------------------
# $meta -- every documented field
# ---------------------------------------------------------------------------

test_that("$meta has every documented field", {
  expected <- c("run_id","p","scenario_name","seed","git_sha",
                "bhcrr_version","r_version","elapsed_seconds",
                "block_map","block_signature","design_settings")
  for (nm in expected)
    expect_true(nm %in% names(.r1$meta), label = paste("$meta field", nm))
})

test_that("$meta$elapsed_seconds is numeric(1), positive, not difftime", {
  es <- .r1$meta$elapsed_seconds
  expect_false(inherits(es, "difftime"), label = "not difftime")
  expect_true(is.numeric(es),           label = "numeric")
  expect_equal(length(es), 1L,          label = "length 1")
  expect_gt(es, 0,                      label = "positive")
})

test_that("$meta$seed is integer(1) when seed = TRUE", {
  expect_true(is.integer(.r1$meta$seed), label = "integer")
  expect_equal(length(.r1$meta$seed), 1L)
  expect_false(is.na(.r1$meta$seed),   label = "not NA")
})

test_that("$meta$git_sha is a non-empty character(1)", {
  expect_true(is.character(.r1$meta$git_sha))
  expect_equal(length(.r1$meta$git_sha), 1L)
  expect_true(nzchar(.r1$meta$git_sha))
})

test_that("$meta$bhcrr_version and $meta$r_version are non-empty character(1)", {
  expect_true(is.character(.r1$meta$bhcrr_version))
  expect_true(nzchar(.r1$meta$bhcrr_version))
  expect_true(is.character(.r1$meta$r_version))
  expect_true(nzchar(.r1$meta$r_version))
})

# ---------------------------------------------------------------------------
# $meta$block_map
# ---------------------------------------------------------------------------

test_that("$meta$block_map has one row per predictor with X_1..X_p names", {
  bm <- .r1$meta$block_map
  p  <- 10L
  expect_equal(nrow(bm), p)
  expect_identical(bm$name, paste0("X_", seq_len(p)))
})

test_that("$meta$block_map has the required columns and correct types", {
  bm <- .r1$meta$block_map
  expect_true(is.character(bm$name))
  expect_true(is.character(bm$block))
  expect_true(is.numeric(bm$block_rho))
  expect_true(is.logical(bm$designed))
  expect_true(is.logical(bm$active1))
  expect_true(is.logical(bm$active2))
})

# ---------------------------------------------------------------------------
# $meta$block_signature
# ---------------------------------------------------------------------------

test_that("$meta$block_signature matches sim_data() for same scenario and p", {
  set.seed(.r1$meta$seed)
  ref <- sim_data(.spec, p = 10L)
  expect_identical(.r1$meta$block_signature, ref$block_signature)
})

# ---------------------------------------------------------------------------
# $meta$design_settings
# ---------------------------------------------------------------------------

test_that("$meta$design_settings has expected fields including event_mix", {
  ds <- .r1$meta$design_settings
  expected <- c("nobs","npredictors","block_props","latent","latent_balanced",
                "unit_variance","layout","event_mix")
  for (nm in expected)
    expect_true(nm %in% names(ds), label = paste("design_settings field", nm))
})

test_that("$meta$design_settings$event_mix sums to 1", {
  em <- .r1$meta$design_settings$event_mix
  expect_true(is.numeric(em))
  expect_equal(sum(em), 1, tolerance = 1e-9)
})

# ---------------------------------------------------------------------------
# $scores
# ---------------------------------------------------------------------------

test_that("$scores$overall has 6 rows in documented model order", {
  ov <- .r1$scores$overall
  expect_equal(nrow(ov), 6L)
  expect_equal(ov$model, c("bhcrr","fg_oracle","lasso","alasso","scad","mcp"))
})

test_that("$scores$by_block has 6 * nblocks rows", {
  nb <- length(unique(.r1$meta$block_map$block))
  expect_equal(nrow(.r1$scores$by_block), 6L * nb)
})

# ---------------------------------------------------------------------------
# $comparators and $fg_oracle are separate top-level fields
# ---------------------------------------------------------------------------

test_that("$comparators retains setup_seconds and standardization", {
  expect_true("setup_seconds"    %in% names(.r1$comparators))
  expect_true("standardization"  %in% names(.r1$comparators))
  for (nm in c("lasso","alasso","scad","mcp"))
    expect_true(nm %in% names(.r1$comparators),
                label = paste("comparators has", nm))
})

test_that("$fg_oracle is a separate field and not inside $comparators", {
  expect_true("fg_oracle" %in% names(.r1),
              label = "$fg_oracle at top level")
  expect_false("fg_oracle" %in% names(.r1$comparators),
               label = "$fg_oracle not inside $comparators")
})

# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------

test_that("same (scenario, p, run_id) produces identical $data$x", {
  r1b <- suppressWarnings(
    run_one(.spec, p = 10L, run_id = 1L,
            fit_control = .tiny_ctrl,
            s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8), max_widen = 0L)
  )
  expect_identical(.r1$data$x, r1b$data$x)
})

test_that("same (scenario, p, run_id) produces identical $scores$overall non-timing columns", {
  r1b <- suppressWarnings(
    run_one(.spec, p = 10L, run_id = 1L,
            fit_control = .tiny_ctrl,
            s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8), max_widen = 0L)
  )
  # fit_seconds is wall clock and varies between calls; exclude it
  drop <- c("fit_seconds", "elapsed_seconds")
  cols <- setdiff(names(.r1$scores$overall), drop)
  expect_equal(.r1$scores$overall[, cols], r1b$scores$overall[, cols])
})

test_that("different run_id produces different $data$x", {
  r2 <- suppressWarnings(
    run_one(.spec, p = 10L, run_id = 2L,
            fit_control = .tiny_ctrl,
            s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8), max_widen = 0L)
  )
  expect_false(identical(.r1$data$x, r2$data$x),
               label = "run_id=1 and run_id=2 produce different x")
})

test_that("same run_id at different p produces different $data$x", {
  r1_p25 <- suppressWarnings(
    run_one(.spec, p = 25L, run_id = 1L,
            fit_control = bhcrr_cv_control(
              nfolds = 3L, ncv = 1L, init_method = "LASSO_cv",
              init_args = list(lambda_path = 10^seq(log10(0.1), log10(0.001), length = 5)),
              fit_args  = list(initial_sparsity = 0.05, theta_a = 1, theta_b = 1,
                               epsilon = 1e-04, inner_maxit_start = 1000, maxit = 20L),
              parallel  = FALSE,
              eval_time = (.spec$u_min + .spec$u_max) / 2),
            s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8), max_widen = 0L)
  )
  # $data$x has different number of columns; compare the shared columns
  shared <- intersect(colnames(.r1$data$x), colnames(r1_p25$data$x))
  expect_false(identical(.r1$data$x[, shared, drop = FALSE],
                          r1_p25$data$x[, shared, drop = FALSE]),
               label = "p=10 and p=25 produce different x for shared columns")
})

# ---------------------------------------------------------------------------
# seed = FALSE
# ---------------------------------------------------------------------------

test_that("seed = FALSE sets $meta$seed to NA_integer_", {
  rna <- suppressWarnings(
    run_one(.spec, p = 10L, run_id = 1L, seed = FALSE,
            fit_control = .tiny_ctrl,
            s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8), max_widen = 0L)
  )
  expect_identical(rna$meta$seed, NA_integer_)
})

test_that("two seed = FALSE calls produce different $data$x", {
  set.seed(NULL)  # reset to unpredictable state
  ra <- suppressWarnings(
    run_one(.spec, p = 10L, run_id = 1L, seed = FALSE,
            fit_control = .tiny_ctrl,
            s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8), max_widen = 0L)
  )
  rb <- suppressWarnings(
    run_one(.spec, p = 10L, run_id = 1L, seed = FALSE,
            fit_control = .tiny_ctrl,
            s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8), max_widen = 0L)
  )
  expect_false(identical(ra$data$x, rb$data$x),
               label = "two unseeded calls produce different x")
})

# ---------------------------------------------------------------------------
# B5: guard against fit_control$seed with seed = TRUE
# ---------------------------------------------------------------------------

test_that("B5 guard errors when fit_control carries non-NULL seed and seed = TRUE", {
  ctrl_with_seed <- .tiny_ctrl
  ctrl_with_seed$seed <- 42L
  expect_error(
    run_one(.spec, p = 10L, run_id = 1L, seed = TRUE,
            fit_control = ctrl_with_seed,
            s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8), max_widen = 0L),
    regexp = "fit_control\\$seed"
  )
})

# ---------------------------------------------------------------------------
# Scenario resolution
# ---------------------------------------------------------------------------

test_that("a character path to a scenarios file is resolved correctly", {
  path <- file.path(.test_root, "sim", "scenarios", "run1.test.R")
  # Write a temp scenario file if it doesn't exist; use .spec values
  if (!file.exists(path)) {
    writeLines(c(
      'scenario <- list(',
      '  name            = "run1.test",',
      '  nobs            = 80L,',
      '  zero_gap_target = 0.05,',
      '  beta1_active    = c(1.0, -1.0),',
      '  beta2_active    = "neg_beta1",',
      '  active_block    = c("strong", "weak"),',
      '  block_props     = c(indep = 0.50, weak = 0.25, strong = 0.25),',
      '  block_rho       = c(indep = 0.00, weak = 0.30, strong = 0.50),',
      '  latent_type     = "gaussian",',
      '  latent_q        = c(indep = 0.50, weak = 0.50, strong = 0.50),',
      '  cens_rate       = 0.05,',
      '  u_min           = 50,',
      '  u_max           = 150',
      ')'),
      path)
  }
  r_path <- suppressWarnings(
    run_one(path, p = 10L, run_id = 1L,
            fit_control = .tiny_ctrl,
            s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8), max_widen = 0L)
  )
  expect_true(inherits(r_path$scenario, "bhcrr_scenario"))
  expect_equal(r_path$scenario$name, "run1.test")
})

test_that("an unvalidated plain list is rejected with an informative error", {
  bare_list <- as.list(.spec)
  class(bare_list) <- "list"   # strip bhcrr_scenario class
  expect_error(
    run_one(bare_list, p = 10L, run_id = 1L, fit_control = .tiny_ctrl,
            s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8), max_widen = 0L),
    regexp = "bhcrr_scenario"
  )
})

# ---------------------------------------------------------------------------
# Error propagation
# ---------------------------------------------------------------------------

test_that("run_one does not swallow errors: impossible p propagates", {
  # p = 0 should cause sim_data to error
  expect_error(
    run_one(.spec, p = 0L, run_id = 1L,
            fit_control = .tiny_ctrl,
            s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8), max_widen = 0L)
  )
})

test_that("run_one does not swallow errors: p smaller than active predictors propagates", {
  # .spec has 2 active predictors; p = 1 < 2 should error in sim_data
  expect_error(
    run_one(.spec, p = 1L, run_id = 1L,
            fit_control = .tiny_ctrl,
            s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8), max_widen = 0L)
  )
})
