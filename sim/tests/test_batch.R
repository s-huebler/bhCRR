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
source(file.path(.test_root, "sim", "local", "batch.R"))

# .parse_rx is defined in sim/R/parse.r; load_sim() already sourced it above.
# Re-read it here so this test file is self-contained when run standalone.
.parse_rx <- local({
  pb <- readLines(file.path(.test_root, "sim", "R", "parse.r"))
  rx_line <- grep('^\\.parse_rx <- ', pb, value = TRUE)[1L]
  eval(parse(text = sub('^\\.parse_rx <- ', '', rx_line)))
})

# ---------------------------------------------------------------------------
# Fixtures -- scenario file and control, shared across all tests
# ---------------------------------------------------------------------------

.scenario_file <- file.path(.test_root, "sim", "scenarios", "run1.test.R")

.tiny_ctrl <- bhcrr_cv_control(
  nfolds      = 3L, ncv = 1L,
  init_method = "LASSO_cv",
  init_args   = list(lambda_path = 10^seq(log10(0.1), log10(0.001), length = 5)),
  fit_args    = list(initial_sparsity = 0.05, theta_a = 1, theta_b = 1,
                     epsilon = 1e-04, inner_maxit_start = 1000, maxit = 20L),
  parallel    = FALSE,
  eval_time   = 100
)

# ---------------------------------------------------------------------------
# dryrun -- fast tests; no fitting, no files written
# ---------------------------------------------------------------------------

test_that("dryrun writes no files and creates no manifest", {
  tmp <- tempdir()
  res <- batch_sweep(.scenario_file,
                     npredictors  = 10L,
                     run_start    = 1L, run_end = 2L,
                     results_root = tmp,
                     fit_control  = .tiny_ctrl,
                     s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8),
                     max_widen = 0L,
                     dryrun = TRUE)
  expect_true(res$dryrun)
  manifest <- file.path(tmp, "run1.test", "scenario.env")
  expect_false(file.exists(manifest), label = "dryrun: no manifest")
  expect_equal(length(res$planned), 2L)
  # No files written
  expect_false(any(vapply(res$planned, file.exists, logical(1))),
               label = "dryrun: no Rdata files")
})

test_that("dryrun returns planned paths matching the canonical filename pattern", {
  tmp <- tempdir()
  res <- batch_sweep(.scenario_file,
                     npredictors  = 10L,
                     run_start    = 1L, run_end = 3L,
                     results_root = tmp,
                     fit_control  = .tiny_ctrl,
                     s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8),
                     max_widen = 0L,
                     dryrun = TRUE)
  for (path in res$planned) {
    fn <- basename(path)
    expect_true(grepl(.parse_rx, fn),
                label = paste("filename matches parse_batch regex:", fn))
  }
})

# ---------------------------------------------------------------------------
# Real sweep -- 1p x 2 runs; writes into a fresh tempdir
# ---------------------------------------------------------------------------

# Build the sweep result once; all real-sweep tests reuse it
.tmp_real <- file.path(tempdir(), paste0("batch_test_", Sys.getpid()))
dir.create(.tmp_real, showWarnings = FALSE)

.sweep <- suppressWarnings(
  batch_sweep(.scenario_file,
              npredictors  = 10L,
              run_start    = 1L, run_end = 2L,
              results_root = .tmp_real,
              fit_control  = .tiny_ctrl,
              s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8),
              max_widen = 0L)
)

test_that("real sweep writes exactly the expected canonical paths", {
  p <- 10L; scenario_name <- "run1.test"
  for (run_id in 1L:2L) {
    tag      <- sprintf("p%d_%s_run%d", p, scenario_name, run_id)
    expected <- file.path(.tmp_real, scenario_name, "runs",
                          sprintf("p%d", p), paste0(tag, ".Rdata"))
    expect_true(file.exists(expected),
                label = paste("exists:", basename(expected)))
  }
})

test_that("every written filename matches parse.r's regex", {
  p <- 10L
  runs_p <- file.path(.tmp_real, "run1.test", "runs", sprintf("p%d", p))
  files <- list.files(runs_p, pattern = "\\.Rdata$")
  for (fn in files)
    expect_true(grepl(.parse_rx, fn),
                label = paste("matches regex:", fn))
})

test_that("manifest is written on first run and check_manifest passes on second", {
  manifest <- file.path(.tmp_real, "run1.test", "scenario.env")
  expect_true(file.exists(manifest), label = "manifest written")
  # A second sweep with identical scenario should not error
  expect_no_error(suppressWarnings(
    batch_sweep(.scenario_file,
                npredictors  = 10L,
                run_start    = 3L, run_end = 3L,
                results_root = .tmp_real,
                fit_control  = .tiny_ctrl,
                s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8),
                max_widen = 0L)
  ))
})

test_that("altered scenario hard-stops on manifest mismatch before any fitting", {
  # Write a scenario file with different nobs into a temp path
  alt_path <- tempfile(fileext = ".R")
  writeLines(c(
    'scenario <- list(',
    '  name            = "run1.test",',
    '  nobs            = 60L,',   # different from 80L in the real spec
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
    alt_path)
  # A manifest for run1.test already exists in .tmp_real (written by .sweep)
  expect_error(
    batch_sweep(alt_path,
                npredictors  = 10L,
                run_start    = 99L, run_end = 99L,
                results_root = .tmp_real,
                fit_control  = .tiny_ctrl,
                s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8),
                max_widen = 0L),
    regexp = "nobs|NOBS|mismatch|differ",
    label  = "altered nobs triggers manifest check error"
  )
})

# ---------------------------------------------------------------------------
# Collision guard
# ---------------------------------------------------------------------------

test_that("collision guard refuses re-run of existing indices and names them", {
  # runs 1 and 2 already exist from .sweep
  expect_error(
    batch_sweep(.scenario_file,
                npredictors  = 10L,
                run_start    = 1L, run_end = 2L,
                results_root = .tmp_real,
                fit_control  = .tiny_ctrl,
                s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8),
                max_widen = 0L,
                overwrite = FALSE),
    regexp = "exist|overwrite|collision",
    ignore.case = TRUE
  )
})

test_that("overwrite = TRUE re-runs existing indices without error", {
  expect_no_error(suppressWarnings(
    batch_sweep(.scenario_file,
                npredictors  = 10L,
                run_start    = 1L, run_end = 1L,
                results_root = .tmp_real,
                fit_control  = .tiny_ctrl,
                s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8),
                max_widen = 0L,
                overwrite = TRUE)
  ))
})

test_that("run_start past the collision proceeds without error", {
  expect_no_error(suppressWarnings(
    batch_sweep(.scenario_file,
                npredictors  = 10L,
                run_start    = 10L, run_end = 10L,
                results_root = .tmp_real,
                fit_control  = .tiny_ctrl,
                s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8),
                max_widen = 0L)
  ))
})

# ---------------------------------------------------------------------------
# Failure record: a run_id that run_one() rejects (p = 1 < n_active)
# Inject via a tiny scenario with 2 active predictors but p = 1
# ---------------------------------------------------------------------------

test_that("failing run produces failure record on disk with expected fields; sweep continues", {
  tmp_fail <- file.path(tempdir(), paste0("batch_fail_", Sys.getpid()))
  dir.create(tmp_fail, showWarnings = FALSE)

  # A scenario with 2 active predictors; p = 1 will cause sim_data() to error.
  # scenario_from_file() enforces name == basename; must name the file fail.test.R
  fail_spec_path <- file.path(tempdir(), "fail.test.R")
  writeLines(c(
    'scenario <- list(',
    '  name            = "fail.test",',
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
    fail_spec_path)

  # p = 1 is too small (2 active predictors); sweep must write failure record
  # and then succeed on a second run with a valid p
  # Use a mixed sweep: run 1 will fail (p=1), run 2 will succeed (p=10)
  # but since npredictors is a scalar here we test p=1 alone
  suppressWarnings(suppressMessages(
    batch_sweep(fail_spec_path,
                npredictors  = 1L,  # too small; will error in sim_data
                run_start    = 1L, run_end = 2L,
                results_root = tmp_fail,
                fit_control  = .tiny_ctrl,
                s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8),
                max_widen = 0L)
  ))

  # Both files should exist (failure records)
  for (run_id in 1L:2L) {
    tag  <- sprintf("p1_fail.test_run%d", run_id)
    path <- file.path(tmp_fail, "fail.test", "runs", "p1",
                      paste0(tag, ".Rdata"))
    expect_true(file.exists(path), label = paste("failure record exists:", tag))
    e <- new.env()
    load(path, envir = e)
    expect_true(isTRUE(e$result$failed),
                label = paste(tag, "$failed == TRUE"))
    expect_true(nzchar(e$result$error),
                label = paste(tag, "has non-empty error string"))
    expect_true(is.integer(e$result$meta$seed) || is.na(e$result$meta$seed),
                label = paste(tag, "seed is integer or NA"))
  }
})

# ---------------------------------------------------------------------------
# Memory: returned value does not contain run results
# ---------------------------------------------------------------------------

test_that("batch_sweep return value does not accumulate run result objects", {
  nms <- names(.sweep)
  expect_false("results" %in% nms, label = "no $results list")
  expect_false("data"    %in% nms, label = "no $data in return value")
  expect_false("fit"     %in% nms, label = "no $fit in return value")
  # Return value is the summary metadata only
  expected_fields <- c("scenario_name","npredictors","run_start","run_end",
                       "n_attempted","n_succeeded","n_failed",
                       "elapsed","runs_dir")
  for (nm in expected_fields)
    expect_true(nm %in% nms, label = paste("return field", nm))
})
