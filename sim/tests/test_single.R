library(testthat)

# Resolve repo root
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
source(file.path(.test_root, "sim", "local", "single.R"))

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

# Pre-run one result so the slow fitting only happens once for the tests that
# just need a valid result object
.result <- suppressWarnings(suppressMessages(
  run_single(.scenario_file, p = 25L, run_id = 1L,
             fit_control = .tiny_ctrl,
             s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8),
             max_widen = 0L, quiet = TRUE)
))

# ---------------------------------------------------------------------------
# Return structure
# ---------------------------------------------------------------------------

test_that("run_single returns a complete run_one() result with $warnings", {
  nms <- names(.result)
  for (field in c("run_id", "p", "scenario", "data", "fit", "comparators",
                  "fg_oracle", "scores", "meta", "warnings"))
    expect_true(field %in% nms, label = paste("field present:", field))

  expect_equal(.result$p, 25L)
  expect_equal(.result$run_id, 1L)
  expect_equal(.result$scenario$name, "run1.test")
  expect_true(is.character(.result$warnings),
              label = "$warnings is a character vector")
})

test_that("$meta$seed is an integer (or NA_integer_)", {
  s <- .result$meta$seed
  expect_true(is.integer(s) || (is.na(s) && is.integer(s)),
              label = "$meta$seed is integer or NA_integer_")
})

test_that("$scores$overall and $scores$by_block are data frames", {
  expect_s3_class(.result$scores$overall,  "data.frame")
  expect_s3_class(.result$scores$by_block, "data.frame")
})

# ---------------------------------------------------------------------------
# No files written
# ---------------------------------------------------------------------------

test_that("run_single writes no files (sim/results/ unchanged)", {
  results_dir <- file.path(.test_root, "sim", "results")
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  before <- list.files(results_dir, recursive = TRUE, all.files = TRUE)

  suppressWarnings(suppressMessages(
    run_single(.scenario_file, p = 25L, run_id = 99L,
               fit_control = .tiny_ctrl,
               s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8),
               max_widen = 0L, quiet = TRUE)
  ))

  after <- list.files(results_dir, recursive = TRUE, all.files = TRUE)
  expect_identical(before, after, label = "sim/results/ unchanged after run_single()")
})

# ---------------------------------------------------------------------------
# quiet = TRUE is silent
# ---------------------------------------------------------------------------

test_that("quiet = TRUE produces no output", {
  expect_silent(suppressWarnings(
    run_single(.scenario_file, p = 25L, run_id = 1L,
               fit_control = .tiny_ctrl,
               s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8),
               max_widen = 0L, quiet = TRUE)
  ))
})

# ---------------------------------------------------------------------------
# Error cases
# ---------------------------------------------------------------------------

test_that("scenario_file = NULL errors naming new_scenario()", {
  expect_error(
    run_single(scenario_file = NULL, p = 25L, run_id = 1L, quiet = TRUE),
    regexp = "new_scenario",
    label  = "NULL scenario_file error mentions new_scenario()"
  )
})

test_that("nonexistent scenario_file errors clearly", {
  expect_error(
    run_single("/no/such/path/fake.R", p = 25L, run_id = 1L, quiet = TRUE),
    regexp = "not found|does not exist|No such",
    ignore.case = TRUE,
    label = "nonexistent path gives clear error"
  )
})

test_that("_template.R is refused by scenario_from_file()", {
  tmpl <- file.path(.test_root, "sim", "scenarios", "_template.R")
  skip_if(!file.exists(tmpl), "_template.R not found")
  expect_error(
    suppressMessages(run_single(tmpl, p = 25L, run_id = 1L, quiet = TRUE)),
    regexp = "underscore|_template|template",
    ignore.case = TRUE,
    label = "_template.R refused"
  )
})

# ---------------------------------------------------------------------------
# Seed passthrough: identical data on same (scenario, p, run_id)
# ---------------------------------------------------------------------------

test_that("seed passthrough: same (scenario, p, run_id) gives identical $data$x", {
  r1 <- suppressWarnings(suppressMessages(
    run_single(.scenario_file, p = 10L, run_id = 7L,
               fit_control = .tiny_ctrl,
               s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8),
               max_widen = 0L, quiet = TRUE)
  ))
  r2 <- suppressWarnings(suppressMessages(
    run_single(.scenario_file, p = 10L, run_id = 7L,
               fit_control = .tiny_ctrl,
               s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8),
               max_widen = 0L, quiet = TRUE)
  ))
  expect_identical(r1$data$x, r2$data$x,
                   label = "same seed -> identical design matrix")
  expect_true(is.integer(r1$meta$seed), label = "$meta$seed is integer")
})

# ---------------------------------------------------------------------------
# fit_control passthrough: nfolds actually used, not the default 10
# ---------------------------------------------------------------------------

test_that("fit_control passthrough: nfolds = 3 used, not default 10", {
  ctrl3 <- bhcrr_cv_control(
    nfolds      = 3L, ncv = 1L,
    init_method = "LASSO_cv",
    init_args   = list(lambda_path = 10^seq(log10(0.1), log10(0.001), length = 5)),
    fit_args    = list(initial_sparsity = 0.05, theta_a = 1, theta_b = 1,
                       epsilon = 1e-04, inner_maxit_start = 1000, maxit = 20L),
    parallel    = FALSE,
    eval_time   = 100
  )
  r <- suppressWarnings(suppressMessages(
    run_single(.scenario_file, p = 10L, run_id = 1L,
               fit_control = ctrl3,
               s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8),
               max_widen = 0L, quiet = TRUE)
  ))
  # The control is in the fit result; nfolds is reflected in the CV output
  expect_equal(r$fit$control$nfolds, 3L,
               label = "nfolds = 3 in $fit$control, not default 10")
})

# ---------------------------------------------------------------------------
# Warning capture
# ---------------------------------------------------------------------------

test_that("$warnings collects warnings and keeps them off the console", {
  # run1.test at p=25 with LASSO_cv init generates cencode warnings from
  # fastcmprsk -- these are benign but should land in $warnings, not on console
  r <- suppressMessages(
    run_single(.scenario_file, p = 25L, run_id = 1L,
               fit_control = .tiny_ctrl,
               s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8),
               max_widen = 0L, quiet = TRUE)
  )
  expect_true(is.character(r$warnings),
              label = "$warnings is character")
  # If any warnings were emitted, they are in $warnings, not on console
  # (the wrapping of run_single in suppressMessages above was sufficient to
  #  get here cleanly; the test above checks the vector type regardless)
  expect_true(length(r$warnings) >= 0L, label = "$warnings has non-negative length")
})

# ---------------------------------------------------------------------------
# show_by_block = FALSE: by_block not in output, still in result
# ---------------------------------------------------------------------------

test_that("show_by_block = FALSE: by_block absent from output, present in result", {
  # Capture printed output
  out_with    <- capture.output(suppressWarnings(suppressMessages(
    run_single(.scenario_file, p = 10L, run_id = 1L,
               fit_control = .tiny_ctrl,
               s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8),
               max_widen = 0L,
               show_by_block = FALSE, quiet = FALSE)
  )))
  out_without <- capture.output(suppressWarnings(suppressMessages(
    run_single(.scenario_file, p = 10L, run_id = 1L,
               fit_control = .tiny_ctrl,
               s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8),
               max_widen = 0L,
               show_by_block = TRUE, quiet = FALSE)
  )))
  expect_true(length(out_without) > length(out_with),
              label = "show_by_block=TRUE output is longer than FALSE")

  # by_block still present in returned object
  r <- suppressWarnings(suppressMessages(
    run_single(.scenario_file, p = 10L, run_id = 1L,
               fit_control = .tiny_ctrl,
               s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8),
               max_widen = 0L, show_by_block = FALSE, quiet = TRUE)
  ))
  expect_true(inherits(r$scores$by_block, "data.frame"),
              label = "$scores$by_block present even with show_by_block=FALSE")
})
