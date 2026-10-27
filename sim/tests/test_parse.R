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

# ---------------------------------------------------------------------------
# Shared fixtures
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
# Regex: parse.r must be identical to parse_batch.R line 147
# ---------------------------------------------------------------------------

.pb_rx <- local({
  pb <- readLines(file.path(.test_root, "sim", "chpc", "R", "parse_batch.R"))
  rx_line <- grep('^rx <- ', pb, value = TRUE)[1L]
  eval(parse(text = sub('^rx <- ', '', rx_line)))
})

test_that("parse.r .parse_rx is identical to parse_batch.R rx", {
  expect_identical(.parse_rx, .pb_rx,
                   label = "parse.r .parse_rx matches parse_batch.R:147 rx")
})

# ---------------------------------------------------------------------------
# Build 1-group fixture: p=10, 2 runs
# ---------------------------------------------------------------------------

.tmp1 <- file.path(tempdir(), paste0("parse_1p_", Sys.getpid()))
dir.create(.tmp1, showWarnings = FALSE)

suppressWarnings(
  batch_sweep(.scenario_file,
              npredictors  = 10L,
              run_start    = 1L, run_end = 2L,
              results_root = .tmp1,
              fit_control  = .tiny_ctrl,
              s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8),
              max_widen = 0L)
)

.runs_root1 <- file.path(.tmp1, "run1.test", "runs")
.out1       <- file.path(.tmp1, "run1.test", "parsed")

.pr1 <- suppressMessages(
  parse_results(.runs_root1, "run1.test", out_dir = .out1)
)

# ---------------------------------------------------------------------------
# 1-group: file set
# ---------------------------------------------------------------------------

test_that("single p-group: expected output files are written", {
  expected <- c(
    "betas_p10_run1.test.csv",
    "betas_p10_run1.test.Rdata",
    "model_performance_p10_run1.test.csv",
    "by_block_p10_run1.test.csv",
    "false_positives_p10_run1.test.csv"
  )
  for (f in expected)
    expect_true(file.exists(file.path(.out1, f)), label = paste("exists:", f))
})

test_that("single p-group: no stacked model_performance file written", {
  expect_false(
    file.exists(file.path(.out1, "model_performance_run1.test.csv")),
    label = "stacked file absent with only one p-group")
})

test_that("return value n_runs == 2, n_failed == 0", {
  expect_equal(.pr1$n_runs,   2L)
  expect_equal(.pr1$n_failed, 0L)
})

# ---------------------------------------------------------------------------
# Betas table
# ---------------------------------------------------------------------------

test_that("betas row count == p * 6_models * n_runs", {
  betas <- read.csv(file.path(.out1, "betas_p10_run1.test.csv"))
  expect_equal(nrow(betas), 10L * 6L * 2L)
})

test_that("betas legacy columns present", {
  betas <- read.csv(file.path(.out1, "betas_p10_run1.test.csv"))
  required <- c("scenario", "predictor", "Block", "BlockRho",
                "True", "Model", "Estimate", "Bias", "nobs", "run_number")
  for (col in required)
    expect_true(col %in% names(betas), label = paste("column:", col))
})

test_that("betas has no FG column -- fg_oracle is a Model row", {
  betas <- read.csv(file.path(.out1, "betas_p10_run1.test.csv"))
  expect_false("FG" %in% names(betas))
})

test_that("betas contains fg_oracle as a Model value", {
  betas <- read.csv(file.path(.out1, "betas_p10_run1.test.csv"),
                    stringsAsFactors = FALSE)
  expect_true("fg_oracle" %in% betas$Model)
})

test_that("betas Bias == Estimate - True for every row", {
  betas <- read.csv(file.path(.out1, "betas_p10_run1.test.csv"))
  expect_equal(betas$Bias, betas$Estimate - betas$True, tolerance = 1e-12)
})

test_that("betas Block/BlockRho match block_map for spot-checked predictor", {
  betas <- read.csv(file.path(.out1, "betas_p10_run1.test.csv"),
                    stringsAsFactors = FALSE)
  run_file <- file.path(.runs_root1, "p10", "p10_run1.test_run1.Rdata")
  e <- new.env(parent = emptyenv())
  load(run_file, envir = e)
  bm   <- e$result$meta$block_map
  pred <- bm$name[1L]
  rows <- betas[betas$predictor == pred & betas$run_number == 1L, ]
  expect_equal(unique(rows$Block),    bm$block[bm$name == pred])
  expect_equal(unique(rows$BlockRho), bm$block_rho[bm$name == pred])
})

test_that("betas Estimate for comparators is coef_raw, not coef_std", {
  run_file <- file.path(.runs_root1, "p10", "p10_run1.test_run1.Rdata")
  e <- new.env(parent = emptyenv())
  load(run_file, envir = e)
  res       <- e$result
  # Find any comparator where coef_raw and coef_std differ (non-trivial test)
  found_diff <- FALSE
  for (mod in c("lasso", "alasso", "scad", "mcp")) {
    raw <- as.numeric(res$comparators[[mod]]$coef_raw)
    std <- as.numeric(res$comparators[[mod]]$coef_std)
    if (!isTRUE(all.equal(raw, std, tolerance = 1e-10))) {
      found_diff <- TRUE
      betas <- read.csv(file.path(.out1, "betas_p10_run1.test.csv"),
                        stringsAsFactors = FALSE)
      rows     <- betas[betas$Model == mod & betas$run_number == 1L, ]
      rows     <- rows[order(rows$predictor), ]
      raw_ord  <- raw[order(names(res$comparators[[mod]]$coef_raw))]
      expect_equal(rows$Estimate, raw_ord,
                   tolerance = 1e-10,
                   label = paste(mod, "Estimate == coef_raw, not coef_std"))
      break
    }
  }
  if (!found_diff)
    skip("coef_raw == coef_std for all comparators in run 1; cannot distinguish")
})

# ---------------------------------------------------------------------------
# model_performance
# ---------------------------------------------------------------------------

test_that("model_performance carries eval_time constant within each run", {
  mp <- read.csv(file.path(.out1, "model_performance_p10_run1.test.csv"),
                 stringsAsFactors = FALSE)
  expect_true("eval_time" %in% names(mp), label = "eval_time column present")
  for (rn in unique(mp$run_number)) {
    et <- unique(mp$eval_time[mp$run_number == rn])
    expect_equal(length(et), 1L,
                 label = paste("eval_time is unique within run", rn))
  }
})

# ---------------------------------------------------------------------------
# false_positives
# ---------------------------------------------------------------------------

test_that("false_positives fp_rate == sum(fp)/sum(fp+tn) verified independently", {
  fp <- read.csv(file.path(.out1, "false_positives_p10_run1.test.csv"),
                 stringsAsFactors = FALSE)
  bb <- read.csv(file.path(.out1, "by_block_p10_run1.test.csv"),
                 stringsAsFactors = FALSE)
  expect_true(nrow(fp) > 0L, label = "false_positives table is non-empty")
  row1   <- fp[1L, ]
  bb_sub <- bb[bb$Model == row1$Model & bb$Block == row1$Block &
                 bb$BlockRho == row1$BlockRho, ]
  expected_rate <- sum(bb_sub$fp) / sum(bb_sub$fp + bb_sub$tn)
  expect_equal(row1$fp_rate, expected_rate, tolerance = 1e-10,
               label = "fp_rate == sum(fp)/sum(fp+tn)")
})

# ---------------------------------------------------------------------------
# scenario.env
# ---------------------------------------------------------------------------

test_that("scenario.env is copied into out_dir", {
  expect_true(file.exists(file.path(.out1, "scenario.env")),
              label = "scenario.env present in out_dir")
})

# ---------------------------------------------------------------------------
# Two p-groups -> stacked file
# ---------------------------------------------------------------------------

.tmp2p <- file.path(tempdir(), paste0("parse_2p_", Sys.getpid()))
dir.create(.tmp2p, showWarnings = FALSE)

suppressWarnings({
  batch_sweep(.scenario_file,
              npredictors  = 10L,
              run_start    = 1L, run_end = 2L,
              results_root = .tmp2p,
              fit_control  = .tiny_ctrl,
              s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8),
              max_widen = 0L)
  batch_sweep(.scenario_file,
              npredictors  = 15L,
              run_start    = 1L, run_end = 2L,
              results_root = .tmp2p,
              fit_control  = .tiny_ctrl,
              s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8),
              max_widen = 0L)
})

.out2p <- file.path(.tmp2p, "run1.test", "parsed")
.pr2p  <- suppressMessages(
  parse_results(file.path(.tmp2p, "run1.test", "runs"),
                "run1.test",
                out_dir = .out2p)
)

test_that("two p-groups produce stacked model_performance file", {
  expect_true(
    file.exists(file.path(.out2p, "model_performance_run1.test.csv")),
    label = "stacked file written with two p-groups")
})

test_that("two p-groups: per-p betas and by_block files present for each p", {
  for (p in c(10L, 15L)) {
    expect_true(
      file.exists(file.path(.out2p, sprintf("betas_p%d_run1.test.csv", p))),
      label = sprintf("betas_p%d exists", p))
    expect_true(
      file.exists(file.path(.out2p, sprintf("by_block_p%d_run1.test.csv", p))),
      label = sprintf("by_block_p%d exists", p))
  }
})

# ---------------------------------------------------------------------------
# HARD STOP: scenario token mismatch in filename
# ---------------------------------------------------------------------------

test_that("scenario token mismatch in filename -> stop naming the offending file", {
  tmp_bs <- file.path(tempdir(), paste0("parse_scen_", Sys.getpid()))
  pdir   <- file.path(tmp_bs, "p10")
  dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
  f1 <- file.path(.runs_root1, "p10", "p10_run1.test_run1.Rdata")
  # Copy under a filename whose scenario token is 'other.scenario', not 'run1.test'
  file.copy(f1, file.path(pdir, "p10_other.scenario_run1.Rdata"))
  expect_error(
    suppressMessages(
      parse_results(tmp_bs, "run1.test", out_dir = tempdir())),
    regexp = "scenario|mismatch|token",
    ignore.case = TRUE,
    label = "scenario token mismatch stops parse"
  )
})

# ---------------------------------------------------------------------------
# HARD STOP: nobs mismatch across runs
# ---------------------------------------------------------------------------

test_that("nobs mismatch across runs -> stop naming both values", {
  tmp_nb <- file.path(tempdir(), paste0("parse_nobs_", Sys.getpid()))
  pdir   <- file.path(tmp_nb, "p10")
  dir.create(pdir, recursive = TRUE, showWarnings = FALSE)

  f1 <- file.path(.runs_root1, "p10", "p10_run1.test_run1.Rdata")
  f2 <- file.path(.runs_root1, "p10", "p10_run1.test_run2.Rdata")
  file.copy(f1, file.path(pdir, "p10_run1.test_run1.Rdata"))

  # Tamper with nobs in run 2
  e2 <- new.env(parent = emptyenv())
  load(f2, envir = e2)
  e2$result$meta$design_settings$nobs <- 999L
  e2$result$scenario$nobs             <- 999L
  save(list = "result",
       file  = file.path(pdir, "p10_run1.test_run2.Rdata"),
       envir = e2)

  expect_error(
    suppressMessages(
      parse_results(tmp_nb, "run1.test", out_dir = tempdir())),
    regexp = "nobs|mismatch",
    ignore.case = TRUE,
    label = "nobs mismatch stops parse"
  )
})

# ---------------------------------------------------------------------------
# HARD STOP: block_signature mismatch within p-group
# ---------------------------------------------------------------------------

test_that("block_signature mismatch within p-group -> stop printing both signatures", {
  tmp_bk <- file.path(tempdir(), paste0("parse_bsig_", Sys.getpid()))
  pdir   <- file.path(tmp_bk, "p10")
  dir.create(pdir, recursive = TRUE, showWarnings = FALSE)

  f1 <- file.path(.runs_root1, "p10", "p10_run1.test_run1.Rdata")
  f2 <- file.path(.runs_root1, "p10", "p10_run1.test_run2.Rdata")
  file.copy(f1, file.path(pdir, "p10_run1.test_run1.Rdata"))

  e2 <- new.env(parent = emptyenv())
  load(f2, envir = e2)
  e2$result$meta$block_signature <- "tampered_block_signature_xyz"
  save(list = "result",
       file  = file.path(pdir, "p10_run1.test_run2.Rdata"),
       envir = e2)

  expect_error(
    suppressMessages(
      parse_results(tmp_bk, "run1.test", out_dir = tempdir())),
    regexp = "block_signature|mismatch|signature",
    ignore.case = TRUE,
    label = "block_signature mismatch stops parse"
  )
})

# ---------------------------------------------------------------------------
# Failure records: skipped, counted, not a hard stop
# ---------------------------------------------------------------------------

test_that("failure record is skipped and counted; sweep continues and betas exclude it", {
  tmp_fl <- file.path(tempdir(), paste0("parse_fail_", Sys.getpid()))
  pdir   <- file.path(tmp_fl, "p10")
  dir.create(pdir, recursive = TRUE, showWarnings = FALSE)

  f1 <- file.path(.runs_root1, "p10", "p10_run1.test_run1.Rdata")
  file.copy(f1, file.path(pdir, "p10_run1.test_run1.Rdata"))

  # Write a synthetic failure record at run 2's canonical path
  result <- list(
    failed = TRUE,
    error  = "synthetic test failure",
    meta   = list(run_id = 2L, p = 10L, scenario_name = "run1.test",
                  seed = 999L, timestamp = Sys.time())
  )
  save(result, file = file.path(pdir, "p10_run1.test_run2.Rdata"))

  out_fl <- file.path(tmp_fl, "parsed")
  pr <- suppressMessages(
    parse_results(tmp_fl, "run1.test", out_dir = out_fl)
  )

  # No hard stop; exactly one succeeded, one skipped
  expect_equal(pr$n_runs,   1L, label = "n_runs == 1 (failure record not counted)")
  expect_equal(pr$n_failed, 1L, label = "n_failed == 1")

  # Error message preserved in summary
  expect_true(nrow(pr$failed_summary) == 1L)
  expect_true(grepl("synthetic test failure", pr$failed_summary$error[1L]),
              label = "failed_summary carries the error message")

  # Betas table has only 1 run * 6 models * 10 predictors rows
  betas <- read.csv(file.path(out_fl, "betas_p10_run1.test.csv"))
  expect_equal(nrow(betas), 10L * 6L * 1L,
               label = "betas row count excludes failure record")
})

# ---------------------------------------------------------------------------
# Non-matching filenames skipped; zero matching -> error
# ---------------------------------------------------------------------------

test_that("non-matching files are skipped and counted in n_skipped", {
  tmp_nm <- file.path(tempdir(), paste0("parse_nm_", Sys.getpid()))
  pdir   <- file.path(tmp_nm, "p10")
  dir.create(pdir, recursive = TRUE, showWarnings = FALSE)

  f1 <- file.path(.runs_root1, "p10", "p10_run1.test_run1.Rdata")
  file.copy(f1, file.path(pdir, "p10_run1.test_run1.Rdata"))
  # Old-format name that does not match the regex
  file.copy(f1, file.path(pdir, "old_n80_p10_run01.Rdata"))

  pr <- suppressMessages(
    parse_results(tmp_nm, "run1.test", out_dir = file.path(tmp_nm, "out"))
  )
  expect_equal(pr$n_skipped, 1L, label = "one non-matching file counted in n_skipped")
  expect_equal(pr$n_runs,    1L, label = "n_runs counts only canonical files")
})

test_that("no files matching canonical pattern -> stop", {
  tmp_all_bad <- file.path(tempdir(), paste0("parse_allbad_", Sys.getpid()))
  pdir        <- file.path(tmp_all_bad, "p10")
  dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
  f1 <- file.path(.runs_root1, "p10", "p10_run1.test_run1.Rdata")
  file.copy(f1, file.path(pdir, "old_n80_p10_run01.Rdata"))
  expect_error(
    suppressMessages(
      parse_results(tmp_all_bad, "run1.test", out_dir = tempdir())),
    regexp = "No files matching",
    label = "zero matching files causes an error"
  )
})

test_that("empty directory (no .Rdata files) -> stop", {
  tmp_empty <- file.path(tempdir(), paste0("parse_empty_", Sys.getpid()))
  dir.create(tmp_empty, showWarnings = FALSE)
  expect_error(
    suppressMessages(
      parse_results(tmp_empty, "run1.test", out_dir = tempdir())),
    regexp = "No .Rdata files",
    label = "empty directory causes an error"
  )
})
