library(testthat)
library(withr)

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
# A2: batch_sweep() runs_root and manifest_dir overrides
# ---------------------------------------------------------------------------

test_that("runs_root override: files land directly under runs_root/p<P>/", {
  tmp <- file.path(tempdir(), paste0("chpc_runs_", Sys.getpid()))
  dir.create(tmp, showWarnings = FALSE)

  suppressWarnings(
    batch_sweep(.scenario_file,
                npredictors  = 10L,
                run_start    = 1L, run_end = 1L,
                results_root = tempdir(),  # ignored when runs_root set
                runs_root    = tmp,
                fit_control  = .tiny_ctrl,
                s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8),
                max_widen = 0L)
  )

  expected <- file.path(tmp, "p10", "p10_run1.test_run1.Rdata")
  expect_true(file.exists(expected), label = "run file under runs_root/p<P>/")

  # No extra "runs" level introduced
  wrong <- file.path(tmp, "runs", "p10", "p10_run1.test_run1.Rdata")
  expect_false(file.exists(wrong), label = "no runs/ sub-level when runs_root given")
})

test_that("manifest_dir override: scenario.env written to manifest_dir, not elsewhere", {
  tmp_runs  <- file.path(tempdir(), paste0("chpc_runs2_", Sys.getpid()))
  tmp_mfst  <- file.path(tempdir(), paste0("chpc_mfst_", Sys.getpid()))
  dir.create(tmp_runs, showWarnings = FALSE)
  dir.create(tmp_mfst, showWarnings = FALSE)

  suppressWarnings(
    batch_sweep(.scenario_file,
                npredictors  = 10L,
                run_start    = 1L, run_end = 1L,
                runs_root    = tmp_runs,
                manifest_dir = tmp_mfst,
                fit_control  = .tiny_ctrl,
                s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8),
                max_widen = 0L)
  )

  expect_true(file.exists(file.path(tmp_mfst, "scenario.env")),
              label = "scenario.env in manifest_dir")
  expect_false(file.exists(file.path(tmp_runs, "scenario.env")),
               label = "scenario.env NOT in runs_root")
})

test_that("default layout (both NULL) unchanged: runs under results_root/<scen>/runs/p<P>/", {
  tmp_root <- file.path(tempdir(), paste0("chpc_default_", Sys.getpid()))
  dir.create(tmp_root, showWarnings = FALSE)

  suppressWarnings(
    batch_sweep(.scenario_file,
                npredictors  = 10L,
                run_start    = 1L, run_end = 1L,
                results_root = tmp_root,
                fit_control  = .tiny_ctrl,
                s0_seq = c(0.04, 0.08), s1_seq = c(0.4, 0.8),
                max_widen = 0L)
  )

  expected <- file.path(tmp_root, "run1.test", "runs", "p10",
                         "p10_run1.test_run1.Rdata")
  expect_true(file.exists(expected), label = "default layout unchanged")
  expect_true(file.exists(file.path(tmp_root, "run1.test", "scenario.env")),
              label = "manifest in scenario_dir by default")
})

# ---------------------------------------------------------------------------
# Helpers for running CHPC wrapper scripts as subprocesses
# ---------------------------------------------------------------------------

.run_wrapper <- function(script, envvars) {
  # Run `Rscript script` with the given environment variable overrides.
  # Returns list(status, output).
  old <- lapply(names(envvars), Sys.getenv)
  names(old) <- names(envvars)
  do.call(Sys.setenv, as.list(envvars))
  on.exit(do.call(Sys.setenv, as.list(old)), add = TRUE)

  out    <- system2("Rscript", script, stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  list(status = if (is.null(status)) 0L else status,
       output = paste(out, collapse = "\n"))
}

.run_array_script <- file.path(.test_root, "sim", "chpc", "R", "run_array.r")
.parse_script     <- file.path(.test_root, "sim", "chpc", "R", "parse.r")

# ---------------------------------------------------------------------------
# run_array.r: missing required env vars -> named error
# ---------------------------------------------------------------------------

.base_env <- c(
  REPO_ROOT       = .test_root,
  SCENARIO        = "run1.test",
  SCRATCH_RUN_DIR = tempdir(),
  NPREDICTORS     = "10",
  RUN_START       = "1",
  RUN_END         = "1"
)

test_that("run_array.r: missing REPO_ROOT errors naming the variable", {
  env <- .base_env[setdiff(names(.base_env), "REPO_ROOT")]
  env["REPO_ROOT"] <- ""   # unset by setting empty
  r <- .run_wrapper(.run_array_script, env)
  expect_false(r$status == 0L, label = "non-zero exit on missing REPO_ROOT")
  expect_true(grepl("REPO_ROOT", r$output, ignore.case = FALSE),
              label = "error names REPO_ROOT")
})

test_that("run_array.r: missing SCENARIO errors naming the variable", {
  env <- .base_env
  env["SCENARIO"] <- ""
  r <- .run_wrapper(.run_array_script, env)
  expect_false(r$status == 0L)
  expect_true(grepl("SCENARIO", r$output), label = "error names SCENARIO")
})

test_that("run_array.r: missing SCRATCH_RUN_DIR errors naming the variable", {
  env <- .base_env
  env["SCRATCH_RUN_DIR"] <- ""
  r <- .run_wrapper(.run_array_script, env)
  expect_false(r$status == 0L)
  expect_true(grepl("SCRATCH_RUN_DIR", r$output),
              label = "error names SCRATCH_RUN_DIR")
})

test_that("run_array.r: missing NPREDICTORS errors naming the variable", {
  env <- .base_env
  env["NPREDICTORS"] <- ""
  r <- .run_wrapper(.run_array_script, env)
  expect_false(r$status == 0L)
  expect_true(grepl("NPREDICTORS", r$output),
              label = "error names NPREDICTORS")
})

# ---------------------------------------------------------------------------
# run_array.r: DGP env vars set -> named error
# ---------------------------------------------------------------------------

test_that("run_array.r: any DGP env var set -> error naming the variable(s)", {
  env <- .base_env
  env["NOBS"] <- "200"   # a DGP variable; must trigger the guard
  r <- .run_wrapper(.run_array_script, env)
  expect_false(r$status == 0L, label = "non-zero exit when DGP var set")
  expect_true(grepl("NOBS", r$output),
              label = "error message names the offending DGP var")
})

test_that("run_array.r: multiple DGP env vars set -> all named in error", {
  env <- .base_env
  env["NOBS"]      <- "200"
  env["CENS_RATE"] <- "0.05"
  r <- .run_wrapper(.run_array_script, env)
  expect_false(r$status == 0L)
  expect_true(grepl("NOBS",      r$output), label = "NOBS named")
  expect_true(grepl("CENS_RATE", r$output), label = "CENS_RATE named")
})

# ---------------------------------------------------------------------------
# run_array.r: missing scenario file -> informative error
# ---------------------------------------------------------------------------

test_that("run_array.r: nonexistent scenario file -> error", {
  env <- .base_env
  env["SCENARIO"] <- "nonexistent.scenario"
  r <- .run_wrapper(.run_array_script, env)
  expect_false(r$status == 0L)
  expect_true(grepl("nonexistent.scenario", r$output, ignore.case = FALSE) ||
                grepl("[Nn]ot found|[Ss]cenario", r$output),
              label = "error mentions the missing scenario")
})

# ---------------------------------------------------------------------------
# End-to-end via run_array.r (cluster mode; requires pre-built SO)
# ---------------------------------------------------------------------------

test_that("run_array.r end-to-end: canonical files appear under scratch/p<P>/", {
  so_path <- file.path(.test_root, "src",
                        paste0("bhCRR", .Platform$dynlib.ext))
  skip_if_not(file.exists(so_path),
              message = "Compiled SO absent; skipping cluster-mode integration test")

  tmp_scratch <- file.path(tempdir(), paste0("chpc_e2e_", Sys.getpid()))
  dir.create(tmp_scratch, showWarnings = FALSE)

  env <- c(
    REPO_ROOT       = .test_root,
    SCENARIO        = "run1.test",
    SCRATCH_RUN_DIR = tmp_scratch,
    NPREDICTORS     = "10",
    RUN_START       = "1",
    RUN_END         = "1"
  )

  r <- .run_wrapper(.run_array_script, env)
  expect_equal(r$status, 0L, label = paste("run_array.r exited 0\n", r$output))

  expected <- file.path(tmp_scratch, "p10", "p10_run1.test_run1.Rdata")
  expect_true(file.exists(expected),
              label = "canonical run file under scratch/p10/")

  fn <- basename(expected)
  rx <- "^p(\\d+)_([A-Za-z0-9][A-Za-z0-9.-]*)_run(\\d+)\\.Rdata$"
  expect_true(grepl(rx, fn),
              label = "filename matches parse_results() regex")
})

# ---------------------------------------------------------------------------
# parse.r (CHPC): missing env vars -> named error
# ---------------------------------------------------------------------------

test_that("parse.r (CHPC): missing REPO_ROOT errors naming the variable", {
  env <- c(REPO_ROOT = "", SCRATCH_RUN_DIR = tempdir(), SCENARIO = "run1.test")
  r <- .run_wrapper(.parse_script, env)
  expect_false(r$status == 0L)
  expect_true(grepl("REPO_ROOT", r$output), label = "error names REPO_ROOT")
})

test_that("parse.r (CHPC): missing SCRATCH_RUN_DIR errors naming the variable", {
  env <- c(REPO_ROOT = .test_root, SCRATCH_RUN_DIR = "", SCENARIO = "run1.test")
  r <- .run_wrapper(.parse_script, env)
  expect_false(r$status == 0L)
  expect_true(grepl("SCRATCH_RUN_DIR", r$output),
              label = "error names SCRATCH_RUN_DIR")
})

test_that("parse.r (CHPC): missing SCENARIO errors naming the variable", {
  env <- c(REPO_ROOT = .test_root, SCRATCH_RUN_DIR = tempdir(), SCENARIO = "")
  r <- .run_wrapper(.parse_script, env)
  expect_false(r$status == 0L)
  expect_true(grepl("SCENARIO", r$output), label = "error names SCENARIO")
})

# ---------------------------------------------------------------------------
# parse.r (CHPC): end-to-end parse of run_array.r output
# ---------------------------------------------------------------------------

test_that("parse.r (CHPC) end-to-end: writes CSVs into out_dir", {
  so_path <- file.path(.test_root, "src",
                        paste0("bhCRR", .Platform$dynlib.ext))
  skip_if_not(file.exists(so_path),
              message = "Compiled SO absent; skipping cluster-mode parse integration test")

  tmp_scratch <- file.path(tempdir(), paste0("chpc_parse_e2e_", Sys.getpid()))
  dir.create(tmp_scratch, showWarnings = FALSE)

  # 1. Run two runs via run_array.r
  run_env <- c(
    REPO_ROOT       = .test_root,
    SCENARIO        = "run1.test",
    SCRATCH_RUN_DIR = tmp_scratch,
    NPREDICTORS     = "10",
    RUN_START       = "1",
    RUN_END         = "2"
  )
  r_run <- .run_wrapper(.run_array_script, run_env)
  skip_if_not(r_run$status == 0L,
              message = paste("run_array.r failed:", r_run$output))

  # 2. Parse via parse.r (CHPC)
  parse_env <- c(
    REPO_ROOT       = .test_root,
    SCRATCH_RUN_DIR = tmp_scratch,
    SCENARIO        = "run1.test"
  )
  r_parse <- .run_wrapper(.parse_script, parse_env)
  expect_equal(r_parse$status, 0L,
               label = paste("parse.r exited 0\n", r_parse$output))

  # 3. Confirm output in sim/results/run1.test/parsed/
  out_dir <- file.path(.test_root, "sim", "results", "run1.test", "parsed")
  expect_true(file.exists(file.path(out_dir, "betas_p10_run1.test.csv")),
              label = "betas CSV written by parse.r")
  expect_true(file.exists(file.path(out_dir, "model_performance_p10_run1.test.csv")),
              label = "model_performance CSV written by parse.r")
})
