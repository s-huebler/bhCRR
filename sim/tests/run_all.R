#!/usr/bin/env Rscript
# sim/tests/run_all.R -- run all sim/ tests in dependency order
#
# Usage (from repo root):
#   Rscript sim/tests/run_all.R
#
# These tests do NOT live in tests/testthat/ because sim/ is excluded from the
# package build by ^sim$ in .Rbuildignore; R CMD check would never find them.
# Each test file is also runnable standalone:
#   testthat::test_file("sim/tests/test_load.R")

library(testthat)

# Resolve repo root from this file's location regardless of working directory
this_file <- normalizePath(
  {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("--file=", args, value = TRUE)
    if (length(file_arg)) sub("--file=", "", file_arg[1]) else "sim/tests/run_all.R"
  },
  mustWork = FALSE
)

sim_tests_dir <- normalizePath(dirname(this_file), mustWork = FALSE)
if (!file.exists(sim_tests_dir))
  sim_tests_dir <- file.path(getwd(), "sim", "tests")

# Ordered list of test files (add new ones here in dependency order)
test_files <- c(
  "test_load.R",
  "test_scenario.R",
  "test_dgp.R",
  "test_fit.R"
)

all_passed <- TRUE
for (tf in test_files) {
  path <- file.path(sim_tests_dir, tf)
  cat("\n--- Running", tf, "---\n")
  result <- tryCatch(
    testthat::test_file(path, reporter = "progress"),
    error = function(e) {
      cat("ERROR running", tf, ":", conditionMessage(e), "\n")
      all_passed <<- FALSE
      NULL
    }
  )
  if (!is.null(result)) {
    failed <- sum(as.data.frame(result)$failed, na.rm = TRUE)
    if (failed > 0) all_passed <- FALSE
  }
}

if (!all_passed) {
  cat("\nSome tests FAILED.\n")
  quit(status = 1)
} else {
  cat("\nAll sim/ tests passed.\n")
}
