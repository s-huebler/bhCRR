library(testthat)

# Resolve the repo root robustly whether called via Rscript --file=, via
# test_file(), or from an interactive session in the repo.
.resolve_root <- function() {
  # 1. REPO_ROOT env var (set by sim/chpc/config.sh on CHPC)
  e <- Sys.getenv("REPO_ROOT")
  if (nzchar(e)) return(normalizePath(e, mustWork = TRUE))
  # 2. --file= arg when invoked as Rscript --file=sim/tests/test_load.R;
  #    the arg may be a relative path, so anchor it to getwd() first
  args <- commandArgs(trailingOnly = FALSE)
  fa <- grep("--file=", args, value = TRUE)
  if (length(fa)) {
    raw <- sub("--file=", "", fa[1])
    if (!grepl("^/", raw)) raw <- file.path(getwd(), raw)
    p <- normalizePath(dirname(raw), mustWork = FALSE)
    # go up two levels: sim/tests -> sim -> repo root
    candidate <- normalizePath(file.path(p, "..", ".."), mustWork = FALSE)
    if (file.exists(file.path(candidate, "DESCRIPTION")))
      return(candidate)
  }
  # 3. Walk up from getwd() (works when run_all.R or test_file() is called
  #    from the repo root, which is the documented usage)
  path <- normalizePath(getwd(), mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "DESCRIPTION"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) break
    path <- parent
  }
  stop("Cannot locate repo root; set REPO_ROOT or run from inside the repo.")
}

root <- .resolve_root()
source(file.path(root, "sim", "R", "load.r"))

# ---------------------------------------------------------------------------
test_that("repo_root() returns a directory whose DESCRIPTION names bhCRR", {
  r <- repo_root()
  expect_true(file.exists(file.path(r, "DESCRIPTION")))
  lines <- readLines(file.path(r, "DESCRIPTION"), warn = FALSE)
  pkg_line <- grep("^Package:", lines, value = TRUE)
  expect_true(length(pkg_line) > 0)
  expect_equal(trimws(sub("^Package:\\s*", "", pkg_line[1])), "bhCRR")
})

# ---------------------------------------------------------------------------
test_that("REPO_ROOT pointing at a dir with no DESCRIPTION errors", {
  empty_dir <- tempfile()
  dir.create(empty_dir)
  on.exit(unlink(empty_dir, recursive = TRUE), add = TRUE)

  withr::with_envvar(c(REPO_ROOT = empty_dir), {
    expect_error(repo_root(), regexp = "DESCRIPTION")
  })
})

test_that("REPO_ROOT pointing at a different package's DESCRIPTION errors", {
  other_dir <- tempfile()
  dir.create(other_dir)
  on.exit(unlink(other_dir, recursive = TRUE), add = TRUE)
  writeLines(c("Package: NotBhCRR", "Version: 1.0.0"),
             file.path(other_dir, "DESCRIPTION"))

  withr::with_envvar(c(REPO_ROOT = other_dir), {
    expect_error(repo_root(), regexp = "Package: bhCRR")
  })
})

# ---------------------------------------------------------------------------
test_that("load_bhcrr('dev') attaches bhCRR and fit_ssl_psdh is callable", {
  # Force a fresh load so we test the actual attachment
  load_bhcrr("dev", force = TRUE)
  expect_true("package:bhCRR" %in% search())
  expect_true(existsFunction("fit_ssl_psdh"))
})

test_that("second consecutive load_bhcrr('dev') returns FALSE invisibly", {
  # Ensure it's already loaded from the previous test
  load_bhcrr("dev", force = FALSE)
  result <- withVisible(load_bhcrr("dev", force = FALSE))
  expect_false(result$value)
  expect_false(result$visible)
})

# ---------------------------------------------------------------------------
test_that("load_bhcrr('nonsense') errors naming both valid modes", {
  expect_error(load_bhcrr("nonsense"), regexp = "dev.*cluster|cluster.*dev")
})

# ---------------------------------------------------------------------------
test_that("load_bhcrr('cluster') succeeds when SO is present", {
  so_path <- file.path(repo_root(), "src",
                       paste0("bhCRR", .Platform$dynlib.ext))
  skip_if(!file.exists(so_path),
          paste("src/bhCRR.so not present; skipping cluster-mode test.",
                "Run: Rscript -e 'pkgbuild::compile_dll()' from repo root."))
  expect_true(load_bhcrr("cluster", force = TRUE))
  expect_true("package:bhCRR" %in% search())
})

# ---------------------------------------------------------------------------
test_that("check_compiled() errors with compile_dll remedy when SO missing", {
  fake_root <- tempfile()
  dir.create(file.path(fake_root, "src"), recursive = TRUE)
  on.exit(unlink(fake_root, recursive = TRUE), add = TRUE)
  expect_error(check_compiled(fake_root), regexp = "compile_dll")
})

# ---------------------------------------------------------------------------
test_that("export_all = FALSE: .cv_fold_path not exported, bhcrr_cv is", {
  load_bhcrr("dev", force = TRUE)
  pkg_env <- as.environment("package:bhCRR")

  expect_false(exists(".cv_fold_path", envir = pkg_env, inherits = FALSE))
  expect_true(exists("bhcrr_cv",       envir = pkg_env, inherits = FALSE))
})

# ---------------------------------------------------------------------------
test_that("load_sim() returns seven module paths and they all exist", {
  paths <- load_sim("dev", force = FALSE)
  expect_length(paths, 7L)
  for (p in paths) {
    expect_true(file.exists(p), info = paste("missing:", p))
  }
})

existsFunction <- function(name) {
  tryCatch({
    f <- get(name, envir = as.environment("package:bhCRR"), inherits = FALSE)
    is.function(f)
  }, error = function(e) FALSE)
}
