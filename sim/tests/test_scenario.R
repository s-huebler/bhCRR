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
load_sim("dev", force = FALSE)   # sources scenario.r (and stubs)

# ---------------------------------------------------------------------------
# A minimal valid spec used in multiple tests below
# ---------------------------------------------------------------------------
.valid_args <- list(
  name            = "test.scenario",
  nobs            = 100L,
  zero_gap_target = 0.1,
  beta1_active    = c(0.5, -0.5, 0.5),
  beta2_active    = "neg_beta1",
  active_block    = c("strong", "weak", "indep"),
  block_props     = c(indep = 0.60, weak = 0.20, strong = 0.20),
  block_rho       = c(indep = 0.00, weak = 0.30, strong = 0.70),
  latent_type     = "gaussian",
  latent_q        = c(indep = 0.50, weak = 0.40, strong = 0.50),
  cens_rate       = 0.05,
  u_min           = 50,
  u_max           = 150
)

.make_spec <- function(...) {
  args <- .valid_args
  extras <- list(...)
  for (nm in names(extras)) args[[nm]] <- extras[[nm]]
  do.call(scenario, args)
}

# ---------------------------------------------------------------------------
test_that("valid spec built by scenario() passes and has class bhcrr_scenario", {
  spec <- do.call(scenario, .valid_args)
  expect_s3_class(spec, "bhcrr_scenario")
  expect_equal(spec$name, "test.scenario")
  expect_equal(spec$nobs, 100L)
  expect_true(is.list(spec$blocks))
  expect_true(all(c("indep", "weak", "strong") %in% names(spec$blocks)))
})

# ---------------------------------------------------------------------------
# Validation rules
# ---------------------------------------------------------------------------

test_that("invalid name (underscore) is rejected with filename-format message", {
  bad <- .valid_args; bad$name <- "bad_name"
  expect_error(do.call(scenario, bad),
               regexp = "underscore|filename|_.*format|p<P>",
               ignore.case = TRUE)
})

test_that("invalid name (empty string) is rejected", {
  bad <- .valid_args; bad$name <- ""
  expect_error(do.call(scenario, bad), regexp = "name")
})

test_that("mismatched active_block length is rejected", {
  bad <- .valid_args; bad$active_block <- c("strong", "weak")  # length 2, not 3
  expect_error(do.call(scenario, bad), regexp = "active_block")
})

test_that("mismatched beta2_active (vector) length is rejected", {
  bad <- .valid_args
  bad$beta2_active <- c(-0.5, 0.5)  # length 2, not 3
  expect_error(do.call(scenario, bad), regexp = "beta2_active")
})

test_that("invalid beta2_active string is rejected", {
  bad <- .valid_args; bad$beta2_active <- "flip_beta1"
  expect_error(do.call(scenario, bad), regexp = "neg_beta1")
})

test_that("block_props without names is rejected", {
  bad <- .valid_args; bad$block_props <- c(0.6, 0.2, 0.2)
  expect_error(do.call(scenario, bad), regexp = "block_props")
})

test_that("block_rho with different label set is rejected", {
  bad <- .valid_args
  bad$block_rho <- c(indep = 0.0, other = 0.3, strong = 0.7)
  expect_error(do.call(scenario, bad), regexp = "block_rho")
})

test_that("latent_q with different label set is rejected", {
  bad <- .valid_args
  bad$latent_q <- c(indep = 0.5, other = 0.4, strong = 0.5)
  expect_error(do.call(scenario, bad), regexp = "latent_q")
})

test_that("active_block entry not in label set is rejected", {
  bad <- .valid_args
  bad$active_block <- c("strong", "weak", "mystery")
  expect_error(do.call(scenario, bad), regexp = "mystery|active_block")
})

test_that("block_props not summing to 1 is rejected", {
  bad <- .valid_args; bad$block_props <- c(indep = 0.50, weak = 0.30, strong = 0.30)
  expect_error(do.call(scenario, bad), regexp = "block_props")
})

test_that("block_rho >= 1 is rejected", {
  bad <- .valid_args; bad$block_rho <- c(indep = 0.00, weak = 0.30, strong = 1.00)
  expect_error(do.call(scenario, bad), regexp = "block_rho")
})

test_that("latent_q outside (0,1) is rejected", {
  bad <- .valid_args; bad$latent_q <- c(indep = 0.5, weak = 0.0, strong = 0.5)
  expect_error(do.call(scenario, bad), regexp = "latent_q")
})

test_that("negative cens_rate is rejected", {
  bad <- .valid_args; bad$cens_rate <- -0.1
  expect_error(do.call(scenario, bad), regexp = "cens_rate")
})

test_that("u_min > u_max is rejected", {
  bad <- .valid_args; bad$u_min <- 200; bad$u_max <- 100
  expect_error(do.call(scenario, bad), regexp = "u_min|u_max")
})

test_that("nobs < 1 is rejected", {
  bad <- .valid_args; bad$nobs <- 0
  expect_error(do.call(scenario, bad), regexp = "nobs")
})

test_that("invalid latent_type is rejected", {
  bad <- .valid_args; bad$latent_type <- "uniform"
  expect_error(do.call(scenario, bad), regexp = "latent_type")
})

# ---------------------------------------------------------------------------
# neg_beta1 sentinel
# ---------------------------------------------------------------------------

test_that("neg_beta1 resolves to -beta1_active and supplied form is recorded", {
  spec <- do.call(scenario, .valid_args)  # beta2_active = "neg_beta1"
  expect_equal(spec$beta2_active, -spec$beta1_active)
  expect_identical(spec$beta2_active_supplied, "neg_beta1")
})

test_that("explicit numeric beta2_active is stored and supplied form recorded", {
  args <- .valid_args
  args$beta2_active <- c(-0.4, 0.4, -0.4)
  spec <- do.call(scenario, args)
  expect_equal(spec$beta2_active, c(-0.4, 0.4, -0.4))
  expect_equal(spec$beta2_active_supplied, c(-0.4, 0.4, -0.4))
})

# ---------------------------------------------------------------------------
# scenario_from_file()
# ---------------------------------------------------------------------------

test_that("scenario_from_file() on _template.R is refused (underscore prefix)", {
  tmpl <- file.path(.test_root, "sim", "scenarios", "_template.R")
  skip_if(!file.exists(tmpl), "_template.R not found")
  expect_error(scenario_from_file(tmpl), regexp = "_|template|underscore|new_scenario",
               ignore.case = TRUE)
})

test_that("new_scenario() creates a file that scenario_from_file() loads", {
  td <- tempfile(); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  path <- new_scenario("smoketest", td)
  expect_true(file.exists(path))
  expect_equal(basename(path), "smoketest.R")

  spec <- scenario_from_file(path)
  expect_s3_class(spec, "bhcrr_scenario")
  expect_equal(spec$name, "smoketest")
})

test_that("new_scenario() refuses to overwrite an existing file", {
  td <- tempfile(); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  new_scenario("smoketest", td)
  expect_error(new_scenario("smoketest", td), regexp = "exists|overwrite")
})

test_that("name-filename mismatch is rejected", {
  td <- tempfile(); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  # Write a file whose internal name doesn't match its filename
  f <- file.path(td, "scenarioX.R")
  writeLines('scenario <- list(name = "scenarioY")', f)
  expect_error(scenario_from_file(f), regexp = "scenarioX|scenarioY|mismatch")
})

test_that("scenario file with library() call is rejected", {
  td <- tempfile(); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  f <- file.path(td, "badlib.R")
  writeLines(c('library(stats)',
               'scenario <- list(name = "badlib")'), f)
  expect_error(scenario_from_file(f), regexp = "library|forbidden")
})

test_that("scenario file with source() call is rejected", {
  td <- tempfile(); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  f <- file.path(td, "badsrc.R")
  writeLines(c('source("something.R")',
               'scenario <- list(name = "badsrc")'), f)
  expect_error(scenario_from_file(f), regexp = "source|forbidden")
})

# ---------------------------------------------------------------------------
# Manifest round-trip
# ---------------------------------------------------------------------------

test_that("write_manifest -> read_manifest -> check_manifest round-trip passes", {
  td <- tempfile(); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  spec <- do.call(scenario, .valid_args)
  write_manifest(spec, td)

  expect_no_error(check_manifest(spec, td))
})

test_that("check_manifest trips on altered nobs", {
  td <- tempfile(); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  spec <- do.call(scenario, .valid_args)
  write_manifest(spec, td)

  env_path <- file.path(td, "scenario.env")
  lines <- readLines(env_path, warn = FALSE)
  lines <- sub("^NOBS=.*$", "NOBS=999", lines)
  writeLines(lines, env_path)

  expect_error(check_manifest(spec, td), regexp = "NOBS")
})

test_that("check_manifest does not trip on 0.40 vs 0.4 formatting difference", {
  td <- tempfile(); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  spec <- do.call(scenario, .valid_args)
  write_manifest(spec, td)

  # Reformat BETA1_ACTIVE to bash-style trailing zeros
  env_path <- file.path(td, "scenario.env")
  lines <- readLines(env_path, warn = FALSE)
  lines <- sub("^BETA1_ACTIVE=.*$",
               "BETA1_ACTIVE=0.50,-0.50,0.50", lines)
  writeLines(lines, env_path)

  expect_no_error(check_manifest(spec, td))
})

# ---------------------------------------------------------------------------
# Manifest key order must match _sig_keys in run_sim.sh
# ---------------------------------------------------------------------------

test_that("manifest key set and order match _sig_keys in run_sim.sh", {
  run_sim_path <- file.path(.test_root, "sim", "chpc", "run_sim.sh")
  skip_if(!file.exists(run_sim_path), "run_sim.sh not found")

  lines <- readLines(run_sim_path, warn = FALSE)
  start_idx <- grep("^_sig_keys=\\(", lines)
  expect_length(start_idx, 1)

  # Collect lines until closing )
  end_idx <- start_idx
  while (!grepl("\\)", lines[end_idx])) end_idx <- end_idx + 1
  sig_block <- paste(lines[start_idx:end_idx], collapse = " ")

  # Extract all-caps words (the key names; include digits and underscores)
  bash_keys <- regmatches(sig_block,
                           gregexpr("[A-Z][A-Z0-9_]*", sig_block))[[1]]
  # Remove the variable name itself (_sig_keys expands to SIG and KEYS fragments;
  # filter to plausible key names using a minimum length and no lone letters)
  bash_keys <- bash_keys[nchar(bash_keys) >= 3]

  # Write a manifest and read its data-line keys in order
  td <- tempfile(); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  spec <- do.call(scenario, .valid_args)
  write_manifest(spec, td)
  man <- read_manifest(file.path(td, "scenario.env"))
  r_keys <- names(man)

  expect_equal(r_keys, bash_keys)
})
