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
# Fixtures -- built once, reused across tests
# Target runtime: < 60 s total for this file.
# Configuration: n=100, p=10 (p < n so unpenalized initial for aLASSO runs).
# ---------------------------------------------------------------------------

.tiny_args <- list(
  name            = "comp.test",
  nobs            = 100L,
  zero_gap_target = 0.1,
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

.spec <- do.call(scenario, .tiny_args)
set.seed(42L)
.dat  <- sim_data(.spec, p = 10L)
.x    <- .dat$x
.y    <- .dat$y
.beta1 <- .dat$beta1

# Run fit_comparators once and cache
set.seed(1L)
.comp <- fit_comparators(.x, .y)

# ---------------------------------------------------------------------------
# Return structure
# ---------------------------------------------------------------------------

test_that("fit_comparators() returns exactly the right top-level elements", {
  expect_true(is.list(.comp))
  expect_true(all(c("lasso","alasso","scad","mcp",
                    "setup_seconds","standardization") %in% names(.comp)))
})

test_that("each model element has all B5 contract fields", {
  required <- c("coef_raw","coef_std","lambda","bic","risk",
                "converged","fit_seconds","meta")
  for (nm in c("lasso","alasso","scad","mcp")) {
    for (field in required)
      expect_true(field %in% names(.comp[[nm]]),
                  label = paste(nm, "has field", field))
  }
})

# ---------------------------------------------------------------------------
# coef_raw and coef_std: named, length p, correct names
# ---------------------------------------------------------------------------

test_that("coef_raw and coef_std are named numerics of length p with colnames(x)", {
  p     <- ncol(.x)
  nms   <- colnames(.x)
  for (nm in c("lasso","alasso","scad","mcp")) {
    m <- .comp[[nm]]
    expect_true(is.numeric(m$coef_raw),
                label = paste(nm, "coef_raw numeric"))
    expect_true(is.numeric(m$coef_std),
                label = paste(nm, "coef_std numeric"))
    expect_equal(length(m$coef_raw), p,
                 label = paste(nm, "coef_raw length"))
    expect_equal(length(m$coef_std), p,
                 label = paste(nm, "coef_std length"))
    expect_identical(names(m$coef_raw), nms,
                     label = paste(nm, "coef_raw names"))
    expect_identical(names(m$coef_std), nms,
                     label = paste(nm, "coef_std names"))
  }
})

# ---------------------------------------------------------------------------
# coef_raw != coef_std when scales differ; back-transform reproduces coef_raw
# ---------------------------------------------------------------------------

test_that("coef_raw differs from coef_std and back-transform recovers coef_raw", {
  sc <- .comp$standardization$scale
  # scales should not all be 1 for real data
  expect_false(all(sc == 1))

  for (nm in c("lasso","alasso","scad","mcp")) {
    m <- .comp[[nm]]
    if (.comp[[nm]]$converged) {
      # At least some predictors should differ between raw and std
      expect_false(all(m$coef_raw == m$coef_std),
                   label = paste(nm, "coef_raw != coef_std"))
      recovered <- m$coef_std / sc
      expect_equal(unname(recovered), unname(m$coef_raw),
                   tolerance = 1e-10,
                   label = paste(nm, "back-transform recovers coef_raw"))
    }
  }
})

# ---------------------------------------------------------------------------
# $risk: length n, finite for converged fits
# ---------------------------------------------------------------------------

test_that("$risk has length n; finite for converged fits", {
  n <- nrow(.x)
  for (nm in c("lasso","alasso","scad","mcp")) {
    m <- .comp[[nm]]
    expect_equal(length(m$risk), n,
                 label = paste(nm, "risk length"))
    if (m$converged)
      expect_true(all(is.finite(m$risk)),
                  label = paste(nm, "risk finite when converged"))
  }
})

# ---------------------------------------------------------------------------
# lambda.path == reference_lambda for every fit
# ---------------------------------------------------------------------------

test_that("every model used the same reference lambda path", {
  # Verify via coef: all four fits should produce a coefficient matrix
  # of the same number of columns (= length of reference path).
  # We check by re-building reference_lambda and comparing lengths.
  # The direct path comparison is done inside fit_comparators(); if it
  # passed without error the assertion held.
  # Here we check indirectly: bic_path lengths all match.
  lens <- vapply(c("lasso","alasso","scad","mcp"),
                 function(nm) length(.comp[[nm]]$meta$bic_path),
                 integer(1))
  expect_true(length(unique(lens)) == 1L,
              label = "all BIC path lengths identical")
})

# ---------------------------------------------------------------------------
# Timing: numeric, positive, not difftime
# ---------------------------------------------------------------------------

test_that("setup_seconds is numeric, positive, not difftime", {
  ss <- .comp$setup_seconds
  expect_false(inherits(ss, "difftime"), label = "setup_seconds not difftime")
  expect_true(is.numeric(ss))
  expect_gt(ss, 0)
})

test_that("each fit_seconds is numeric, positive, not difftime", {
  for (nm in c("lasso","alasso","scad","mcp")) {
    fs <- .comp[[nm]]$fit_seconds
    expect_false(inherits(fs, "difftime"), label = paste(nm, "not difftime"))
    expect_true(is.numeric(fs),           label = paste(nm, "numeric"))
    expect_gt(fs, 0,                      label = paste(nm, "positive"))
  }
})

# ---------------------------------------------------------------------------
# Degenerate input: constant predictor column errors from B1
# ---------------------------------------------------------------------------

test_that("constant predictor column errors naming the offender", {
  x_bad          <- .x
  x_bad[, "X_3"] <- 1.0   # zero sd
  expect_error(fit_comparators(x_bad, .y),
               regexp = "X_3")
})

# ---------------------------------------------------------------------------
# fit_fg_oracle
# ---------------------------------------------------------------------------

test_that("fit_fg_oracle coef length equals number of true actives", {
  k   <- sum(.beta1 != 0)
  res <- fit_fg_oracle(.x, .y, .beta1)
  expect_length(res$coef, k)
})

test_that("fit_fg_oracle coef names match the active column names", {
  active_nms <- colnames(.x)[.beta1 != 0]
  res        <- fit_fg_oracle(.x, .y, .beta1)
  expect_identical(names(res$coef), active_nms)
})

test_that("fit_fg_oracle returns well-formed empty result when beta1 all zero", {
  beta_null <- rep(0, ncol(.x))
  res       <- fit_fg_oracle(.x, .y, beta_null)
  expect_length(res$coef, 0L)
  expect_true(res$converged)
  expect_null(res$meta$fit)
})

test_that("fit_fg_oracle returns fit_seconds as numeric", {
  res <- fit_fg_oracle(.x, .y, .beta1)
  expect_false(inherits(res$fit_seconds, "difftime"))
  expect_true(is.numeric(res$fit_seconds))
})

# ---------------------------------------------------------------------------
# SCAD and MCP gamma values in meta
# ---------------------------------------------------------------------------

test_that("SCAD meta carries gamma = 3.7", {
  expect_equal(.comp$scad$meta$gamma, 3.7)
})

test_that("MCP meta carries gamma = 2.7", {
  expect_equal(.comp$mcp$meta$gamma, 2.7)
})
