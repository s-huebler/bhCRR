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
# Configuration: n=100, p=10.
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
# Lambda paths: LASSO/SCAD/MCP share one ordinary path; aLASSO has its own
# ---------------------------------------------------------------------------

test_that("LASSO, SCAD and MCP share one lambda path; aLASSO's differs", {
  lp <- lapply(c(lasso = "lasso", alasso = "alasso", scad = "scad", mcp = "mcp"),
               function(nm) .comp[[nm]]$meta$lambda_path)
  expect_identical(lp$scad, lp$lasso)
  expect_identical(lp$mcp,  lp$lasso)
  expect_equal(length(lp$alasso), length(lp$lasso))
  expect_false(isTRUE(all.equal(lp$alasso, lp$lasso)))
})

test_that("aLASSO path starts at the weighted KKT bound", {
  sc <- .comp$standardization
  xs <- sweep(sweep(.x, 2L, sc$center, "-"), 2L, sc$scale, "/")
  ns <- .fg_null_score(xs, .y[, "time"], .y[, "status"])
  w  <- .comp$alasso$meta$adaptive_weights
  expect_equal(.comp$alasso$meta$lambda_path[1],
               (1 + 1e-6) * max(abs(ns$grad) / w) / nrow(.x), tolerance = 1e-10)
})

test_that("every path is all-zero at its first lambda and not at its second", {
  # Checks that each lambda.max sits exactly at the KKT boundary for the
  # penalty actually fitted (weighted for aLASSO, ordinary otherwise).
  for (nm in c("lasso","alasso","scad","mcp")) {
    dp <- unname(.comp[[nm]]$meta$df_path)
    expect_equal(dp[1], 0, label = paste(nm, "df at lambda.max"))
    expect_gte(dp[2], 1,   label = paste(nm, "df at second lambda"))
  }
})

# ---------------------------------------------------------------------------
# Path diagnostics retained
# ---------------------------------------------------------------------------

test_that("meta carries per-lambda solver diagnostics of full length", {
  fields <- c("lambda_path","logLik_path","converged_path","iter_path",
              "bic_path","df_path")
  for (nm in c("lasso","alasso","scad","mcp")) {
    m <- .comp[[nm]]$meta
    for (f in fields)
      expect_equal(length(m[[f]]), 100L, label = paste(nm, f, "length"))
    expect_true(is.logical(m$converged_path), label = paste(nm, "converged_path"))
    expect_equal(m$selected_df, sum(.comp[[nm]]$coef_std != 0),
                 label = paste(nm, "selected_df"))
    expect_true(is.logical(m$selected_at_last_valid),
                label = paste(nm, "selected_at_last_valid"))
    expect_true(is.logical(m$selection_truncated) &&
                  length(m$selection_truncated) == 1L,
                label = paste(nm, "selection_truncated"))
  }
})

# ---------------------------------------------------------------------------
# aLASSO initializer and weights
# ---------------------------------------------------------------------------

test_that("aLASSO uses the ridge initializer and bounded, median-1 weights", {
  m <- .comp$alasso$meta
  expect_identical(m$initial_estimator, "ridge")
  expect_true(is.finite(m$ridge_lambda) && m$ridge_lambda > 0)
  expect_equal(length(m$initial_beta_std), ncol(.x))
  w <- m$adaptive_weights
  expect_identical(names(w), colnames(.x))
  expect_true(all(w >= 0.1 & w <= 25))
  # Clipping is monotone and 1 lies inside the bounds, so the median survives.
  expect_equal(stats::median(w), 1, tolerance = 1e-12)
  # Larger first-stage coefficients get smaller penalties.
  ord <- order(abs(m$initial_beta_std))
  expect_true(all(diff(w[ord]) <= 1e-12))
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

test_that("fit_fg_oracle coef_raw and coef_std are length p with colnames(x)", {
  p   <- ncol(.x)
  res <- fit_fg_oracle(.x, .y, .beta1)
  expect_equal(length(res$coef_raw), p,
               label = "coef_raw length p")
  expect_equal(length(res$coef_std), p,
               label = "coef_std length p")
  expect_identical(names(res$coef_raw), colnames(.x),
                   label = "coef_raw names")
  expect_identical(names(res$coef_std), colnames(.x),
                   label = "coef_std names")
  expect_equal(length(res$risk), nrow(.x),
               label = "risk length n")
})

test_that("fit_fg_oracle inactive entries are zero; active nonzero when converged; meta$active_idx correct", {
  res          <- fit_fg_oracle(.x, .y, .beta1)
  active_idx   <- which(.beta1 != 0)
  inactive_idx <- which(.beta1 == 0)
  expect_equal(res$meta$active_idx, active_idx,
               label = "meta$active_idx matches which(beta1 != 0)")
  expect_true(all(res$coef_raw[inactive_idx] == 0),
              label = "inactive coef_raw are zero")
  if (isTRUE(res$converged))
    expect_true(all(res$coef_raw[active_idx] != 0),
                label = "active coef_raw nonzero when converged")
})

test_that("fit_fg_oracle returns well-formed empty result when beta1 all zero", {
  beta_null <- rep(0, ncol(.x))
  res       <- fit_fg_oracle(.x, .y, beta_null)
  expect_equal(length(res$coef_raw), ncol(.x),
               label = "coef_raw length p when no actives")
  expect_true(all(res$coef_raw == 0),
              label = "coef_raw all zero when no actives")
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
