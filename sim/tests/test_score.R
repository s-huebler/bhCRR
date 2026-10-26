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
# Fixtures -- built once; intended total runtime for this file < 90 s.
# n=100, p=10, 2 active predictors; tiny CV grid for speed.
# ---------------------------------------------------------------------------

.sargs <- list(
  name            = "score.test",
  nobs            = 100L,
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

set.seed(77L)
.dat     <- sim_data(.spec, p = 10L)
.x       <- .dat$x
.y       <- .dat$y
.beta1   <- .dat$beta1
.block_map <- .dat$block_map

.tiny_ctrl <- bhcrr_cv_control(
  nfolds      = 3L, ncv = 1L,
  init_method = "LASSO_cv",
  init_args   = list(lambda_path = 10^seq(log10(0.1), log10(0.001), length = 5)),
  fit_args    = list(initial_sparsity = 0.05, theta_a = 1, theta_b = 1,
                     epsilon = 1e-04, inner_maxit_start = 1000, maxit = 20L),
  parallel    = FALSE,
  eval_time   = (.spec$u_min + .spec$u_max) / 2
)

set.seed(2L)
.fit_r  <- fit_bhcrr(.x, .y, .spec, run_id = 1L,
                     control = .tiny_ctrl,
                     s0_seq = c(0.04, 0.08),
                     s1_seq = c(0.4, 0.8),
                     max_widen = 0L)

set.seed(3L)
.comp_r <- fit_comparators(.x, .y)
.oracle_r <- fit_fg_oracle(.x, .y, .beta1)

# Assemble the five-element list score_run expects
.all_comp <- c(.comp_r[c("lasso","alasso","scad","mcp",
                          "setup_seconds","standardization")],
               list(fg_oracle = .oracle_r))

set.seed(4L)
.scores <- score_run(.fit_r, .all_comp, .x, .y, .spec, .beta1, .block_map)

# ---------------------------------------------------------------------------
# Basic structure
# ---------------------------------------------------------------------------

test_that("$overall has exactly 6 rows in documented model order", {
  ov <- .scores$overall
  expect_equal(nrow(ov), 6L)
  expect_equal(ov$model, c("bhcrr","fg_oracle","lasso","alasso","scad","mcp"))
})

test_that("$by_block has 6 * nblocks rows", {
  nb <- length(unique(.block_map$block))
  expect_equal(nrow(.scores$by_block), 6L * nb)
})

test_that("neither $overall nor $by_block contains factors", {
  for (nm in names(.scores$overall))
    expect_false(is.factor(.scores$overall[[nm]]),
                 label = paste("overall$", nm, "not factor"))
  for (nm in names(.scores$by_block))
    expect_false(is.factor(.scores$by_block[[nm]]),
                 label = paste("by_block$", nm, "not factor"))
})

# ---------------------------------------------------------------------------
# Column classes
# ---------------------------------------------------------------------------

test_that("$overall has the documented column classes", {
  ov <- .scores$overall
  expect_true(is.character(ov$model))
  for (col in c("wolbers_c","eval_time","deviance","precision","recall",
                "specificity","mse_coef","fit_seconds"))
    expect_true(is.numeric(ov[[col]]),
                label = paste("overall$", col, "numeric"))
  for (col in c("tp","fp","tn","fn"))
    expect_true(is.integer(ov[[col]]),
                label = paste("overall$", col, "integer"))
})

test_that("$by_block has the documented column classes", {
  bb <- .scores$by_block
  expect_true(is.character(bb$model))
  expect_true(is.character(bb$block))
  expect_true(is.numeric(bb$block_rho))
  for (col in c("tp","fp","tn","fn"))
    expect_true(is.integer(bb[[col]]),
                label = paste("by_block$", col, "integer"))
  for (col in c("precision","recall","mse_coef"))
    expect_true(is.numeric(bb[[col]]),
                label = paste("by_block$", col, "numeric"))
})

# ---------------------------------------------------------------------------
# wolbers_c
# ---------------------------------------------------------------------------

test_that("wolbers_c is in [0,1] or NA for every model", {
  wc <- .scores$overall$wolbers_c
  ok <- is.na(wc) | (wc >= 0 & wc <= 1)
  expect_true(all(ok))
})

test_that("wolbers_c is finite for converged models", {
  ov <- .scores$overall
  for (i in seq_len(nrow(ov))) {
    # bhcrr convergence
    conv <- if (ov$model[i] == "bhcrr") isTRUE(.fit_r$converged)
            else {
              el <- .all_comp[[ov$model[i]]]
              isTRUE(el$converged)
            }
    if (conv)
      expect_true(is.finite(ov$wolbers_c[i]),
                  label = paste(ov$model[i], "wolbers_c finite"))
  }
})

# ---------------------------------------------------------------------------
# fg_oracle: tp == k and fp == 0 by construction
# ---------------------------------------------------------------------------

test_that("fg_oracle tp equals number of true actives and fp is 0", {
  k   <- sum(.beta1 != 0)
  row <- .scores$overall[.scores$overall$model == "fg_oracle", ]
  # Oracle is fitted only on active columns; those coef are non-zero (unpenalised)
  # so tp must equal k if the oracle converged
  if (isTRUE(.oracle_r$converged)) {
    expect_equal(row$tp, as.integer(k),
                 label = "fg_oracle tp == k")
    expect_equal(row$fp, 0L,
                 label = "fg_oracle fp == 0")
  }
})

# ---------------------------------------------------------------------------
# precision / recall / specificity return NA_real_ on zero denominator
# ---------------------------------------------------------------------------

test_that(".confusion_rates returns NA_real_ (not NaN, not 0) on zero denominators", {
  # precision: tp+fp == 0
  r1 <- .confusion_rates(0L, 0L, 5L, 2L)
  expect_identical(r1$precision, NA_real_)
  expect_false(is.nan(r1$precision))

  # recall: tp+fn == 0 (no true positives, no false negatives)
  r2 <- .confusion_rates(0L, 3L, 5L, 0L)
  expect_identical(r2$recall, NA_real_)
  expect_false(is.nan(r2$recall))

  # specificity: tn+fp == 0
  r3 <- .confusion_rates(2L, 0L, 0L, 1L)
  expect_identical(r3$specificity, NA_real_)
  expect_false(is.nan(r3$specificity))
})

# ---------------------------------------------------------------------------
# confusion counts on all-NA coefficient vector are all zero; no error
# ---------------------------------------------------------------------------

test_that(".confusion_counts on all-NA estimate produces all-zero counts", {
  cc <- .confusion_counts(true_beta = c(0, 1, -1, 0),
                          estimate  = c(NA_real_, NA_real_, NA_real_, NA_real_))
  expect_equal(unname(cc), c(0L, 0L, 0L, 0L))
})

# ---------------------------------------------------------------------------
# mse_coef uses coef_raw, not coef_std
# ---------------------------------------------------------------------------

test_that("mse_coef matches coef_raw computation, not coef_std", {
  # Take the LASSO result and confirm mse is computed from coef_raw
  lasso_row  <- .scores$overall[.scores$overall$model == "lasso", ]
  lasso_elem <- .all_comp$lasso
  expected_mse <- mean((lasso_elem$coef_raw - .beta1)^2)
  expect_equal(lasso_row$mse_coef, expected_mse, tolerance = 1e-10)

  # Confirm it differs from what coef_std would give (scales should differ)
  mse_from_std <- mean((lasso_elem$coef_std - .beta1)^2)
  if (!all(lasso_elem$coef_std == lasso_elem$coef_raw))
    expect_false(isTRUE(all.equal(lasso_row$mse_coef, mse_from_std)))
})

# ---------------------------------------------------------------------------
# eval_time propagation
# ---------------------------------------------------------------------------

test_that("scoring uses fit_result$eval_time: removing it causes an error", {
  bad_fit <- .fit_r
  bad_fit$eval_time <- NULL
  expect_error(
    score_run(bad_fit, .all_comp, .x, .y, .spec, .beta1, .block_map),
    regexp = "eval_time"
  )
})

test_that("$overall carries eval_time equal to fit_result$eval_time in every row", {
  # A value-change test on wolbers_c cannot distinguish "scored at the
  # supplied tau" from "scored at a quantile fallback that happened to
  # differ"; the column-presence test below is the direct observable.
  # wolbers_c is rank-based and its value need not change when tau changes
  # (identical comparable-pair sets are possible at n=100).
  expect_true("eval_time" %in% names(.scores$overall))
  expect_identical(unique(.scores$overall$eval_time), .fit_r$eval_time)
})

# ---------------------------------------------------------------------------
# by_block counts sum to overall counts for every model
# ---------------------------------------------------------------------------

test_that("by_block tp/fp/tn/fn sum to $overall counts for every model", {
  for (nm in c("bhcrr","fg_oracle","lasso","alasso","scad","mcp")) {
    ov_row <- .scores$overall[.scores$overall$model == nm, ]
    bb_sub <- .scores$by_block[.scores$by_block$model == nm, ]

    expect_equal(sum(bb_sub$tp), ov_row$tp,
                 label = paste(nm, "tp sum"))
    expect_equal(sum(bb_sub$fp), ov_row$fp,
                 label = paste(nm, "fp sum"))
    expect_equal(sum(bb_sub$tn), ov_row$tn,
                 label = paste(nm, "tn sum"))
    expect_equal(sum(bb_sub$fn), ov_row$fn,
                 label = paste(nm, "fn sum"))
  }
})

# ---------------------------------------------------------------------------
# Missing block_map predictor causes an error
# ---------------------------------------------------------------------------

test_that("block_map missing a predictor causes an error", {
  bad_map <- .block_map[-1L, ]   # drop first predictor
  expect_error(
    score_run(.fit_r, .all_comp, .x, .y, .spec, .beta1, bad_map),
    regexp = colnames(.x)[1L]
  )
})

# ---------------------------------------------------------------------------
# wolbers_c, deviance, fit_seconds absent from $by_block
# ---------------------------------------------------------------------------

test_that("wolbers_c, eval_time, deviance, fit_seconds are absent from $by_block", {
  bb_names <- names(.scores$by_block)
  expect_false("wolbers_c"   %in% bb_names)
  expect_false("eval_time"   %in% bb_names)
  expect_false("deviance"    %in% bb_names)
  expect_false("fit_seconds" %in% bb_names)
})
