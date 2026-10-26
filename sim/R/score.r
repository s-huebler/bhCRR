# sim/R/score.r -- compute all evaluation metrics for one run
#
# OWNS
#   score_run(fit_result, comparator_results, x, y, scenario, beta1, block_map)
#                       compute all metrics for bhCRR and each comparator;
#                       all metrics are computed on every run because they are
#                       cheap relative to fitting, eliminating the separate
#                       "wolbers", "wolbers_cv", etc. script variants from the
#                       legacy code; per-block metrics are computed here so
#                       parse.r is a pure stacker
#
#   .confusion_counts(true_beta, estimate, tol = 0)
#                       internal; port of parse_batch_run.R:80; NA-tolerant
#
#   .confusion_rates(tp, fp, tn, fn)
#                       internal; precision, recall, specificity with explicit
#                       NA_real_ on zero denominators (never NaN, never 0)
#
# INPUTS
#   fit_result          list    output of fit_bhcrr(); must carry $eval_time,
#                               $fit_seconds, $coef (named, length p), $fit
#   comparator_results  list    output of fit_comparators() PLUS an element
#                               named "fg_oracle" from fit_fg_oracle(); must
#                               have elements: fg_oracle, lasso, alasso, scad,
#                               mcp; each shaped per comparators.r contract
#   x                   matrix n x p  raw design matrix
#   y                   matrix n x 2  (time, status)
#   scenario            list    validated bhcrr_scenario (not used directly;
#                               reserved for future metrics)
#   beta1               numeric p  full cause-1 true coefficient vector; zeros
#                               for inactive predictors
#   block_map           data.frame  from sim_data()$block_map; one row per
#                               predictor; columns: name, block, block_rho,
#                               designed, active1, active2; MUST cover every
#                               predictor (a missing predictor is an error,
#                               not a silently dropped row -- dgp.r guarantees
#                               one row per predictor)
#
# RETURNS
#   list(overall = data.frame, by_block = data.frame), stringsAsFactors = FALSE
#
#   $overall  six rows, one per model, in this order:
#               bhcrr, fg_oracle, lasso, alasso, scad, mcp
#     model         character
#     wolbers_c     numeric   IPCW C-index at eval_time from fit_result$eval_time
#     eval_time     numeric   evaluation horizon used; same value in all six rows;
#                             equals fit_result$eval_time; present so downstream
#                             code can verify what tau was in effect without
#                             reconstructing it from the scenario
#     deviance      numeric   -2 * partial log-pseudo-likelihood at selected coef
#     tp, fp, tn, fn   integer  confusion counts; threshold: |coef| > 0
#     precision     numeric   tp / (tp+fp); NA_real_ when tp+fp == 0
#     recall        numeric   tp / (tp+fn); NA_real_ when tp+fn == 0
#     specificity   numeric   tn / (tn+fp); NA_real_ when tn+fp == 0
#     mse_coef      numeric   mean((coef_raw - beta1)^2) over all p predictors
#     fit_seconds   numeric   wall seconds for THIS model only; NOT including
#                             comparator_results$setup_seconds (shared across all
#                             four penalised comparators); run_one.r records
#                             setup_seconds separately; to compare one comparator
#                             against bhCRR's total fairly, add setup_seconds
#                             once: comparator_total = setup_seconds + fit_seconds
#
#   $by_block  6 * nblocks rows, model-major, blocks in block_map label order:
#     model, block, block_rho, tp, fp, tn, fn, precision, recall, mse_coef
#     wolbers_c, eval_time, deviance, fit_seconds are ABSENT: they are whole-model
#     quantities; a C-index over a block subset is meaningless
#
# TRUTH SET NOTE
#   The cause-1 truth set is which(beta1 != 0). It is NOT block_map$designed.
#   A designed active predictor may carry beta1 == 0 by construction in the
#   source papers; treating $designed as truth would inflate false negatives
#   for every model.
#
# ASYMMETRY NOTE
#   bhCRR risk: predict(fit_result$fit, newx=x, prediction_time=tau) -- a CIF
#   Comparator risk: each model element's $risk field -- a linear predictor
#   This is deliberate and harmless: wolbers_c() is rank-based, and the CIF is
#   a monotone transform of the linear predictor, so both give identical orderings
#
# MSE_COEF SCALE NOTE
#   mse_coef ALWAYS uses coef_raw, never coef_std. bhCRR is fitted on raw x;
#   comparators on standardized x_std. coef_std is on a different scale and
#   comparing it to beta1 would be comparing different quantities. coef_raw for
#   all six models is in the same units as beta1.
#
# DEVIANCE NOTE
#   deviance = -2 * partial log-pseudo-likelihood, IN-SAMPLE, NO complexity
#   penalty. A denser model and the oracle look better by construction; this is
#   NOT a predictive metric. Fine-Gray partial log-pseudo-likelihood depends
#   only on the linear predictor, so values are comparable across raw and
#   standardized scales. Sources per model:
#     bhCRR      fit_result$fit$final_model_object$logLik (fcrrp from final M-step)
#     comparators  recovered from BIC: deviance = bic_path[index] - log(n)*df_path[index]
#                  (because bic = -2*logLik + log(n)*df)
#     fg_oracle    meta$fit$logLik (fastCrr object, single value)
#
# TIMING NOTE
#   For a fair bhCRR-vs-comparator wall-time comparison, the comparator total
#   for one model is: comparator_results$setup_seconds + model$fit_seconds.
#   aLASSO's fit_seconds already includes its own initial estimator (fastCrr
#   or ridge), since those weights are part of its cost.
#
# DEPENDS ON
#   sim/R/load.r
#   bhCRR::wolbers_c, bhCRR::predict.ssl_psdh
#
# REPLACES
#   Scoring logic scattered across sim/_legacy/code/batch_run_complex.R,
#   batch_run_wolbers.R, batch_run_cv_wolbers.R, single_run.R (all variants)
#   sim/chpc/R/run_batch.R scoring section
#   sim/_legacy/code/parse_batch_run.R confusion_counts() (lines 80-101)

# ---------------------------------------------------------------------------
# .confusion_counts -- port of parse_batch_run.R:80
# ---------------------------------------------------------------------------

.confusion_counts <- function(true_beta, estimate, tol = 0) {
  # NA handling: a failed fit with NA coefficients yields all-zero counts,
  # not an error.
  keep <- !is.na(estimate) & !is.na(true_beta)
  sel  <- keep & abs(estimate) > tol
  act  <- keep & true_beta != 0
  c(TP = sum(sel & act),
    FP = sum(sel & !act),
    TN = sum(!sel & !act & keep),
    FN = sum(!sel & act))
}

# ---------------------------------------------------------------------------
# .confusion_rates -- derived rates with explicit zero-denominator handling
# ---------------------------------------------------------------------------

.confusion_rates <- function(tp, fp, tn, fn) {
  precision   <- if ((tp + fp) == 0L)  NA_real_ else tp / (tp + fp)
  recall      <- if ((tp + fn) == 0L)  NA_real_ else tp / (tp + fn)
  specificity <- if ((tn + fp) == 0L)  NA_real_ else tn / (tn + fp)
  list(precision = precision, recall = recall, specificity = specificity)
}

# ---------------------------------------------------------------------------
# score_run
# ---------------------------------------------------------------------------

score_run <- function(fit_result, comparator_results, x, y, scenario,
                      beta1, block_map) {

  n <- nrow(x)
  p <- ncol(x)

  # ---- tau: fixed from fit_result; never recomputed ----
  tau <- fit_result$eval_time
  if (is.null(tau) || !is.finite(tau))
    stop("fit_result$eval_time is missing or non-finite. ",
         "fit_bhcrr() must record the evaluation horizon.")

  # ---- Verify block_map covers every predictor ----
  pred_names <- colnames(x)
  missing_preds <- setdiff(pred_names, block_map$name)
  if (length(missing_preds) > 0L)
    stop("block_map is missing entries for predictor(s): ",
         paste(missing_preds, collapse = ", "),
         ". dgp.r guarantees one row per predictor.")

  # ---- Build per-model description table ----
  # model label, coef_raw (length p), risk (length n), fit_seconds, converged
  model_names <- c("bhcrr", "fg_oracle", "lasso", "alasso", "scad", "mcp")

  bhcrr_coef_raw <- fit_result$coef   # already raw scale

  models <- list(
    bhcrr     = list(
      coef_raw    = bhcrr_coef_raw,
      converged   = fit_result$converged,
      fit_seconds = fit_result$fit_seconds
    ),
    fg_oracle = comparator_results$fg_oracle,
    lasso     = comparator_results$lasso,
    alasso    = comparator_results$alasso,
    scad      = comparator_results$scad,
    mcp       = comparator_results$mcp
  )

  # ---- Risk scores ----
  # bhCRR: CIF via predict.ssl_psdh; comparators: linear predictor $risk
  bhcrr_risk <- tryCatch(
    predict(fit_result$fit, newx = x, prediction_time = tau),
    error = function(e) {
      warning("predict.ssl_psdh failed: ", conditionMessage(e))
      rep(NA_real_, n)
    }
  )
  models$bhcrr$risk <- as.numeric(bhcrr_risk)

  # ---- Deviance ----
  # bhCRR: from final_model_object (fcrrp); single logLik value
  bhcrr_loglik <- tryCatch(
    as.numeric(fit_result$fit$final_model_object$logLik),
    error = function(e) NA_real_
  )

  .comparator_deviance <- function(m) {
    bp  <- m$meta$bic_path
    dp  <- m$meta$df_path
    idx <- m$meta$index
    if (is.null(bp) || is.null(dp) || is.null(idx) ||
        !is.finite(idx) || idx < 1L || idx > length(bp))
      return(NA_real_)
    # bic = -2*logLik + log(n)*df  =>  deviance = -2*logLik = bic - log(n)*df
    as.numeric(bp[[idx]] - log(n) * dp[[idx]])
  }

  oracle_deviance <- tryCatch({
    fit_obj <- comparator_results$fg_oracle$meta$fit
    if (is.null(fit_obj)) NA_real_ else -2 * as.numeric(fit_obj$logLik)
  }, error = function(e) NA_real_)

  deviance_vec <- c(
    bhcrr     = if (is.na(bhcrr_loglik)) NA_real_ else -2 * bhcrr_loglik,
    fg_oracle = oracle_deviance,
    lasso     = .comparator_deviance(comparator_results$lasso),
    alasso    = .comparator_deviance(comparator_results$alasso),
    scad      = .comparator_deviance(comparator_results$scad),
    mcp       = .comparator_deviance(comparator_results$mcp)
  )

  # ---- Assemble $overall ----
  overall_rows <- lapply(model_names, function(nm) {
    m    <- models[[nm]]
    coef <- m$coef_raw
    risk <- m$risk

    # wolbers_c
    wc <- tryCatch(
      wolbers_c(y_true = y, risk_score = risk, evaluation_time = tau),
      error = function(e) NA_real_
    )

    # confusion counts (uses cause-1 truth, NOT block_map$designed)
    cc <- .confusion_counts(true_beta = beta1, estimate = coef)
    rates <- .confusion_rates(cc["TP"], cc["FP"], cc["TN"], cc["FN"])

    # mse_coef: always coef_raw vs beta1
    mse <- if (any(is.na(coef))) NA_real_
           else mean((coef - beta1)^2)

    data.frame(
      model       = nm,
      wolbers_c   = as.numeric(wc),
      eval_time   = tau,
      deviance    = deviance_vec[[nm]],
      tp          = as.integer(cc["TP"]),
      fp          = as.integer(cc["FP"]),
      tn          = as.integer(cc["TN"]),
      fn          = as.integer(cc["FN"]),
      precision   = rates$precision,
      recall      = rates$recall,
      specificity = rates$specificity,
      mse_coef    = mse,
      fit_seconds = as.numeric(m$fit_seconds),
      stringsAsFactors = FALSE
    )
  })

  overall <- do.call(rbind, overall_rows)
  rownames(overall) <- NULL

  # ---- Assemble $by_block ----
  # Join block_map to predictor names; ordered by block_map's label order
  bm_ordered  <- block_map[match(pred_names, block_map$name), ]
  block_labels <- unique(bm_ordered$block)   # in block_map label order

  by_block_rows <- lapply(model_names, function(nm) {
    m    <- models[[nm]]
    coef <- m$coef_raw

    lapply(block_labels, function(b) {
      idx <- which(bm_ordered$block == b)
      if (length(idx) == 0L) return(NULL)

      b_rho <- bm_ordered$block_rho[idx[1L]]

      cc    <- .confusion_counts(beta1[idx], coef[idx])
      rates <- .confusion_rates(cc["TP"], cc["FP"], cc["TN"], cc["FN"])
      mse   <- if (any(is.na(coef[idx]))) NA_real_
               else mean((coef[idx] - beta1[idx])^2)

      data.frame(
        model       = nm,
        block       = b,
        block_rho   = b_rho,
        tp          = as.integer(cc["TP"]),
        fp          = as.integer(cc["FP"]),
        tn          = as.integer(cc["TN"]),
        fn          = as.integer(cc["FN"]),
        precision   = rates$precision,
        recall      = rates$recall,
        mse_coef    = mse,
        stringsAsFactors = FALSE
      )
    })
  })

  by_block <- do.call(rbind, c(unlist(by_block_rows, recursive = FALSE),
                                list(stringsAsFactors = FALSE)))
  rownames(by_block) <- NULL

  list(overall = overall, by_block = by_block)
}
