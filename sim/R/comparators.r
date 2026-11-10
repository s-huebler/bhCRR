# sim/R/comparators.r -- LASSO, aLASSO, SCAD, MCP, and oracle comparator fits
#
# OWNS  (public functions)
#   fit_comparators(x, y)
#                    standardise x once, build the null score and lambda
#                    paths, fit the four penalised Fine-Gray models (LASSO,
#                    aLASSO, SCAD, MCP), select by BIC, and return the
#                    assembled result list (B5). LASSO, SCAD and MCP share one
#                    ordinary path; aLASSO gets its own weight-adjusted path.
#                    This is the only public entry point for the four comparators.
#
#   fit_fg_oracle(x, y, beta1)
#                    unpenalised Fine-Gray fit restricted to the TRUE active
#                    columns (where beta1 != 0); an UPPER BOUND on achievable
#                    performance, not a fair competitor. Kept separate from
#                    fit_comparators() because it consumes ground truth.
#                    Migration step 5 (fit.r) established this handoff;
#                    run_one.r (step 9) is responsible for calling both and
#                    assembling the full comparator block.
#                    Source: batch_run_complex.R:173 --
#                      fastcmprsk::fastCrr(
#                        fastcmprsk::Crisk(time, status) ~
#                        x[, which(beta1 != 0)], variance = FALSE)
#                    If the true active set is empty, returns a well-formed
#                    result with a zero-length coefficient vector.
#
# INTERNALS  (module-level, prefixed with .)
#   .fg_null_score(x_std, ftime, status, failcode, cencode)
#                    null Fine-Gray score (summed over events) and diagonal
#                    null information / n, from one pass over the events
#
#   .fg_lambda_path(grad, n, penalty_factor, nlambda, lambda_min_ratio)
#                    clean geometric lambda grid from the KKT condition of the
#                    null model, max_j |grad_j| / (n * penalty_factor_j);
#                    penalty_factor = NULL gives the ordinary path. Needed
#                    because fastCrrp()'s built-in grid degenerates when there
#                    are no censored observations
#
#   .path_diag(fit) per-lambda lambda, logLik, converged, iter from a fit
#
#   .fit_fastcrrp(x_std, ftime, status, predictor_names, nlambda,
#                 lambda_min_ratio, algorithm_tolerance, max_iterations,
#                 penalty, penalty_factor, gamma, lambda)
#                    shared fastCrrp wrapper; handles the getBreslowJumps
#                    version guard and coefficient-matrix orientation
#
#   .select_by_bic(fit, method_label, predictor_names, nobs, x_scale)
#                    BIC selection over a fitted lambda path; returns
#                    beta_std, beta_raw, and the full BIC path for
#                    diagnostics
#
# INPUTS (fit_comparators)
#   x      matrix n x p   raw design matrix from sim_data(); predictors are
#                         standardised ONCE inside fit_comparators() -- see
#                         standardization note below
#   y      matrix n x 2   outcome (colnames time, status); status in {0,1,2}
#
# INPUTS (fit_fg_oracle)
#   x      matrix n x p   same raw design matrix
#   y      matrix n x 2   same outcome matrix
#   beta1  numeric p      full cause-1 coefficient vector; used only to
#                         identify the true active set (which(beta1 != 0))
#
# RETURNS (fit_comparators)
#   named list:
#     $lasso / $alasso / $scad / $mcp  -- each:
#       list(
#         coef_raw    numeric p  named; coefficients in ORIGINAL predictor
#                                units; use this for mse_coef against beta_true
#                                (bhCRR is fitted on raw x; comparing coef_std
#                                to beta_true would be comparing different scales)
#         coef_std    numeric p  named; coefficients in standardized units;
#                                stored for diagnostics
#         lambda      numeric(1) BIC-selected lambda
#         bic         numeric(1) BIC at the selected model
#         risk        numeric n  linear predictor x_std %*% coef_std; a monotone
#                                risk score sufficient for wolbers_c() which is
#                                rank-based; no baseline hazard reconstruction
#                                needed (this is a deliberate asymmetry with
#                                bhCRR, which uses predict.ssl_psdh() -- the two
#                                give identical rankings)
#         converged   logical(1)
#         fit_seconds numeric(1) wall seconds for THIS model ONLY; excludes
#                                setup_seconds (standardization + shared lambda
#                                path), which is shared -- folding it into each
#                                model would count it four times.
#                                aLASSO's fit_seconds DOES include its own
#                                ridge initializer and weighted path, since
#                                those weights are part
#                                of its cost. To compare one comparator against
#                                bhCRR's fit_seconds total, add setup_seconds
#                                once: comparator_total = setup_seconds +
#                                fit_seconds["alasso"] (for example)
#         meta        list(...)  all methods:
#                                  method, index, bic_path, df_path,
#                                  selected_df, selected_at_last_valid (TRUE
#                                  when BIC chose the smallest lambda with a
#                                  converged solution, i.e. a boundary),
#                                  selection_truncated (TRUE when an
#                                  unconverged lambda had lower BIC than the
#                                  chosen one),
#                                  lambda_path, logLik_path, converged_path,
#                                  iter_path (iter == 5000 means max.iter hit)
#                                SCAD / MCP: gamma
#                                aLASSO: initial_estimator ("ridge"),
#                                  ridge_lambda, ridge_converged, ridge_iter,
#                                  initial_beta_std, weight_eps,
#                                  adaptive_weights, weight_quantiles,
#                                  n_weight_lower, n_weight_upper
#       )
#       a failed fit produces coef_raw = coef_std = rep(NA_real_, p),
#       risk = rep(NA_real_, n), and converged = FALSE; its meta carries
#       error = TRUE plus whatever path or weight diagnostics were computed
#       before the failure
#     $setup_seconds   numeric(1)  standardization + shared lambda path
#     $standardization list(center = numeric p, scale = numeric p)
#
# RETURNS (fit_fg_oracle)
#   SAME element shape as fit_comparators() entries (score_run treats all six
#   models identically):
#   list(
#     coef_raw    numeric p  named with colnames(x); fitted values at active
#                            positions, 0 elsewhere. The oracle is fitted on
#                            RAW active columns, so coef_raw is the fitted
#                            quantity. coef_std = coef_raw * x_scale is the
#                            exact reparameterization: x %*% coef_raw and
#                            x_std %*% coef_std have identical rank ordering
#     coef_std    numeric p  named; coef_raw * apply(x, 2, sd)
#     risk        numeric n  linear predictor as.numeric(x %*% coef_raw)
#     converged   logical(1) or NA if fastCrr does not report convergence
#     fit_seconds numeric(1)
#     meta        list(active_idx, fit = fastCrr object or NULL,
#                      degenerate = TRUE when active set is empty)
#   )
#   Zero-actives case: coef_raw = coef_std = rep(0, p), risk = rep(0, n),
#   converged = TRUE (vacuously), meta$degenerate = TRUE, meta$fit = NULL.
#
# ASSEMBLY NOTE
#   run_one.r (step 9) is responsible for assembling the five-element list
#   that score_run() consumes:
#     c(fit_comparators(x, y),
#       list(fg_oracle = fit_fg_oracle(x, y, beta1)))
#
# FUTURE WORK
#   Lambda is currently selected by BIC over a geometric path.
#   A cross-validated selection exists at
#   sim/_legacy/code/other_comparators_cv_section.R (799 lines, Wolbers C-index
#   criterion) and is deliberately not ported; comparator tuning is an open
#   problem to revisit. The BIC path has at least two fragile points:
#     1. BIC penalises df = number of nonzero coefficients (model size), not
#        the effective-df trace; this overpenalises on correlated predictors.
#     2. The null-score lambda.max is computed from cause-1 events only; under
#        high censoring the grid may be too wide, leaving the BIC minimum at a
#        nearly-saturated model (df close to p).
#     3. fastCrrp's convergence test is RELATIVE per coefficient
#        (|b_new - b_old| / |b_old| < eps, utils.c checkConvergence) and eps
#        here is 1e-7, ten times stricter than the package default. A
#        coefficient that has just entered at a tiny magnitude must settle to
#        seven significant digits, so small-lambda LASSO/SCAD fits hit
#        max.iter, and BIC silently skips those lambdas.
#   When returning to comparator tuning, examine bic_path, df_path,
#   converged_path and selected_at_last_valid from stored results to check
#   whether the selected model is consistently near the endpoints of the grid.
#   BIC targets neither held-out prediction nor false-discovery control; if
#   support recovery is the endpoint, EBIC or stability selection fits better,
#   and if prediction is, nested CV does.
#
# STANDARDIZATION NOTE
#   Predictors are centred and scaled to unit variance ONCE before any
#   comparator is called. The same x_std goes to every comparator. This puts
#   the adaptive weights on the same coefficient scale as the penalised fits
#   and avoids method-to-method preprocessing differences (mirrors
#   batch_run_complex.R lines 223-232).
#   bhCRR (fit_bhcrr) is fitted on the raw x; comparators use x_std. Comparing
#   coef_raw (back-transformed from x_std) to bhCRR's coef is apples-to-apples.
#
# DEPENDS ON
#   fastcmprsk::fastCrrp     (LASSO, aLASSO, SCAD, MCP penalised fits)
#   fastcmprsk::fastCrr      (fg_oracle unpenalised fit on active subset)
#   survival::survfit        (Kaplan-Meier of censoring in fg_lambda_path)
#   survival::Surv
#
# REPLACES
#   sim/_legacy/code/other_comparators_cv_section.R  (parked; not ported)
#   sim/_legacy/code/batch_run_complex.R lines 217-589 (comparator block)
#   sim/_legacy/code/batch_run_comparison.R          (entire file)
#   sim/_legacy/code/single_run_comparison.R         (entire file)

# ---------------------------------------------------------------------------
# .fg_null_score -- score and diagonal information of the null Fine-Gray model
# ---------------------------------------------------------------------------
# fastCrrp()'s built-in lambda grid degenerates when there are no censored
# observations (the IPCW/score projection collapses to roundoff noise), so the
# grid is built here from the KKT condition of the null Fine-Gray model. Port
# of batch_run_complex.R line 447.
#
# Returns
#   grad  numeric p  null score summed over cause-1 events (NOT divided by n)
#   info  numeric p  diagonal of the null observed information divided by n,
#                    i.e. the per-observation curvature fastCrrp's coordinate
#                    descent sees at beta = 0. Used only to put the aLASSO
#                    first-stage ridge penalty on a scale-free footing.
# Both come from one pass over the events; computing them together costs one
# extra colSums per event.

.fg_null_score <- function(x_std, ftime, status,
                           failcode = 1L,
                           cencode  = 0L) {
  n <- nrow(x_std)
  p <- ncol(x_std)

  cen <- as.integer(status == cencode)
  if (any(cen == 1L)) {
    kmC  <- survival::survfit(survival::Surv(ftime, cen) ~ 1)
    Gfun <- stats::stepfun(kmC$time, c(1, kmC$surv))
  } else {
    Gfun <- function(t) rep(1, length(t))
  }

  ev <- which(status == failcode)
  if (length(ev) == 0L)
    stop("No cause-", failcode, " events; cannot construct a lambda path.")

  grad <- numeric(p)
  info <- numeric(p)
  for (i in ev) {
    ti          <- ftime[i]
    at_risk     <- ftime >= ti
    comp_before <- ftime < ti & status != failcode & status != cencode
    w           <- numeric(n)
    w[at_risk]  <- 1
    if (any(comp_before))
      w[comp_before] <- Gfun(ti) / Gfun(ftime[comp_before])
    S    <- sum(w)
    xbar <- colSums(w * x_std) / S
    grad <- grad + (x_std[i, ] - xbar)
    info <- info + (colSums(w * x_std^2) / S - xbar^2)
  }

  list(grad = grad, info = info / n)
}

# ---------------------------------------------------------------------------
# .fg_lambda_path -- geometric lambda grid from the KKT null-score condition
# ---------------------------------------------------------------------------
# fastCrrp's coordinate update thresholds coefficient j at
# lambda * penalty.factor[j] (densePen.c: l1 = lam[l] * m[j]), with no
# rescaling of penalty.factor. The all-zero solution is therefore optimal for
#   lambda >= max_j |grad_j| / (n * penalty_factor_j).
# With penalty_factor = NULL this is the ordinary LASSO/SCAD/MCP lambda.max
# (SCAD and MCP have the same derivative at zero as the LASSO, so one path
# serves all three). The weighted form is used for the adaptive LASSO only.

.fg_lambda_path <- function(grad, n,
                            penalty_factor   = NULL,
                            nlambda          = 100L,
                            lambda_min_ratio = 0.001) {
  if (is.null(penalty_factor)) penalty_factor <- rep(1, length(grad))
  if (length(penalty_factor) != length(grad))
    stop("penalty_factor must have one value per predictor.")

  # The (1 + 1e-6) keeps the first grid point strictly above the threshold.
  # At exactly max|g|/n, rounding in fastCrrp's soft-threshold test can let the
  # top predictor in at ~1e-16, and the solver's RELATIVE convergence test
  # then burns max.iter sweeps on that coefficient (observed at n = 500,
  # p = 1000). The null score here agrees with fastCrrp's internal score to
  # ~1e-15, so this is rounding, not a mismatch.
  lambda_max <- (1 + 1e-6) * max(abs(grad) / penalty_factor) / n
  if (!is.finite(lambda_max) || lambda_max < 1e-8)
    stop("Degenerate lambda.max (", format(lambda_max, digits = 3),
         "): null-score projection collapsed. Inspect outcome/design.")

  10^seq(log10(lambda_max),
         log10(lambda_min_ratio * lambda_max),
         length.out = nlambda)
}

# ---------------------------------------------------------------------------
# .fit_fastcrrp -- shared fastCrrp wrapper
# ---------------------------------------------------------------------------
# Port of batch_run_complex.R line 253. Handles version guard (getBreslowJumps)
# and coefficient-matrix orientation normalisation.

.fit_fastcrrp <- function(x_std, ftime, status, predictor_names,
                           nlambda, lambda_min_ratio,
                           algorithm_tolerance, max_iterations,
                           penalty, penalty_factor = NULL,
                           gamma = NULL, lambda = NULL) {
  n          <- nrow(x_std)
  npredictors <- length(predictor_names)

  fit_data   <- data.frame(
    TTE    = ftime,
    Status = as.integer(status),
    x_std,
    check.names = FALSE
  )
  fg_formula <- stats::as.formula(
    paste0("fastcmprsk::Crisk(TTE, Status, cencode = 0, failcode = 1) ~ ",
           paste(predictor_names, collapse = " + "))
  )

  args <- list(
    formula          = fg_formula,
    data             = fit_data,
    eps              = algorithm_tolerance,
    max.iter         = max_iterations,
    standardize      = FALSE,
    penalty          = penalty,
    lambda.min.ratio = lambda_min_ratio,
    nlambda          = nlambda
  )

  if (!is.null(penalty_factor)) {
    if (length(penalty_factor) != npredictors ||
        any(!is.finite(penalty_factor)) ||
        any(penalty_factor <= 0))
      stop("penalty_factor must have one positive finite value per predictor.")
    args$penalty.factor <- as.numeric(penalty_factor)
  }
  if (!is.null(gamma))  args$gamma  <- gamma
  if (!is.null(lambda)) {
    if (any(!is.finite(lambda)) || any(lambda < 0))
      stop("lambda must contain finite non-negative values.")
    args$lambda <- as.numeric(lambda)
  }
  if ("getBreslowJumps" %in% names(formals(fastcmprsk::fastCrrp)))
    args$getBreslowJumps <- FALSE

  fit       <- do.call(fastcmprsk::fastCrrp, args)
  beta_path <- as.matrix(stats::coef(fit))

  if (nrow(beta_path) != npredictors && ncol(beta_path) == npredictors) {
    beta_path  <- t(beta_path)
    fit$coef   <- beta_path
  }
  if (nrow(beta_path) != npredictors)
    stop("Unexpected coef-matrix dimensions from fastCrrp().")

  rownames(fit$coef) <- predictor_names
  fit
}

# ---------------------------------------------------------------------------
# .select_by_bic -- BIC selection over a lambda path
# ---------------------------------------------------------------------------
# Port of batch_run_complex.R line 318. df = model size (number of nonzero
# coefs); see FUTURE WORK in the header for the fragility this introduces.
# Returns beta_std and beta_raw so score.r can choose the right scale.

.select_by_bic <- function(fit, method_label, predictor_names, nobs, x_scale) {
  npredictors <- length(predictor_names)
  beta_path   <- as.matrix(stats::coef(fit))
  if (nrow(beta_path) != npredictors && ncol(beta_path) == npredictors)
    beta_path <- t(beta_path)

  df_path  <- colSums(beta_path != 0)
  bic_path <- -2 * as.numeric(fit$logLik) + log(nobs) * df_path

  if (length(bic_path) != ncol(beta_path))
    stop("BIC path and coefficient path have incompatible lengths.")

  valid <- is.finite(bic_path)
  if (!is.null(fit$converged) && length(fit$converged) == length(bic_path))
    valid <- valid & as.logical(fit$converged)
  if (!any(valid))
    stop("No finite, converged solution available for ", method_label, ".")

  bic_sel          <- bic_path
  bic_sel[!valid]  <- Inf
  index            <- which.min(bic_sel)

  beta_std          <- as.numeric(beta_path[, index])
  names(beta_std)   <- predictor_names
  # x_std = (x - center) / scale, so divide by scale to recover original units.
  # Centering is absorbed by the baseline hazard.
  beta_raw          <- beta_std / x_scale
  names(beta_raw)   <- predictor_names

  list(
    method    = method_label,
    index     = index,
    lambda    = as.numeric(fit$lambda.path[index]),
    bic       = as.numeric(bic_path[index]),
    df        = as.integer(df_path[index]),
    converged = if (!is.null(fit$converged) &&
                    length(fit$converged) == length(bic_path))
                  as.logical(fit$converged[index]) else NA,
    beta_std  = beta_std,
    beta_raw  = beta_raw,
    selected  = beta_std != 0,
    bic_path  = bic_path,
    df_path   = df_path,
    # Smallest lambda with a usable solution. BIC selecting this index means
    # the minimum sits at the effective end of the path (either the grid ran
    # out or every smaller lambda failed to converge), so the "optimum" is
    # really a boundary.
    last_valid_index = max(which(valid)),
    # TRUE when some finite but unconverged lambda has a LOWER BIC than the
    # one chosen: the selection was cut short by solver failure rather than
    # settled by the criterion. Seen for SCAD/MCP at n = 500, p = 1000, where
    # BIC keeps falling as the fit approaches divergence.
    truncated     = any(is.finite(bic_path) & !valid & bic_path < bic_path[index])
  )
}

# ---------------------------------------------------------------------------
# .path_diag -- per-lambda solver diagnostics from a fastCrrp fit
# ---------------------------------------------------------------------------
# Kept separate from .select_by_bic so that a fit whose selection step fails
# (e.g. no converged lambda) still returns its path diagnostics.

.path_diag <- function(fit) {
  list(
    lambda_path    = as.numeric(fit$lambda.path),
    logLik_path    = as.numeric(fit$logLik),
    converged_path = as.logical(fit$converged),
    iter_path      = as.integer(fit$iter)
  )
}

# ---------------------------------------------------------------------------
# fit_comparators -- public entry point
# ---------------------------------------------------------------------------

fit_comparators <- function(x, y) {
  n           <- nrow(x)
  p           <- ncol(x)
  predictor_names <- colnames(x)
  if (is.null(predictor_names)) predictor_names <- paste0("X_", seq_len(p))

  ftime  <- y[, "time"]
  status <- y[, "status"]
  nobs   <- n

  # ---- Shared constants (mirror batch_run_complex.R line 248) ----
  nlambda             <- 100L
  lambda_min_ratio    <- 0.001
  algorithm_tolerance <- 1e-7
  max_iterations      <- 5000L

  # ---- aLASSO constants (see the B3b comment block for the rationale) ----
  # Stage 1, ridge initializer. ridge_lambda = ridge_info_multiplier *
  # median(null diagonal information). Multiplier 1 means that in an
  # orthogonal design every first-stage coefficient is shrunk to about half
  # of its unpenalized value: moderate, and fixed in advance rather than
  # chosen from the data.
  ridge_info_multiplier <- 1
  ridge_nlambda         <- 10L     # short warm-start path down to ridge_lambda
  init_tolerance        <- 1e-5    # weights need ~2 significant digits, not 7
  # Stage 2, weights.
  adaptive_power        <- 1
  weight_eps_fraction   <- 0.05    # eps = fraction * median nonzero |beta_init|
  weight_bounds         <- c(0.1, 25)  # applied AFTER median normalisation

  # ---- B1: Standardise once; share x_std across all comparators ----
  # Same block as batch_run_complex.R lines 224-232. Puts adaptive weights on
  # the same coefficient scale as the penalised fits and avoids
  # method-to-method preprocessing differences.
  t_setup_start <- proc.time()[["elapsed"]]

  x_center <- colMeans(x)
  x_scale  <- apply(x, 2L, stats::sd)
  bad_scale <- !is.finite(x_scale) | x_scale <= sqrt(.Machine$double.eps)
  if (any(bad_scale))
    stop("Predictors with zero or invalid standard deviation: ",
         paste(predictor_names[bad_scale], collapse = ", "))

  x_std            <- sweep(x, 2L, x_center, FUN = "-")
  x_std            <- sweep(x_std, 2L, x_scale,  FUN = "/")
  colnames(x_std)  <- predictor_names

  # ---- B2: Null score and the ordinary (unweighted) lambda path ----
  # LASSO, SCAD and MCP share reference_lambda: their all-zero threshold is
  # the same. aLASSO builds its own weighted path from the same null score in
  # B3b, so the score is computed once here.
  null_score <- .fg_null_score(x_std, ftime, status,
                               failcode = 1L, cencode = 0L)
  reference_lambda <- .fg_lambda_path(
    null_score$grad, n,
    nlambda          = nlambda,
    lambda_min_ratio = lambda_min_ratio
  )

  t_setup_end   <- proc.time()[["elapsed"]]
  setup_seconds <- as.numeric(t_setup_end - t_setup_start)

  # ---- Helper: wrap a model result into the B5 contract ----
  .make_result <- function(sel, fit_obj, secs, extra_meta = list()) {
    risk <- if (sel$converged %in% TRUE)
      as.numeric(x_std %*% sel$beta_std)
    else
      rep(NA_real_, n)

    list(
      coef_raw    = sel$beta_raw,
      coef_std    = sel$beta_std,
      lambda      = sel$lambda,
      bic         = sel$bic,
      risk        = risk,
      converged   = isTRUE(sel$converged),
      fit_seconds = secs,
      meta        = c(
        list(method      = sel$method,
             index       = sel$index,
             bic_path    = sel$bic_path,
             df_path     = sel$df_path,
             selected_df = sel$df,
             selected_at_last_valid = sel$index == sel$last_valid_index,
             selection_truncated    = sel$truncated),
        .path_diag(fit_obj),
        extra_meta
      )
    )
  }

  # extra_meta carries whatever diagnostics were gathered before the error,
  # so a failed fit still shows where along the path it broke down.
  .failed_result <- function(method, secs, extra_meta = list()) {
    list(
      coef_raw    = stats::setNames(rep(NA_real_, p), predictor_names),
      coef_std    = stats::setNames(rep(NA_real_, p), predictor_names),
      lambda      = NA_real_,
      bic         = NA_real_,
      risk        = rep(NA_real_, n),
      converged   = FALSE,
      fit_seconds = secs,
      meta        = c(list(method = method, error = TRUE), extra_meta)
    )
  }

  # Assignments inside a tryCatch() expression land in this frame, so each
  # *_diag list below is visible to the matching error handler.

  # ---- B3a: LASSO ----
  t0 <- proc.time()[["elapsed"]]
  lasso_diag <- list()
  lasso_result <- tryCatch({
    fit_lasso <- .fit_fastcrrp(
      x_std, ftime, status, predictor_names,
      nlambda, lambda_min_ratio, algorithm_tolerance, max_iterations,
      penalty = "LASSO", lambda = reference_lambda
    )
    lasso_diag <- .path_diag(fit_lasso)
    if (!isTRUE(all.equal(fit_lasso$lambda.path, reference_lambda, tolerance = 1e-12)))
      stop("LASSO fit did not use the reference lambda path.")
    sel <- .select_by_bic(fit_lasso, "LASSO", predictor_names, nobs, x_scale)
    .make_result(sel, fit_lasso, as.numeric(proc.time()[["elapsed"]] - t0))
  }, error = function(e) {
    warning("LASSO fit failed: ", conditionMessage(e))
    .failed_result("LASSO", as.numeric(proc.time()[["elapsed"]] - t0),
                   lasso_diag)
  })

  # ---- B3b: adaptive LASSO ----
  # Three stages, each a separate decision:
  #
  # 1. Initializer: ridge at ONE pre-specified penalty, for every p. The
  #    earlier scheme used unpenalized fastCrr when p < n and otherwise kept
  #    the ridge solution with the largest in-sample log-likelihood. Ridge
  #    log-likelihood rises monotonically as lambda falls, so that rule always
  #    picked the smallest lambda on the path, i.e. a near-unpenalized fit,
  #    which is unstable under correlated predictors whether or not p < n.
  #    Selecting the ridge penalty by inner CV would be better in principle but
  #    multiplies the aLASSO cost by the number of folds; the ridge path alone
  #    ran ~100 s at n = 500, p = 200. ridge_lambda is tied to the null
  #    information so it scales with the data rather than being a bare number.
  #    A short geometric path from 10 * ridge_lambda supplies warm starts.
  #
  # 2. Weights: w_j = 1 / (|beta_j| + eps)^power, eps a small fraction of the
  #    median nonzero |beta_j|, normalised by the median weight and THEN
  #    clipped to weight_bounds, so the largest-to-smallest penalty ratio is at
  #    most 250. (Previously the cap was applied before dividing by the minimum
  #    weight and could not bound the ratio; it was also inert, because the
  #    1e-4 floor already bounded raw weights at 1e4.) The lower bound keeps a
  #    few large first-stage coefficients from becoming effectively
  #    unpenalized, which is where coordinate descent stalls in correlated
  #    designs.
  #
  # 3. Lambda path: fastCrrp does not rescale penalty.factor, so the all-zero
  #    threshold for the weighted penalty is max_j |g_j| / (n w_j). The path is
  #    built from that value. Median normalisation alone does not change the
  #    fit (it only relabels lambda); the bounds and the weighted path do.
  t0 <- proc.time()[["elapsed"]]
  alasso_diag <- list()
  alasso_result <- tryCatch({
    # -- Stage 1: ridge initializer --
    ridge_lambda <- ridge_info_multiplier * stats::median(null_score$info)
    if (!is.finite(ridge_lambda) || ridge_lambda <= 0)
      stop("Non-positive ridge penalty from null information (",
           format(ridge_lambda, digits = 3), ").")
    ridge_path <- ridge_lambda * 10^seq(1, 0, length.out = ridge_nlambda)

    ridge_fit <- .fit_fastcrrp(
      x_std, ftime, status, predictor_names,
      ridge_nlambda, lambda_min_ratio, init_tolerance, max_iterations,
      penalty = "RIDGE", lambda = ridge_path
    )
    ridge_last       <- length(ridge_fit$lambda.path)
    initial_beta_std <- as.numeric(ridge_fit$coef[, ridge_last])
    names(initial_beta_std) <- predictor_names
    ridge_converged  <- isTRUE(as.logical(ridge_fit$converged[ridge_last]))
    alasso_diag <- list(
      initial_estimator = "ridge",
      ridge_lambda      = ridge_lambda,
      ridge_converged   = ridge_converged,
      ridge_iter        = as.integer(ridge_fit$iter[ridge_last]),
      initial_beta_std  = initial_beta_std
    )
    if (!all(is.finite(initial_beta_std)))
      stop("Ridge initializer returned non-finite coefficients.")
    # An unconverged ridge after max_iterations sweeps is still a usable
    # source of weight magnitudes; it is recorded, not treated as fatal.

    # -- Stage 2: stabilised, bounded weights --
    nonzero_init <- abs(initial_beta_std[initial_beta_std != 0])
    weight_eps   <- weight_eps_fraction * stats::median(nonzero_init)
    if (!is.finite(weight_eps) || weight_eps <= 0) weight_eps <- 0.01
    adaptive_weights <- 1 / (abs(initial_beta_std) + weight_eps)^adaptive_power
    adaptive_weights <- adaptive_weights / stats::median(adaptive_weights)
    adaptive_weights <- pmin(pmax(adaptive_weights, weight_bounds[1L]),
                             weight_bounds[2L])
    names(adaptive_weights) <- predictor_names
    alasso_diag <- c(alasso_diag, list(
      weight_eps       = weight_eps,
      adaptive_weights = adaptive_weights,
      weight_quantiles = stats::quantile(adaptive_weights,
                                         c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1)),
      n_weight_lower   = sum(adaptive_weights <= weight_bounds[1L]),
      n_weight_upper   = sum(adaptive_weights >= weight_bounds[2L])
    ))

    # -- Stage 3: weighted lambda path and final fit --
    alasso_lambda <- .fg_lambda_path(
      null_score$grad, n,
      penalty_factor   = adaptive_weights,
      nlambda          = nlambda,
      lambda_min_ratio = lambda_min_ratio
    )
    fit_alasso <- .fit_fastcrrp(
      x_std, ftime, status, predictor_names,
      nlambda, lambda_min_ratio, algorithm_tolerance, max_iterations,
      penalty = "LASSO", penalty_factor = adaptive_weights,
      lambda = alasso_lambda
    )
    alasso_diag <- c(alasso_diag, .path_diag(fit_alasso))
    if (!isTRUE(all.equal(fit_alasso$lambda.path, alasso_lambda, tolerance = 1e-12)))
      stop("aLASSO fit did not use its weighted lambda path.")
    sel <- .select_by_bic(fit_alasso, "adaptive LASSO", predictor_names, nobs,
                          x_scale)
    # .make_result adds the path diagnostics itself; pass only stage 1-2 here.
    .make_result(sel, fit_alasso,
                 as.numeric(proc.time()[["elapsed"]] - t0),
                 extra_meta = alasso_diag[setdiff(names(alasso_diag),
                                                  names(.path_diag(fit_alasso)))])
  }, error = function(e) {
    warning("aLASSO fit failed: ", conditionMessage(e))
    .failed_result("aLASSO", as.numeric(proc.time()[["elapsed"]] - t0),
                   alasso_diag)
  })

  # ---- B3c: SCAD (gamma = 3.7) ----
  # paper's fixed concavity parameter
  t0 <- proc.time()[["elapsed"]]
  scad_diag <- list(gamma = 3.7)
  scad_result <- tryCatch({
    fit_scad <- .fit_fastcrrp(
      x_std, ftime, status, predictor_names,
      nlambda, lambda_min_ratio, algorithm_tolerance, max_iterations,
      penalty = "SCAD", gamma = 3.7, lambda = reference_lambda
    )
    scad_diag <- c(scad_diag, .path_diag(fit_scad))
    if (!isTRUE(all.equal(fit_scad$lambda.path, reference_lambda, tolerance = 1e-12)))
      stop("SCAD fit did not use the reference lambda path.")
    sel <- .select_by_bic(fit_scad, "SCAD", predictor_names, nobs, x_scale)
    .make_result(sel, fit_scad, as.numeric(proc.time()[["elapsed"]] - t0),
                 extra_meta = list(gamma = 3.7))
  }, error = function(e) {
    warning("SCAD fit failed: ", conditionMessage(e))
    .failed_result("SCAD", as.numeric(proc.time()[["elapsed"]] - t0),
                   scad_diag)
  })

  # ---- B3d: MCP (gamma = 2.7) ----
  # paper's fixed concavity parameter
  t0 <- proc.time()[["elapsed"]]
  mcp_diag <- list(gamma = 2.7)
  mcp_result <- tryCatch({
    fit_mcp <- .fit_fastcrrp(
      x_std, ftime, status, predictor_names,
      nlambda, lambda_min_ratio, algorithm_tolerance, max_iterations,
      penalty = "MCP", gamma = 2.7, lambda = reference_lambda
    )
    mcp_diag <- c(mcp_diag, .path_diag(fit_mcp))
    if (!isTRUE(all.equal(fit_mcp$lambda.path, reference_lambda, tolerance = 1e-12)))
      stop("MCP fit did not use the reference lambda path.")
    sel <- .select_by_bic(fit_mcp, "MCP", predictor_names, nobs, x_scale)
    .make_result(sel, fit_mcp, as.numeric(proc.time()[["elapsed"]] - t0),
                 extra_meta = list(gamma = 2.7))
  }, error = function(e) {
    warning("MCP fit failed: ", conditionMessage(e))
    .failed_result("MCP", as.numeric(proc.time()[["elapsed"]] - t0),
                   mcp_diag)
  })

  list(
    lasso           = lasso_result,
    alasso          = alasso_result,
    scad            = scad_result,
    mcp             = mcp_result,
    setup_seconds   = setup_seconds,
    standardization = list(center = x_center, scale = x_scale)
  )
}

# ---------------------------------------------------------------------------
# fit_fg_oracle -- upper bound; NOT a penalised comparator
# ---------------------------------------------------------------------------

fit_fg_oracle <- function(x, y, beta1) {
  n              <- nrow(x)
  p              <- ncol(x)
  predictor_names <- colnames(x)
  if (is.null(predictor_names)) predictor_names <- paste0("X_", seq_len(p))

  ftime  <- y[, "time"]
  status <- as.integer(y[, "status"])

  # x_scale needed to compute coef_std = coef_raw * scale
  x_scale <- apply(x, 2L, stats::sd)

  active_idx <- which(beta1 != 0)

  t0 <- proc.time()[["elapsed"]]

  # ---- Zero-actives: well-formed degenerate result ----
  if (length(active_idx) == 0L) {
    zero_coef <- stats::setNames(rep(0.0, p), predictor_names)
    return(list(
      coef_raw    = zero_coef,
      coef_std    = zero_coef,
      risk        = rep(0.0, n),
      converged   = TRUE,
      fit_seconds = as.numeric(proc.time()[["elapsed"]] - t0),
      meta        = list(active_idx = integer(0), fit = NULL, degenerate = TRUE)
    ))
  }

  x_active     <- x[, active_idx, drop = FALSE]
  active_names <- predictor_names[active_idx]

  fit_args <- list(
    formula  = stats::as.formula(
      paste0("fastcmprsk::Crisk(ftime, status) ~ ",
             paste(active_names, collapse = " + "))
    ),
    variance = FALSE
  )
  if ("getBreslowJumps" %in% names(formals(fastcmprsk::fastCrr)))
    fit_args$getBreslowJumps <- FALSE

  oracle_env        <- new.env(parent = parent.frame())
  oracle_env$ftime  <- ftime
  oracle_env$status <- status
  for (nm in active_names)
    oracle_env[[nm]] <- x_active[, nm]
  fit_args$data     <- oracle_env

  fit <- do.call(fastcmprsk::fastCrr, fit_args)

  active_coef <- as.numeric(stats::coef(fit))
  if (length(active_coef) != length(active_idx))
    active_coef <- active_coef[seq_len(length(active_idx))]

  # Scatter into length-p vectors; non-active positions stay zero
  coef_raw          <- stats::setNames(rep(0.0, p), predictor_names)
  coef_raw[active_idx] <- active_coef

  coef_std          <- coef_raw * x_scale

  converged <- if (!is.null(fit$converged)) isTRUE(all(as.logical(fit$converged)))
               else NA

  list(
    coef_raw    = coef_raw,
    coef_std    = coef_std,
    risk        = as.numeric(x %*% coef_raw),
    converged   = converged,
    fit_seconds = as.numeric(proc.time()[["elapsed"]] - t0),
    meta        = list(active_idx = active_idx, fit = fit)
  )
}
