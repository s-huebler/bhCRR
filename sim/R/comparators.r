# sim/R/comparators.r -- LASSO, aLASSO, SCAD, MCP, and oracle comparator fits
#
# OWNS  (public functions)
#   fit_comparators(x, y)
#                    standardise x once, build a shared lambda path, fit the
#                    four penalised Fine-Gray models (LASSO, aLASSO, SCAD, MCP),
#                    select by BIC, and return the assembled result list (B5).
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
#   .fg_lambda_path(x_std, ftime, status, failcode, cencode, nlambda,
#                   lambda_min_ratio)
#                    build a clean geometric lambda grid from the KKT condition
#                    of the null Fine-Gray model; shared across all four
#                    penalised fits; needed because fastCrrp()'s built-in
#                    grid degenerates when there are no censored observations
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
#                                initial estimator (unpenalized fastCrr or
#                                ridge endpoint), since those weights are part
#                                of its cost. To compare one comparator against
#                                bhCRR's fit_seconds total, add setup_seconds
#                                once: comparator_total = setup_seconds +
#                                fit_seconds["alasso"] (for example)
#         meta        list(...)  method, index, bic_path, df_path, gamma,
#                                initial_estimator (aLASSO only)
#       )
#       a failed fit produces coef_raw = coef_std = rep(NA_real_, p),
#       risk = rep(NA_real_, n), and converged = FALSE
#     $setup_seconds   numeric(1)  standardization + shared lambda path
#     $standardization list(center = numeric p, scale = numeric p)
#
# RETURNS (fit_fg_oracle)
#   list(
#     coef        numeric k  named; coefficients on the active columns only
#                            (length = number of true actives, NOT p)
#     converged   logical(1) or NA if fastCrr does not report convergence
#     fit_seconds numeric(1)
#     meta        list(active_idx = integer k, fit = fastCrr object or NULL)
#   )
#   empty-active case: coef = named numeric(0); converged = TRUE vacuously;
#   fit = NULL
#
# FUTURE WORK
#   Lambda is currently selected by BIC over a shared geometric path.
#   A cross-validated selection exists at
#   sim/_legacy/code/other_comparators_cv_section.R (799 lines, Wolbers C-index
#   criterion) and is deliberately not ported; comparator tuning is an open
#   problem to revisit. The BIC path has at least two fragile points:
#     1. BIC penalises df = number of nonzero coefficients (model size), not
#        the effective-df trace; this overpenalises on correlated predictors.
#     2. The null-score lambda.max is computed from cause-1 events only; under
#        high censoring the grid may be too wide, leaving the BIC minimum at a
#        nearly-saturated model (df close to p).
#   When returning to comparator tuning, examine bic_path and df_path from
#   stored results to check whether the selected model is consistently near
#   the endpoints of the grid.
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
# .fg_lambda_path -- shared geometric lambda grid from KKT null score
# ---------------------------------------------------------------------------
# fastCrrp()'s built-in lambda grid degenerates when there are no censored
# observations (the IPCW/score projection collapses to roundoff noise). Build
# one clean grid from the KKT condition of the null Fine-Gray model and share
# it across all comparators. Port of batch_run_complex.R line 447.

.fg_lambda_path <- function(x_std, ftime, status,
                             failcode         = 1L,
                             cencode          = 0L,
                             nlambda          = 100L,
                             lambda_min_ratio = 0.001) {
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
  }

  lambda_max <- max(abs(grad)) / n
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
    df_path   = df_path
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

  # ---- B2: Shared reference lambda path ----
  reference_lambda <- .fg_lambda_path(
    x_std, ftime, status,
    failcode         = 1L,
    cencode          = 0L,
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
        list(method   = sel$method,
             index    = sel$index,
             bic_path = sel$bic_path,
             df_path  = sel$df_path),
        extra_meta
      )
    )
  }

  .failed_result <- function(method, secs) {
    list(
      coef_raw    = stats::setNames(rep(NA_real_, p), predictor_names),
      coef_std    = stats::setNames(rep(NA_real_, p), predictor_names),
      lambda      = NA_real_,
      bic         = NA_real_,
      risk        = rep(NA_real_, n),
      converged   = FALSE,
      fit_seconds = secs,
      meta        = list(method = method, error = TRUE)
    )
  }

  # ---- B3a: LASSO ----
  t0 <- proc.time()[["elapsed"]]
  lasso_result <- tryCatch({
    fit_lasso <- .fit_fastcrrp(
      x_std, ftime, status, predictor_names,
      nlambda, lambda_min_ratio, algorithm_tolerance, max_iterations,
      penalty = "LASSO", lambda = reference_lambda
    )
    if (!isTRUE(all.equal(fit_lasso$lambda.path, reference_lambda, tolerance = 1e-12)))
      stop("LASSO fit did not use the reference lambda path.")
    sel <- .select_by_bic(fit_lasso, "LASSO", predictor_names, nobs, x_scale)
    .make_result(sel, fit_lasso, as.numeric(proc.time()[["elapsed"]] - t0))
  }, error = function(e) {
    warning("LASSO fit failed: ", conditionMessage(e))
    .failed_result("LASSO", as.numeric(proc.time()[["elapsed"]] - t0))
  })

  # ---- B3b: adaptive LASSO ----
  # Initial estimator: unpenalised fastCrr (when p < n and convergence OK),
  # falling back to ridge fastCrrp endpoint. Port of batch_run_complex.R
  # lines 532-576. The initial estimator used here is the UNPENALIZED
  # Fine-Gray fit (fit_initial_finegray), NOT a LASSO first stage; the CV
  # variant in other_comparators_cv_section.R used a different scheme.
  t0 <- proc.time()[["elapsed"]]
  alasso_result <- tryCatch({
    initial_beta_std  <- NULL
    initial_estimator <- NULL

    if (p < n) {
      afit_data <- data.frame(
        TTE    = ftime,
        Status = as.integer(status),
        x_std,
        check.names = FALSE
      )
      init_args <- list(
        formula    = stats::as.formula(
          paste0("fastcmprsk::Crisk(TTE, Status, cencode = 0, failcode = 1) ~ ",
                 paste(predictor_names, collapse = " + "))
        ),
        data       = afit_data,
        eps        = algorithm_tolerance,
        max.iter   = max_iterations,
        standardize = FALSE,
        variance   = FALSE
      )
      if ("getBreslowJumps" %in% names(formals(fastcmprsk::fastCrr)))
        init_args$getBreslowJumps <- FALSE

      init_fit <- tryCatch(
        do.call(fastcmprsk::fastCrr, init_args),
        error = function(e) {
          warning("Unpenalized fastCrr() failed: ", conditionMessage(e))
          NULL
        }
      )

      if (!is.null(init_fit)) {
        candidate <- as.numeric(stats::coef(init_fit))
        init_conv <- is.null(init_fit$converged) ||
          all(as.logical(init_fit$converged))
        if (init_conv && length(candidate) == p && all(is.finite(candidate))) {
          initial_beta_std  <- candidate
          initial_estimator <- "unpenalized fastCrr"
        }
      }
    }

    if (is.null(initial_beta_std)) {
      warning("Using ridge-regularised initial estimate for adaptive-LASSO weights.")
      ridge_fit   <- .fit_fastcrrp(
        x_std, ftime, status, predictor_names,
        nlambda, lambda_min_ratio, algorithm_tolerance, max_iterations,
        penalty = "RIDGE", lambda = reference_lambda
      )
      ridge_beta  <- as.matrix(stats::coef(ridge_fit))
      ridge_valid <- rep(TRUE, ncol(ridge_beta))
      if (!is.null(ridge_fit$converged) &&
          length(ridge_fit$converged) == ncol(ridge_beta))
        ridge_valid <- as.logical(ridge_fit$converged)
      ridge_valid <- ridge_valid &
        apply(ridge_beta, 2L, function(z) all(is.finite(z)))
      if (!any(ridge_valid))
        stop("No usable ridge solution for adaptive-LASSO weights.")
      ridge_idx        <- which.max(
        replace(as.numeric(ridge_fit$logLik), !ridge_valid, -Inf)
      )
      initial_beta_std  <- as.numeric(ridge_beta[, ridge_idx])
      initial_estimator <- "ridge fastCrrp endpoint"
    }
    names(initial_beta_std) <- predictor_names

    adaptive_power    <- 1
    aLASSO_floor      <- 1e-4
    aLASSO_weight_cap <- 1e4
    adaptive_weights  <- 1 / pmax(abs(initial_beta_std), aLASSO_floor)^adaptive_power
    adaptive_weights  <- pmin(adaptive_weights, aLASSO_weight_cap)
    adaptive_weights  <- adaptive_weights / min(adaptive_weights)
    names(adaptive_weights) <- predictor_names

    fit_alasso <- .fit_fastcrrp(
      x_std, ftime, status, predictor_names,
      nlambda, lambda_min_ratio, algorithm_tolerance, max_iterations,
      penalty = "LASSO", penalty_factor = adaptive_weights,
      lambda = reference_lambda
    )
    if (!isTRUE(all.equal(fit_alasso$lambda.path, reference_lambda, tolerance = 1e-12)))
      stop("aLASSO fit did not use the reference lambda path.")
    sel <- .select_by_bic(fit_alasso, "adaptive LASSO", predictor_names, nobs,
                          x_scale)
    .make_result(sel, fit_alasso,
                 as.numeric(proc.time()[["elapsed"]] - t0),
                 extra_meta = list(initial_estimator = initial_estimator))
  }, error = function(e) {
    warning("aLASSO fit failed: ", conditionMessage(e))
    .failed_result("aLASSO", as.numeric(proc.time()[["elapsed"]] - t0))
  })

  # ---- B3c: SCAD (gamma = 3.7) ----
  # paper's fixed concavity parameter
  t0 <- proc.time()[["elapsed"]]
  scad_result <- tryCatch({
    fit_scad <- .fit_fastcrrp(
      x_std, ftime, status, predictor_names,
      nlambda, lambda_min_ratio, algorithm_tolerance, max_iterations,
      penalty = "SCAD", gamma = 3.7, lambda = reference_lambda
    )
    if (!isTRUE(all.equal(fit_scad$lambda.path, reference_lambda, tolerance = 1e-12)))
      stop("SCAD fit did not use the reference lambda path.")
    sel <- .select_by_bic(fit_scad, "SCAD", predictor_names, nobs, x_scale)
    .make_result(sel, fit_scad, as.numeric(proc.time()[["elapsed"]] - t0),
                 extra_meta = list(gamma = 3.7))
  }, error = function(e) {
    warning("SCAD fit failed: ", conditionMessage(e))
    .failed_result("SCAD", as.numeric(proc.time()[["elapsed"]] - t0))
  })

  # ---- B3d: MCP (gamma = 2.7) ----
  # paper's fixed concavity parameter
  t0 <- proc.time()[["elapsed"]]
  mcp_result <- tryCatch({
    fit_mcp <- .fit_fastcrrp(
      x_std, ftime, status, predictor_names,
      nlambda, lambda_min_ratio, algorithm_tolerance, max_iterations,
      penalty = "MCP", gamma = 2.7, lambda = reference_lambda
    )
    if (!isTRUE(all.equal(fit_mcp$lambda.path, reference_lambda, tolerance = 1e-12)))
      stop("MCP fit did not use the reference lambda path.")
    sel <- .select_by_bic(fit_mcp, "MCP", predictor_names, nobs, x_scale)
    .make_result(sel, fit_mcp, as.numeric(proc.time()[["elapsed"]] - t0),
                 extra_meta = list(gamma = 2.7))
  }, error = function(e) {
    warning("MCP fit failed: ", conditionMessage(e))
    .failed_result("MCP", as.numeric(proc.time()[["elapsed"]] - t0))
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
  ftime  <- y[, "time"]
  status <- as.integer(y[, "status"])

  active_idx <- which(beta1 != 0)

  t0 <- proc.time()[["elapsed"]]

  if (length(active_idx) == 0L) {
    return(list(
      coef        = stats::setNames(numeric(0), character(0)),
      converged   = TRUE,
      fit_seconds = as.numeric(proc.time()[["elapsed"]] - t0),
      meta        = list(active_idx = integer(0), fit = NULL)
    ))
  }

  x_active     <- x[, active_idx, drop = FALSE]
  active_names <- colnames(x_active)

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

  coef_raw <- as.numeric(stats::coef(fit))
  if (length(coef_raw) != length(active_idx))
    coef_raw <- coef_raw[seq_len(length(active_idx))]
  names(coef_raw) <- active_names

  converged <- if (!is.null(fit$converged)) isTRUE(all(as.logical(fit$converged)))
               else NA

  list(
    coef        = coef_raw,
    converged   = converged,
    fit_seconds = as.numeric(proc.time()[["elapsed"]] - t0),
    meta        = list(active_idx = active_idx, fit = fit)
  )
}
