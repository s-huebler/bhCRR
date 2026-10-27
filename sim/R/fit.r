# sim/R/fit.r -- tuning + final fit for one simulated dataset
#
# OWNS
#   fit_bhcrr(x, y, scenario, run_id,
#             control = NULL, s0_seq = NULL, s1_seq = NULL, max_widen = 1L,
#             maxit_folds = 50L, maxit_final = 100L)
#                    run bhcrr_autotune() on an untuned bhCRR fit, refit at
#                    the chosen (s0, s1), and return everything score.r needs
#
# INPUTS
#   x           matrix n x p   design matrix from sim/R/dgp.r, column names
#                              X_1..X_p
#   y           matrix n x 2   outcome from sim/R/dgp.r (colnames time, status)
#   scenario    list           validated bhcrr_scenario; used for
#                                $zero_gap_target -> beta_min
#                                $u_min, $u_max    -> eval_time (fixed horizon)
#   run_id      integer(1)     run index; used only in error messages
#   control     bhcrr_cv_control or NULL
#                              NULL builds the default (nfolds=10, ncv=2).
#                              Supplied controls are MERGED: the five mirrored
#                              hyperparameters are filled from defaults into
#                              control$fit_args so fold fits and direct fits are
#                              aligned by construction. See MIRRORING NOTE.
#   s0_seq      numeric or NULL  override the default s0 grid from preflight
#   s1_seq      numeric or NULL  override the default s1 grid from preflight
#   max_widen   integer(1)     max grid-widening attempts (default 1L)
#   maxit_folds integer(1)     EM iterations for the initial fit and each fold
#                              fit; default 50L. If control$fit_args$maxit is
#                              present it takes precedence over this formal.
#   maxit_final integer(1)     EM iterations for the final refit; default 100L.
#                              Kept independent of maxit_folds so a cheap
#                              interactive run can truncate the CV folds for
#                              speed while the final model still converges.
#
# RETURNS
#   list with fields:
#     $coef         numeric p     final coefficient vector; named with
#                                 colnames(x), length ncol(x)
#     $s0           numeric(1)    tuned spike scale
#     $s1           numeric(1)    tuned slab scale
#     $autotune     bhcrr_autotune  full autotune object (diagnostics)
#     $fit          ssl_psdh      final fit at chosen (s0, s1)
#     $mod_initial  ssl_psdh      untuned fit used to estimate theta and to
#                                 seed preflight; saved per the save-everything-
#                                 for-now decision in run_one.r's header
#     $converged    logical(1)    convergence flag of the final fit
#     $all_pass     logical(1)    from autotune$all_pass
#     $flags        character     from autotune$flags
#     $n_attempts   integer(1)    from autotune$n_attempts
#     $eval_time    numeric(1)    horizon at which CV scored each fold; also
#                                 the tau score.r must pass to wolbers_c() so
#                                 that tuning and reporting use the same metric
#     $fit_seconds  numeric(1)    total wall seconds for all three stages
#                                 (init_fit + autotune + final_fit); numeric,
#                                 never difftime. Total-vs-total is the honest
#                                 comparison: comparators' timings include their
#                                 lambda CV, so bhCRR's timing includes its CV.
#     $timing       list          per-stage wall seconds, each numeric:
#                                   $init_fit   untuned fit
#                                   $autotune   bhcrr_autotune() including CV
#                                   $final_fit  final refit at best (s0, s1)
#     $control      bhcrr_cv_control  the resolved, merged control actually
#                                 used; contains the five mirrored fit_args,
#                                 eval_time, init_method, and init_args as
#                                 applied to every fit path
#
# MIRRORING NOTE
#   Five hyperparameters must be identical in the fold fits (inside bhcrr_cv)
#   and the two direct fit_ssl_psdh() calls:
#     initial_sparsity = 0.05, theta_a = 1, theta_b = 1 (legacy choice; NOT
#     ncol(x)), epsilon = 1e-04, inner_maxit_start = 1000.
#   theta_b in particular: fit_ssl_psdh's default is ncol(x), so omitting it
#   from fit_args silently changes the prior. With the old "accept wholesale"
#   approach, a supplied control whose fit_args omitted theta_b would cause
#   fold fits to tune at theta_b=ncol(x) while the final model ran at theta_b=1.
#   The merge makes consistency structural: a caller can override any mirrored
#   value via control$fit_args and it propagates everywhere.
#
# EVAL_TIME NOTE
#   Legacy bhcrr_autotune() (sim/_legacy/code/batch_run_complex.R:192-198) had
#   no control argument and no eval_time field; the commented-out
#   select_by_wolbers used eval_quantile = 0.5 (median of cause-1 event times)
#   as a quantile fallback. A data-dependent horizon is not acceptable for
#   manuscript runs because the same (s0,s1) pair is evaluated on different
#   scales across simulation replicates. This file fixes the horizon at
#     eval_time = (scenario$u_min + scenario$u_max) / 2
#   This value is scenario-determined (not run-determined), lies at the centre
#   of the Fine-Gray event window, and coincides approximately with the median
#   cause-1 event time under the uniform scenario construction. score.r must
#   use the same value (returned as $eval_time).
#
# INIT_METHOD NOTE
#   init_method and init_args are read from the resolved control for ALL THREE
#   fit paths (init fit, fold fits, final fit). Before this refactor, the direct
#   fits used file-level locals and only the fold fits read from the control,
#   making bhcrr_cv_control(init_method="LASSO_bic") apply only to the folds.
#   Now the control is the single authoritative source for initialization.
#
# FG_ORACLE NOTE
#   fg_oracle is NOT fitted here. It is fitted in sim/R/comparators.r by
#   fit_fg_oracle(x, y, beta1), which is kept separate from fit_comparators()
#   precisely because it consumes ground truth (the true active set). Migration
#   step 6 implements it. Legacy's fg_true object appears at
#   sim/_legacy/code/batch_run_complex.R:173.
#
# ERROR POLICY
#   Errors propagate to the caller. Per-run failure handling belongs in the
#   driver (sim/local/batch.R for local runs; sim/chpc/R/run_array.r for
#   CHPC array jobs). No tryCatch wrapping here.
#
# DEPENDS ON
#   bhCRR::fit_ssl_psdh, bhCRR::bhcrr_autotune, bhCRR::bhcrr_cv_control
#
# STABLE SURFACE CONTRACT
#   This file reads exactly five fields from the bhcrr_autotune object:
#     $best$s0, $best$s1, $all_pass, $flags, $n_attempts
#   All other fields are internal to the bridge and may change without notice.
#   See design/autotune-bridge.md §4 (STABLE SURFACE) for the rationale.
#
# REPLACES
#   sim/_legacy/code/batch_run_complex.R lines ~176-215 (SSL timing + fit)
#   sim/_legacy/code/single_run_complex.R (entire file; same logic, one run)

# ---------------------------------------------------------------------------
# fit_bhcrr
# ---------------------------------------------------------------------------

fit_bhcrr <- function(x, y, scenario, run_id,
                      control     = NULL,
                      s0_seq      = NULL,
                      s1_seq      = NULL,
                      max_widen   = 1L,
                      maxit_folds = 50L,
                      maxit_final = 100L) {

  # ---- Mirrored hyperparameter defaults ----
  # Applied to all three fit paths. A supplied control$fit_args overrides these
  # selectively; unset entries get the values below.
  .mirror_defaults <- list(
    initial_sparsity  = 0.05,
    theta_a           = 1,
    theta_b           = 1,      # deliberate legacy choice; fit_ssl_psdh default is
                                 # ncol(x), which changes the prior relative to
                                 # batch_run_complex.R -- do not "fix" to default
    epsilon           = 1e-04,
    inner_maxit_start = 1000
  )

  # ---- eval_time: scenario-fixed, not data-dependent (see EVAL_TIME NOTE) ----
  tau <- (scenario$u_min + scenario$u_max) / 2

  # ---- Resolve control: build default or merge into a supplied one ----
  # After this block, control$fit_args carries all five mirrored entries plus
  # maxit = maxit_folds, so fold fits and direct fits read from the same source.
  if (is.null(control)) {
    control <- bhcrr_cv_control(
      nfolds      = 10L,
      ncv         = 2L,
      init_method = "LASSO_cv",
      init_args   = list(lambda_path =
                           10^seq(log10(0.1), log10(0.001), length = 10)),
      fit_args    = c(.mirror_defaults, list(maxit = maxit_folds)),
      parallel    = FALSE,   # forking unreliable in Positron; CHPC array tasks
                              # use one cpu so there is nothing to fork onto
      eval_time   = tau
    )
  } else {
    # Merge: mirrored defaults < control$fit_args (minus maxit).
    # maxit in control$fit_args overrides maxit_folds for fold fits;
    # maxit_final stays independent so a cheap run can truncate folds for speed
    # while the final model still converges properly.
    fa <- if (is.null(control$fit_args)) list() else control$fit_args
    if (!is.null(fa$maxit)) maxit_folds <- fa$maxit
    mirrored         <- modifyList(.mirror_defaults, fa[setdiff(names(fa), "maxit")])
    control$fit_args <- c(mirrored, list(maxit = maxit_folds))

    # init_method / init_args: apply defaults if not set in the supplied control
    if (is.null(control$init_method))
      control$init_method <- "LASSO_cv"
    if (!("lambda_path" %in% names(control$init_args)))
      control$init_args[["lambda_path"]] <-
        10^seq(log10(0.1), log10(0.001), length = 10)

    # eval_time: honor the supplied value; fall back to tau
    if (is.null(control$eval_time))
      control$eval_time <- tau
  }

  # ---- Stage 1: untuned fit to estimate theta for preflight ----
  t0 <- proc.time()[["elapsed"]]

  mod_initial <- fit_ssl_psdh(
    x, y,
    ss                = c(0.04, 0.6),
    initial_sparsity  = control$fit_args$initial_sparsity,
    theta_a           = control$fit_args$theta_a,
    theta_b           = control$fit_args$theta_b,
    maxit             = maxit_folds,
    epsilon           = control$fit_args$epsilon,
    init_method       = control$init_method,
    init_args         = control$init_args,
    inner_maxit_start = control$fit_args$inner_maxit_start
  )

  t1 <- proc.time()[["elapsed"]]

  # ---- Stage 2: autotune ----
  autotune <- bhcrr_autotune(
    mod_initial,
    beta_min   = scenario$zero_gap_target,
    beta_floor = 0.01,
    s0_seq     = s0_seq,
    s1_seq     = s1_seq,
    control    = control,
    max_widen  = max_widen
  )

  t2 <- proc.time()[["elapsed"]]

  # ---- Guard: autotune must return a usable (s0, s1) ----
  best_s0 <- autotune$best$s0
  best_s1 <- autotune$best$s1

  if (is.null(autotune$best) || is.null(best_s0) || is.null(best_s1) ||
      !is.finite(best_s0) || !is.finite(best_s1))
    stop("bhcrr_autotune() returned no usable (s0, s1) for run_id = ", run_id,
         ". Check autotune$flags for diagnostics.")

  # ---- Stage 3: final refit at chosen (s0, s1) ----
  final_fit <- fit_ssl_psdh(
    x, y,
    ss                = c(best_s0, best_s1),
    initial_sparsity  = control$fit_args$initial_sparsity,
    theta_a           = control$fit_args$theta_a,
    theta_b           = control$fit_args$theta_b,
    maxit             = maxit_final,
    epsilon           = control$fit_args$epsilon,
    init_method       = control$init_method,
    init_args         = control$init_args,
    inner_maxit_start = control$fit_args$inner_maxit_start
  )

  t3 <- proc.time()[["elapsed"]]

  # ---- Timing ----
  timing <- list(
    init_fit  = as.numeric(t1 - t0),
    autotune  = as.numeric(t2 - t1),
    final_fit = as.numeric(t3 - t2)
  )
  fit_seconds <- as.numeric(t3 - t0)

  # ---- Assemble return ----
  # ssl_psdh stores coefficients as a data.frame($Variable, $Estimate); extract
  # as a named numeric vector so score.r and comparators see the same interface.
  coef_df <- final_fit$coefficients
  coef    <- stats::setNames(as.numeric(coef_df$Estimate),
                              as.character(coef_df$Variable))

  list(
    coef         = coef,
    s0           = best_s0,
    s1           = best_s1,
    autotune     = autotune,
    fit          = final_fit,
    mod_initial  = mod_initial,
    converged    = final_fit$conv,   # ssl_psdh uses $conv, not $converged
    all_pass     = autotune$all_pass,
    flags        = autotune$flags,
    n_attempts   = autotune$n_attempts,
    eval_time    = control$eval_time,
    fit_seconds  = fit_seconds,
    timing       = timing,
    control      = control
  )
}
