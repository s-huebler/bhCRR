# sim/R/fit.r -- tuning + final fit for one simulated dataset
#
# OWNS
#   fit_bhcrr(x, y, scenario, run_id,
#             control = NULL, s0_seq = NULL, s1_seq = NULL, max_widen = 1L)
#                    run bhcrr_autotune() on an untuned bhCRR fit, refit at
#                    the chosen (s0, s1), and return everything score.r needs
#
# INPUTS
#   x          matrix n x p   design matrix from sim/R/dgp.r, column names
#                             X_1..X_p
#   y          matrix n x 2   outcome from sim/R/dgp.r (colnames time, status)
#   scenario   list           validated bhcrr_scenario; used for
#                               $zero_gap_target -> beta_min
#                               $u_min, $u_max    -> eval_time (fixed horizon)
#   run_id     integer(1)     run index; used only in error messages
#   control    bhcrr_cv_control or NULL
#                             NULL builds the default; pass an overriding
#                             object in tests or when changing CV parameters
#   s0_seq     numeric or NULL  override the default s0 grid from preflight
#   s1_seq     numeric or NULL  override the default s1 grid from preflight
#   max_widen  integer(1)     max grid-widening attempts (default 1L)
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
#     $control      bhcrr_cv_control  the control actually used; stored so
#                                 tests can verify hyperparameter mirroring
#                                 without inspecting internals
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
#   init_method must be kept identical in TWO places:
#     1. The direct fit_ssl_psdh() calls (init fit and final fit)
#     2. control$init_method -- the CV fold fits read it from there
#   They are separate code paths that feed the same initializer. A mismatch
#   means the untuned fit and the fold fits start from different initial
#   estimates, which biases the surface that autotune searches.
#
# FG_ORACLE NOTE
#   fg_oracle is NOT fitted here. It is fitted in sim/R/comparators.r by
#   fit_fg_oracle(x, y, beta1), which is kept separate from fit_comparators()
#   precisely because it consumes ground truth (the true active set). Migration
#   step 6 implements it. Legacy's fg_true object appears at
#   sim/_legacy/code/batch_run_complex.R:173.
#
# ERROR POLICY
#   Errors propagate to the caller. Per-run failure handling belongs in
#   run_one.r and the batch driver (migration steps 9 and 10). No tryCatch
#   wrapping here.
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
                      control   = NULL,
                      s0_seq    = NULL,
                      s1_seq    = NULL,
                      max_widen = 1L) {

  # ---- Shared hyperparameters (define once; reused in all three fit paths) ----
  initial_sparsity  <- 0.05
  theta_a           <- 1
  theta_b           <- 1      # deliberate legacy choice; fit_ssl_psdh default is
                               # ncol(x), which changes the prior relative to
                               # batch_run_complex.R -- do not "fix" to default
  epsilon           <- 1e-04
  inner_maxit_start <- 1000
  init_method       <- "LASSO_cv"
  init_args         <- list(lambda_path =
                              10^seq(log10(0.1), log10(0.001), length = 10))

  # ---- eval_time: scenario-fixed, not data-dependent (see EVAL_TIME NOTE) ----
  tau <- (scenario$u_min + scenario$u_max) / 2

  # ---- Control (build default if not supplied) ----
  if (is.null(control)) {
    control <- bhcrr_cv_control(
      nfolds      = 10L,
      ncv         = 2L,
      init_method = init_method,
      init_args   = init_args,
      # fold fits receive every hyperparameter through fit_args; without
      # theta_b = 1 here, folds run at theta_b = ncol(x) and autotune would
      # select (s0, s1) under a different prior than the final model uses
      fit_args    = list(
        initial_sparsity  = initial_sparsity,
        theta_a           = theta_a,
        theta_b           = theta_b,
        epsilon           = epsilon,
        inner_maxit_start = inner_maxit_start,
        maxit             = 50L
      ),
      parallel    = FALSE,   # forking is unreliable in Positron; CHPC array
                              # tasks use one cpu (SB_CPUS=1 in config.sh) so
                              # there is nothing to fork onto; if that ever
                              # changes, set workers from the SLURM allocation
                              # rather than detectCores()
      eval_time   = tau
    )
  }

  # ---- Stage 1: untuned fit to estimate theta for preflight ----
  t0 <- proc.time()[["elapsed"]]

  mod_initial <- fit_ssl_psdh(
    x, y,
    ss                = c(0.04, 0.6),
    initial_sparsity  = initial_sparsity,
    theta_a           = theta_a,
    theta_b           = theta_b,
    maxit             = 50L,
    epsilon           = epsilon,
    init_method       = init_method,
    init_args         = init_args,
    inner_maxit_start = inner_maxit_start
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
    initial_sparsity  = initial_sparsity,
    theta_a           = theta_a,
    theta_b           = theta_b,
    maxit             = 100L,    # legacy used 100 here vs 50 for the untuned fit
    epsilon           = epsilon,
    init_method       = init_method,
    init_args         = init_args,
    inner_maxit_start = inner_maxit_start
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
    eval_time    = tau,
    fit_seconds  = fit_seconds,
    timing       = timing,
    control      = control
  )
}
