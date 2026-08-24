#################### CHPC Batch Runner ################
# CHPC/SLURM-portable adaptation of Sim/code/batch_run_complex.R.
#
# All configuration comes from environment variables (set by config.sh and
# passed through by run_array.slurm). Nothing is hardcoded to the Mac.
#
#   REPO_ROOT        path to the repo checkout on CHPC
#   SCRATCH_RUN_DIR  shared-scratch PARENT dir; each npredictors value gets its
#                    own p<P> subdirectory beneath it (see the sweep at the end)
#   NOBS             sample size
#   NPREDICTORS      comma list, e.g. "25,50,221"
#   RUN_SUBSTART     first run index this array task handles (from run_array.slurm)
#   RUN_SUBEND       last  run index this array task handles
#   ZERO_GAP_TARGET  clinically-relevant minimum treatment effect
#   BETA1_ACTIVE     comma list of active cause-1 coefficients
#   BETA2_ACTIVE     comma list OR sentinel "neg_beta1" (default); cause-2 coefficients
#   ACTIVE_BLOCK     comma list of block labels aligned with BETA1_ACTIVE
#                      e.g. "strong,strong,weak,weak,indep"
#   BLOCK_PROPS      named num "name=value,...", block proportions summing to 1
#                      e.g. "indep=0.50,weak=0.25,strong=0.25"
#   BLOCK_RHO        named num "name=value,...", within-block correlations
#                      e.g. "indep=0.00,weak=0.30,strong=0.60"
#   LATENT_TYPE      latent variable distribution, default "gaussian"
#   LATENT_Q         named num "name=value,...", latent quantile per block
#                      e.g. "indep=0.5,weak=0.4,strong=0.5"
#   CENS_RATE        Exp rate for random censoring (0 = none), default "0.05"
#   U_MIN            lower bound for uniform observation window, default "100"
#   U_MAX            upper bound for uniform observation window, default "100"
#
# Files are named n<nobs>_p<p>_run<run>.Rdata, one object `result` each -- the
# exact contract chpc/R/parse_batch.R (and the original parse_batch_run.R)
# expects. set.seed() is intentionally NOT called: each run draws fresh data.


#################### Config from environment ################
get_env <- function(name, default = NULL) {
  v <- Sys.getenv(name, unset = NA_character_)
  if (is.na(v) || v == "") {
    if (is.null(default)) stop("Required environment variable not set: ", name)
    return(default)
  }
  v
}
num_list <- function(s) as.numeric(strsplit(s, ",", fixed = TRUE)[[1]])
int_list <- function(s) as.integer(strsplit(s, ",", fixed = TRUE)[[1]])

# Parse "name=value,name=value" into a named numeric vector.
named_num_list <- function(s) {
  pairs <- strsplit(s, ",", fixed = TRUE)[[1]]
  kv    <- strsplit(pairs, "=", fixed = TRUE)
  nms   <- vapply(kv, `[[`, character(1), 1L)
  vals  <- as.numeric(vapply(kv, `[[`, character(1), 2L))
  setNames(vals, nms)
}
# Parse "a,b,c" into a character vector.
chr_list <- function(s) strsplit(s, ",", fixed = TRUE)[[1]]

repo_root        <- get_env("REPO_ROOT")
out_dir          <- get_env("SCRATCH_RUN_DIR")
nobs             <- as.integer(get_env("NOBS", "200"))
npredictors_grid <- int_list(get_env("NPREDICTORS", "221"))
run_start        <- as.integer(get_env("RUN_SUBSTART", get_env("RUN_START", "1")))
run_end          <- as.integer(get_env("RUN_SUBEND",   get_env("RUN_END",   "10")))
zero_gap_target  <- as.numeric(get_env("ZERO_GAP_TARGET", "0.1"))
beta1_active     <- num_list(get_env("BETA1_ACTIVE", "0.40,-0.50,0.60,0.75,-0.80"))

beta2_active_raw <- get_env("BETA2_ACTIVE", "neg_beta1")
if (beta2_active_raw == "neg_beta1") {
  message("  BETA2_ACTIVE: sentinel 'neg_beta1' -- will use -beta1_active each run")
} else {
  beta2_active <- num_list(beta2_active_raw)
  message("  BETA2_ACTIVE: explicit list -- ", paste(beta2_active, collapse = ", "))
}

active_block <- chr_list(get_env("ACTIVE_BLOCK",
                                 "strong,strong,weak,weak,indep"))
block_props  <- named_num_list(get_env("BLOCK_PROPS",
                                       "indep=0.50,weak=0.25,strong=0.25"))
block_rho    <- named_num_list(get_env("BLOCK_RHO",
                                       "indep=0.00,weak=0.30,strong=0.60"))
latent_type  <- get_env("LATENT_TYPE", "gaussian")
latent_q     <- named_num_list(get_env("LATENT_Q",
                                       "indep=0.5,weak=0.4,strong=0.5"))
cens_rate    <- as.numeric(get_env("CENS_RATE", "0.05"))
u_min        <- as.numeric(get_env("U_MIN", "100"))
u_max        <- as.numeric(get_env("U_MAX", "100"))

# Validate block consistency: BLOCK_PROPS, BLOCK_RHO and LATENT_Q must carry
# the same block names, and every entry in ACTIVE_BLOCK must be one of them.
block_names_ref <- sort(names(block_props))
if (!identical(sort(names(block_rho)), block_names_ref)) {
  stop("BLOCK_RHO names (", paste(sort(names(block_rho)), collapse = ", "),
       ") do not match BLOCK_PROPS names (", paste(block_names_ref, collapse = ", "), ").")
}
if (!identical(sort(names(latent_q)), block_names_ref)) {
  stop("LATENT_Q names (", paste(sort(names(latent_q)), collapse = ", "),
       ") do not match BLOCK_PROPS names (", paste(block_names_ref, collapse = ", "), ").")
}
bad_active <- setdiff(active_block, block_names_ref)
if (length(bad_active) > 0) {
  stop("ACTIVE_BLOCK contains unknown block name(s): ",
       paste(bad_active, collapse = ", "),
       ". Known blocks: ", paste(block_names_ref, collapse = ", "), ".")
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message(sprintf("run_batch.R  repo=%s", repo_root))
message(sprintf("  out_dir=%s", out_dir))
message(sprintf("  nobs=%d  npredictors=%s  runs=%d..%d  zero_gap=%.3g",
                nobs, paste(npredictors_grid, collapse = ","),
                run_start, run_end, zero_gap_target))
message(sprintf("  active_block=%s  cens_rate=%.3g  u_min=%.3g  u_max=%.3g",
                paste(active_block, collapse = ","), cens_rate, u_min, u_max))


#################### Set Up (load everything once) ################
suppressPackageStartupMessages({
  library(tidyverse)
  library(fastcmprsk)   # provided by the renv library (renv.lock pins it)
  library(survival)
})

source(file.path(repo_root, "R/helpers.r"))
source(file.path(repo_root, "R/utils.r"))
source(file.path(repo_root, "R/dedupe_warnings.r"))

source(file.path(repo_root, "R/update_betas.r"))
source(file.path(repo_root, "R/expected_inclusion_probs.r"))
source(file.path(repo_root, "R/expected_penalty_weights.r"))
source(file.path(repo_root, "R/update_mixture_prob.r"))

Rcpp::sourceCpp(file.path(repo_root, "src/cv_fastcrrp.cpp"))
Rcpp::sourceCpp(file.path(repo_root, "src/RcppExports.cpp"))
source(file.path(repo_root, "R/cv_fastCrrp_cpp.r"))
source(file.path(repo_root, "R/fit_ssl_psdh.r"))

source(file.path(repo_root, "R/generate_foldid.r"))
source(file.path(repo_root, "R/predict_from_ssl_psdh.r"))
source(file.path(repo_root, "R/wolbers_c.r"))
source(file.path(repo_root, "R/cv_ssl_psdh.r"))
source(file.path(repo_root, "R/tune_ssl_psdh.r"))
source(file.path(repo_root, "R/threshold.R"))
source(file.path(repo_root, "R/tuning_diagnostics.r"))

# Block-correlated design matrix generator. Lives outside R/ on purpose:
# it is simulation-only code, not package code.
source(file.path(repo_root, "Sim/code/generate_block_X.R"))

# Start from a clean, non-deterministic RNG state.
if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)


#################### One Run ################
# Port of run_once() from Sim/code/batch_run_complex.R -- generates fresh data,
# fits the initial model, autotunes, refits, and returns a bundled result list.
run_once <- function(nobs, npredictors, beta1_active, beta2_active,
                     active_block, block_props, block_rho, latent_type,
                     latent_q, cens_rate, u_min, u_max,
                     zero_gap_target) {

  ## ----- Data generation -----
  design <- generate_block_X(nobs          = nobs,
                             npredictors   = npredictors,
                             beta1_active  = beta1_active,
                             beta2_active  = beta2_active,
                             active_block  = active_block,
                             block_props   = block_props,
                             block_rho     = block_rho,
                             latent        = latent_type,
                             latent_q      = latent_q,
                             unit_variance = TRUE,
                             layout        = "active_first")

  beta1    <- design$beta1
  beta2    <- design$beta2
  X_design <- design$X

  sim <- simulateTwoCauseFineGrayModel(nobs = nobs,
                                       beta1 = beta1,
                                       beta2 = beta2,
                                       X = X_design,
                                       u.min = u_min,
                                       u.max = u_max,
                                       p = 0.5,
                                       returnX = TRUE)

  # Guard against the generator substituting its own design matrix, which would
  # make the event times unrelated to the covariates we are analysing.
  if (!isTRUE(all.equal(unname(as.matrix(sim$X)), unname(X_design),
                        tolerance = 1e-12))) {
    stop("simulateTwoCauseFineGrayModel() did not return the supplied design ",
         "matrix unchanged; the event times were not generated from X_design.")
  }

  if (!is.finite(cens_rate) || cens_rate < 0) {
    stop("cens_rate must be a finite, non-negative Exponential rate.")
  }
  if (cens_rate > 0) {
    cens_time <- stats::rexp(nobs, rate = cens_rate)
    is_censored <- cens_time < sim$ftime
    sim$ftime[is_censored]   <- cens_time[is_censored]
    sim$fstatus[is_censored] <- 0
  }
  event_mix <- c(censored = mean(sim$fstatus == 0),
                 cause1   = mean(sim$fstatus == 1),
                 cause2   = mean(sim$fstatus == 2))

  sim_data <- cbind(data.frame(ID = 1:nobs,
                               TTE = sim$ftime,
                               Status = as.integer(sim$fstatus)),
                    X_design)

  names(sim_data) <- c("ID", "TTE", "Status", design$meta$name)

  x <- as.matrix(sim_data %>% select(starts_with("X")))
  y <- as.matrix(sim_data %>%
                   mutate(Status = as.numeric(Status)) %>%
                   select(TTE, Status))

  ## ----- FG model -----
  active_idx <- which(beta1 != 0)
  active_x <- x[, active_idx, drop = FALSE]

  fg_true <- fastcmprsk::fastCrr(fastcmprsk::Crisk(sim$ftime, as.integer(sim$fstatus)) ~ active_x,
                                 variance = FALSE)

  ## ----- SSL model -----
  ssl_psdh_time1 <- Sys.time()

  mod <- fit_ssl_psdh(x, y,
                      ss = c(0.04, 0.6),
                      initial_sparsity = 0.05,
                      theta_a = 1,
                      theta_b = 1,
                      maxit = 50,
                      epsilon = 1e-04,
                      init = NULL,
                      init_lam_path = 10^seq(log10(0.1), log10(0.001), length = 10),
                      inner_maxit_start = 1000)

  autotune <- bhcrr_autotune(mod,
                             beta_min   = zero_gap_target,
                             beta_floor = 0.01,
                             nfolds     = 10,
                             ncv        = 2,
                             foldid     = NULL,
                             max_widen  = 1)

  best_s0 <- autotune$best$s0
  best_s1 <- autotune$best$s1

  final_mod <- fit_ssl_psdh(x, y,
                            ss = c(best_s0, best_s1),
                            initial_sparsity = 0.05,
                            theta_a = 1,
                            theta_b = 1,
                            maxit = 100,
                            epsilon = 1e-04,
                            init = NULL,
                            init_lam_path = 10^seq(log10(0.1), log10(0.001), length = 10),
                            inner_maxit_start = 1000)

  ssl_psdh_time2 <- Sys.time()
  ssl_psdh_time  <- ssl_psdh_time2 - ssl_psdh_time1

  ## ----- Comparators: LASSO only (aLASSO / SCAD / MCP commented out) -----
  predictor_names <- design$meta$name
  colnames(x) <- predictor_names

  x_center <- colMeans(x)
  x_scale  <- apply(x, 2L, stats::sd)
  bad_scale <- !is.finite(x_scale) | x_scale <= sqrt(.Machine$double.eps)
  if (any(bad_scale)) {
    stop("Predictors with zero or invalid standard deviation: ",
         paste(predictor_names[bad_scale], collapse = ", "))
  }
  x_std <- sweep(x, 2L, x_center, FUN = "-")
  x_std <- sweep(x_std, 2L, x_scale, FUN = "/")
  colnames(x_std) <- predictor_names

  fit_data <- data.frame(
    TTE    = sim_data$TTE,
    Status = sim_data$Status,
    x_std,
    check.names = FALSE
  )

  fg_formula <- stats::as.formula(
    paste0("Crisk(TTE, Status, cencode = 0, failcode = 1) ~ ",
           paste(predictor_names, collapse = " + "))
  )

  nlambda             <- 100L
  lambda_min_ratio    <- 0.001
  algorithm_tolerance <- 1e-7
  max_iterations      <- 5000L

  fit_fastcrrp <- function(penalty, penalty_factor = NULL, gamma = NULL,
                           lambda = NULL) {
    args <- list(
      formula = fg_formula,
      data = fit_data,
      eps = algorithm_tolerance,
      max.iter = max_iterations,
      standardize = FALSE,
      penalty = penalty,
      lambda.min.ratio = lambda_min_ratio,
      nlambda = nlambda
    )

    if (!is.null(penalty_factor)) {
      if (length(penalty_factor) != npredictors ||
          any(!is.finite(penalty_factor)) ||
          any(penalty_factor <= 0)) {
        stop("penalty_factor must contain one positive finite value per predictor.")
      }
      args$penalty.factor <- as.numeric(penalty_factor)
    }

    if (!is.null(gamma)) args$gamma <- gamma

    if (!is.null(lambda)) {
      if (any(!is.finite(lambda)) || any(lambda < 0)) {
        stop("lambda must contain finite non-negative values.")
      }
      args$lambda <- as.numeric(lambda)
    }

    if ("getBreslowJumps" %in% names(formals(fastcmprsk::fastCrrp))) {
      args$getBreslowJumps <- FALSE
    }

    fit <- do.call(fastcmprsk::fastCrrp, args)
    beta_path <- as.matrix(stats::coef(fit))

    if (nrow(beta_path) != npredictors && ncol(beta_path) == npredictors) {
      beta_path <- t(beta_path)
      fit$coef <- beta_path
    }
    if (nrow(beta_path) != npredictors) {
      stop("Unexpected coefficient-matrix dimensions returned by fastCrrp().")
    }

    rownames(fit$coef) <- predictor_names
    fit
  }

  # LASSO-ONLY: fit_initial_finegray() needed only for aLASSO initial weights
  # fit_initial_finegray <- function() {
  #   args <- list(
  #     formula = fg_formula,
  #     data = fit_data,
  #     eps = algorithm_tolerance,
  #     max.iter = max_iterations,
  #     standardize = FALSE,
  #     variance = FALSE
  #   )
  #   if ("getBreslowJumps" %in% names(formals(fastcmprsk::fastCrr))) {
  #     args$getBreslowJumps <- FALSE
  #   }
  #   do.call(fastcmprsk::fastCrr, args)
  # }

  select_by_bic <- function(fit, method_label) {
    beta_path <- as.matrix(stats::coef(fit))
    if (nrow(beta_path) != npredictors && ncol(beta_path) == npredictors) {
      beta_path <- t(beta_path)
    }

    df_path  <- colSums(beta_path != 0)
    bic_path <- -2 * as.numeric(fit$logLik) + log(nobs) * df_path

    if (length(bic_path) != ncol(beta_path)) {
      stop("BIC path and coefficient path have incompatible lengths.")
    }

    valid <- is.finite(bic_path)
    if (!is.null(fit$converged) && length(fit$converged) == length(bic_path)) {
      valid <- valid & as.logical(fit$converged)
    }
    if (!any(valid)) {
      stop("No finite, converged solution was available for ", method_label, ".")
    }

    bic_for_selection <- bic_path
    bic_for_selection[!valid] <- Inf
    index <- which.min(bic_for_selection)

    beta_std <- as.numeric(beta_path[, index])
    names(beta_std) <- predictor_names

    beta_raw <- beta_std / x_scale
    names(beta_raw) <- predictor_names

    selected <- beta_std != 0

    list(
      method    = method_label,
      index     = index,
      lambda    = as.numeric(fit$lambda.path[index]),
      bic       = as.numeric(bic_path[index]),
      df        = as.integer(df_path[index]),
      converged = if (length(fit$converged) == length(bic_path)) {
        as.logical(fit$converged[index])
      } else {
        NA
      },
      beta_std  = beta_std,
      beta_raw  = beta_raw,
      selected  = selected,
      bic_path  = bic_path,
      df_path   = df_path
    )
  }

  fg_lambda_path <- function(x_std, ftime, status,
                             failcode = 1L, cencode = 0L,
                             nlambda = 100L, lambda_min_ratio = 0.001) {
    n <- length(ftime)
    p <- ncol(x_std)

    cen <- as.integer(status == cencode)
    if (any(cen == 1L)) {
      kmC  <- survival::survfit(survival::Surv(ftime, cen) ~ 1)
      Gfun <- stats::stepfun(kmC$time, c(1, kmC$surv))
    } else {
      Gfun <- function(t) rep(1, length(t))
    }

    ev <- which(status == failcode)
    if (length(ev) == 0L) {
      stop("No cause-", failcode, " events; cannot construct a lambda path.")
    }

    grad <- numeric(p)
    for (i in ev) {
      ti          <- ftime[i]
      at_risk     <- ftime >= ti
      comp_before <- ftime < ti & status != failcode & status != cencode
      w           <- numeric(n)
      w[at_risk]  <- 1
      if (any(comp_before)) {
        w[comp_before] <- Gfun(ti) / Gfun(ftime[comp_before])
      }
      S    <- sum(w)
      xbar <- colSums(w * x_std) / S
      grad <- grad + (x_std[i, ] - xbar)
    }

    lambda_max <- max(abs(grad)) / n
    if (!is.finite(lambda_max) || lambda_max < 1e-8) {
      stop("Degenerate lambda.max (", format(lambda_max, digits = 3),
           "): the null-score projection collapsed. Inspect the outcome/design.")
    }

    10^seq(log10(lambda_max),
           log10(lambda_min_ratio * lambda_max),
           length.out = nlambda)
  }

  other_time1 <- Sys.time()

  reference_lambda <- fg_lambda_path(
    x_std, fit_data$TTE, fit_data$Status,
    failcode = 1L, cencode = 0L,
    nlambda = nlambda, lambda_min_ratio = lambda_min_ratio
  )
  fit_lasso <- fit_fastcrrp(penalty = "LASSO", lambda = reference_lambda)

  # LASSO-ONLY: aLASSO initial-estimate machinery commented out
  # initial_fit       <- NULL
  # initial_beta_std  <- NULL
  # initial_estimator <- NULL
  #
  # if (npredictors < nobs) {
  #   initial_fit <- tryCatch(
  #     fit_initial_finegray(),
  #     error = function(e) {
  #       warning("Unpenalized fastCrr() failed: ", conditionMessage(e))
  #       NULL
  #     }
  #   )
  #
  #   if (!is.null(initial_fit)) {
  #     candidate <- as.numeric(stats::coef(initial_fit))
  #     initial_converged <- is.null(initial_fit$converged) ||
  #       all(as.logical(initial_fit$converged))
  #     if (initial_converged &&
  #         length(candidate) == npredictors &&
  #         all(is.finite(candidate))) {
  #       initial_beta_std  <- candidate
  #       initial_estimator <- "unpenalized fastCrr"
  #     }
  #   }
  # }
  #
  # if (is.null(initial_beta_std)) {
  #   warning("Using a ridge-regularized initial estimate for adaptive-LASSO weights.")
  #   initial_fit <- fit_fastcrrp(penalty = "RIDGE", lambda = reference_lambda)
  #   ridge_beta  <- as.matrix(stats::coef(initial_fit))
  #   ridge_valid <- rep(TRUE, ncol(ridge_beta))
  #   if (!is.null(initial_fit$converged) &&
  #       length(initial_fit$converged) == ncol(ridge_beta)) {
  #     ridge_valid <- as.logical(initial_fit$converged)
  #   }
  #   ridge_valid <- ridge_valid & apply(ridge_beta, 2L, function(z) all(is.finite(z)))
  #   if (!any(ridge_valid)) {
  #     stop("No usable ridge solution was available for adaptive-LASSO weights.")
  #   }
  #   ridge_index <- which.max(replace(as.numeric(initial_fit$logLik), !ridge_valid, -Inf))
  #   initial_beta_std  <- as.numeric(ridge_beta[, ridge_index])
  #   initial_estimator <- "ridge fastCrrp endpoint"
  # }
  #
  # names(initial_beta_std) <- predictor_names
  #
  # adaptive_power    <- 1
  # aLASSO_floor      <- 1e-4
  # aLASSO_weight_cap <- 1e4
  # adaptive_weights  <- 1 / pmax(abs(initial_beta_std), aLASSO_floor)^adaptive_power
  # adaptive_weights  <- pmin(adaptive_weights, aLASSO_weight_cap)
  # adaptive_weights  <- adaptive_weights / min(adaptive_weights)
  # names(adaptive_weights) <- predictor_names
  #
  # fit_alasso <- fit_fastcrrp(
  #   penalty = "LASSO",
  #   penalty_factor = adaptive_weights,
  #   lambda = reference_lambda
  # )

  # LASSO-ONLY: fit_scad and fit_mcp commented out
  # fit_scad <- fit_fastcrrp(penalty = "SCAD", gamma = 3.7, lambda = reference_lambda)
  # fit_mcp  <- fit_fastcrrp(penalty = "MCP",  gamma = 2.7, lambda = reference_lambda)

  fits <- list(LASSO = fit_lasso) # LASSO-ONLY: was list(LASSO, aLASSO, SCAD, MCP)

  same_lambda_grid <- vapply(
    fits,
    function(z) isTRUE(all.equal(z$lambda.path, reference_lambda, tolerance = 1e-12)),
    logical(1)
  )
  if (!all(same_lambda_grid)) {
    stop("The fitted models did not use the same lambda grid: ",
         paste(names(same_lambda_grid)[!same_lambda_grid], collapse = ", "))
  }

  selected_models <- list( # LASSO-ONLY: was list(LASSO, aLASSO, SCAD, MCP)
    LASSO = select_by_bic(fit_lasso, "LASSO")
  )

  other_time2 <- Sys.time()
  other_time  <- other_time2 - other_time1

  list(
    meta = list(nobs = nobs,
                npredictors = npredictors,
                beta1 = unname(beta1),
                beta2 = unname(beta2),
                zero_gap_target = zero_gap_target,
                timestamp = Sys.time(),
                ssl_psdh_time = ssl_psdh_time,
                other_time = other_time,
                block_map = design$meta[, c("name", "block", "rho",
                                            "designed", "active1", "active2")],
                block_size      = design$block_size,
                block_rho       = design$block_rho,
                block_signature = design$signature,
                design_settings = design$settings,
                cens_rate       = cens_rate,
                event_mix       = event_mix),
    sim_data  = sim_data,
    x         = x,
    y         = y,
    mod       = mod,
    autotune  = autotune,
    best      = list(s0 = best_s0, s1 = best_s1),
    final_mod = final_mod,
    fg_true   = fg_true,
    selected_models = selected_models
  )
}


#################### Sweep (this task's runs only) ################
for (npredictors in npredictors_grid) {

  # Per-npredictors subdirectory. The column-to-block map is derived from
  # npredictors, so runs at different p describe different strata and must not be
  # pooled into one summary table. Keeping them in separate directories makes
  # that structural. parse_batch.R recurses, so pointing it at the parent still
  # picks up every p (grouped, as before, by the p in each FILENAME), while
  # pointing it at one p<P> subdir scopes the parse to that scenario alone.
  out_dir_p <- file.path(out_dir, sprintf("p%d", npredictors))
  dir.create(out_dir_p, showWarnings = FALSE, recursive = TRUE)

  for (run in run_start:run_end) {

    # Resolve neg_beta1 sentinel so every run uses -beta1_active as cause-2 effect.
    beta2_run <- if (beta2_active_raw == "neg_beta1") -beta1_active else beta2_active

    tag <- sprintf("n%d_p%d_run%d", nobs, npredictors, run)
    message(sprintf("[%s] %s starting...", format(Sys.time(), "%H:%M:%S"), tag))

    result <- tryCatch(
      run_once(nobs, npredictors, beta1_active, beta2_run,
               active_block, block_props, block_rho, latent_type,
               latent_q, cens_rate, u_min, u_max, zero_gap_target),
      error = function(e) {
        message(sprintf("  !! %s FAILED: %s", tag, conditionMessage(e)))
        list(meta = list(nobs = nobs, npredictors = npredictors, run = run,
                         timestamp = Sys.time(), error = conditionMessage(e)))
      }
    )

    out_file <- file.path(out_dir_p, paste0(tag, ".Rdata"))
    save(result, file = out_file, compress = "xz")
    message(sprintf("  -> saved %s", out_file))
  }
}

message("Batch (task) complete.")
