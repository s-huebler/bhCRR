#################### Batch Runner ################
# Runs the single_run workflow repeatedly across a set of npredictors values,
# saving one result file per run to Sim/Local_Testing/.
#
# NOTE: set.seed() is intentionally NOT called anywhere. Each call to
# simulateTwoCauseFineGrayModel() draws fresh data, so every run produces a
# different dataset and therefore a different fitted model.
#
# The predictor VALUES therefore change from run to run, but the column-to-block
# assignment does not: generate_block_X() derives it from npredictors,
# block_props, block_rho and active_block with no RNG involved. X_17 is in the
# same block in run 1 and run 100, which is what makes per-block pooling in
# parse_batch_run.R valid. meta$block_signature records the map so a mismatch
# across pooled runs is caught rather than averaged over.
#
# Files are named like: n200_p25_run1.Rdata
# Each file contains a single object `result` (a list); load with:
#   load("Sim/Local_Testing/n200_p25_run1.Rdata"); str(result, max.level = 1)


#################### Set Up (load everything once) ################
library(tidyverse)
library(fastcmprsk, lib.loc = "/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library")
library(survival)

repo_root <- "/Users/sophiehuebler/Documents/bhCRR"

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

# Block-correlated design matrix generator, shared with single_run_complex.R so
# the two scripts cannot drift apart.
source(file.path(repo_root, "Sim/code/generate_block_X.R"))

# Start from a clean, non-deterministic RNG state. Remove any existing global
# seed so each run draws fresh data (belt-and-suspenders with the seed-leak
# fixes now in the CV/sim functions).
if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)


#################### Batch Configuration ################
# ---- edit these to control the sweep ----
nobs             <- 200                 # sample size (fixed across runs)
npredictors_grid <- c(221)      # set of npredictors to sweep over
run_start        <- 1           # e.g., set to 1 for first batch, 11 for second
run_end          <- 10           # e.g., set to 10 for first batch, 20 for second

# Effects of the designed active predictors. beta1_active, beta2_active and
# active_block are aligned elementwise: entry k is one predictor, with its
# cause-1 effect, its cause-2 effect, and the block it is planted in.
beta1_active <- c(0.40, -0.50, 0.60, 0.75, -0.80)
#beta2_active <- c(0,     0.3,  0,    0,   -0.2)
beta2_active <- -beta1_active
active_block <- c("strong", "strong", "weak", "weak", "indep")

# Predictor correlation blocks. Must match single_run_complex.R for the two to
# be comparable. The column-to-block map is deterministic given these settings
# plus npredictors, so a given predictor sits in the same block in every run of
# a sweep; meta$block_signature records the map so parse_batch_run.R can refuse
# to pool runs whose maps differ.
block_props <- c(indep = 0.50, weak = 0.25, strong = 0.25)
block_rho   <- c(indep = 0.00, weak = 0.30, strong = 0.60)
latent_type <- "gaussian"
latent_q    <- c(indep = 0.5, weak = 0.4, strong = 0.5)

# Random exponential censoring, ported from single_run_complex.R so the batch
# generates the same kind of data. Set cens_rate = 0 for no random censoring.
cens_rate <- 0.05
u_min <- 100
u_max <- 100

zero_gap_target <- 0.1   # clinically relevant minimum treatment effect

out_dir <- file.path(repo_root, "Sim/Local_Testing")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


#################### One Run ################
# Generates fresh data, fits the initial model, autotunes, refits, and returns
# a bundled result list. No set.seed() here on purpose.
run_once <- function(nobs, npredictors, beta1_active, beta2_active,
                     active_block, block_props, block_rho, latent_type,
                     latent_q, cens_rate, u_min, u_max,
                     zero_gap_target) {

  ## ----- Data generation -----
  # Block-correlated design matrix. The column-to-block map uses no RNG, so it
  # is identical across every run of this sweep at a given npredictors.
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

  beta1 <- design$beta1
  beta2 <- design$beta2
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

  # Random exponential censoring: C_i ~ Exp(cens_rate), censored where C_i lands
  # before the simulated event time. Overwrite sim$ftime / sim$fstatus so the
  # oracle Fine-Gray fit below sees the censored data too.
  if (!is.finite(cens_rate) || cens_rate < 0) {
    stop("cens_rate must be a finite, non-negative Exponential rate.")
  }
  if (cens_rate > 0) {
    cens_time <- stats::rexp(nobs, rate = cens_rate)
    is_censored <- cens_time < sim$ftime
    sim$ftime[is_censored] <- cens_time[is_censored]
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

  # Initial model fit
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

  # Tuning
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

  ## ----- Other comparators (LASSO / aLASSO / SCAD / MCP) -----
  # Standardize ONCE and feed this same matrix to every comparator, mirroring
  # single_run.R. This puts the adaptive weights on the same coefficient scale
  # as the penalized fits and avoids method-to-method preprocessing differences.
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

  # fastcmprsk releases differ in whether fastCrrp() exposes getBreslowJumps.
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

  fit_initial_finegray <- function() {
    args <- list(
      formula = fg_formula,
      data = fit_data,
      eps = algorithm_tolerance,
      max.iter = max_iterations,
      standardize = FALSE,
      variance = FALSE
    )
    if ("getBreslowJumps" %in% names(formals(fastcmprsk::fastCrr))) {
      args$getBreslowJumps <- FALSE
    }
    do.call(fastcmprsk::fastCrr, args)
  }

  select_by_bic <- function(fit, method_label) {
    beta_path <- as.matrix(stats::coef(fit))
    if (nrow(beta_path) != npredictors && ncol(beta_path) == npredictors) {
      beta_path <- t(beta_path)
    }

    # Paper's simplified df choice: model size rather than the effective-df
    # trace (also what AIC.fcrrp uses internally).
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

    # x_std = (x - center) / scale, so divide by scale to return coefficients to
    # original predictor units. Centering is absorbed by the baseline hazard.
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

  # select_by_wolbers <- function(fit, method_label, eval_quantile = 0.5){
  #
  #   beta_path <- as.matrix(stats::coef(fit))
  #   if (nrow(beta_path) != npredictors && ncol(beta_path) == npredictors) {
  #     beta_path <- t(beta_path)
  #   }
  #
  #   # Paper's simplified df choice: model size rather than the effective-df
  #   # trace (also what AIC.fcrrp uses internally).
  #   df_path  <- colSums(beta_path != 0)
  #   bic_path <- -2 * as.numeric(fit$logLik) + log(nobs) * df_path
  #
  #   if (length(bic_path) != ncol(beta_path)) {
  #     stop("BIC path and coefficient path have incompatible lengths.")
  #   }
  #
  #   valid <- is.finite(bic_path)
  #   if (!is.null(fit$converged) && length(fit$converged) == length(bic_path)) {
  #     valid <- valid & as.logical(fit$converged)
  #   }
  #   if (!any(valid)) {
  #     stop("No finite, converged solution was available for ", method_label, ".")
  #   }
  #
  #   bic_for_selection <- bic_path
  #   bic_for_selection[!valid] <- Inf
  #   #index <- which.min(bic_for_selection)
  #
  #   original_data <- fit$call$data
  #
  #   dummy_model <- coxph(Surv(TTE, Status) ~ .,
  #                        data = original_data)
  #   eval_time <- as.numeric(quantile(y[y[,2] == 1, 1], eval_quantile))
  #
  #   bh <- basehaz(dummy_model, centered = F)
  #
  #   idx <- findInterval(eval_time, bh$time)
  #
  #   ###INSERT THE PREDICTED RISK HERE
  #
  #   beta_std <- as.numeric(beta_path[, index])
  #   names(beta_std) <- predictor_names
  #
  #   # x_std = (x - center) / scale, so divide by scale to return coefficients to
  #   # original predictor units. Centering is absorbed by the baseline hazard.
  #   beta_raw <- beta_std / x_scale
  #   names(beta_raw) <- predictor_names
  #
  #   selected <- beta_std != 0
  #
  #   list(
  #     method    = method_label,
  #     index     = index,
  #     lambda    = as.numeric(fit$lambda.path[index]),
  #     bic       = as.numeric(bic_path[index]),
  #     df        = as.integer(df_path[index]),
  #     converged = if (length(fit$converged) == length(bic_path)) {
  #       as.logical(fit$converged[index])
  #     } else {
  #       NA
  #     },
  #     beta_std  = beta_std,
  #     beta_raw  = beta_raw,
  #     selected  = selected,
  #     bic_path  = bic_path,
  #     df_path   = df_path
  #   )
  # }

  # fastCrrp()'s built-in lambda grid degenerates when there are no censored
  # observations (the IPCW/score projection collapses to roundoff). Build one
  # clean geometric grid from the KKT condition of the null Fine-Gray model and
  # share it across all comparators. See single_run.R for the full rationale.
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

  # LASSO: build the robust shared lambda grid, then fit on it.
  reference_lambda <- fg_lambda_path(
    x_std, fit_data$TTE, fit_data$Status,
    failcode = 1L, cencode = 0L,
    nlambda = nlambda, lambda_min_ratio = lambda_min_ratio
  )
  fit_lasso <- fit_fastcrrp(penalty = "LASSO", lambda = reference_lambda)

  # Adaptive LASSO: theta_j = 1 / abs(beta_hat_j) from the unpenalized FG fit
  # (p < n), falling back to a ridge endpoint if that fit is unavailable.
  initial_fit       <- NULL
  initial_beta_std  <- NULL
  initial_estimator <- NULL

  if (npredictors < nobs) {
    initial_fit <- tryCatch(
      fit_initial_finegray(),
      error = function(e) {
        warning("Unpenalized fastCrr() failed: ", conditionMessage(e))
        NULL
      }
    )

    if (!is.null(initial_fit)) {
      candidate <- as.numeric(stats::coef(initial_fit))
      initial_converged <- is.null(initial_fit$converged) ||
        all(as.logical(initial_fit$converged))
      if (initial_converged &&
          length(candidate) == npredictors &&
          all(is.finite(candidate))) {
        initial_beta_std  <- candidate
        initial_estimator <- "unpenalized fastCrr"
      }
    }
  }

  if (is.null(initial_beta_std)) {
    warning("Using a ridge-regularized initial estimate for adaptive-LASSO weights.")
    initial_fit <- fit_fastcrrp(penalty = "RIDGE", lambda = reference_lambda)
    ridge_beta  <- as.matrix(stats::coef(initial_fit))
    ridge_valid <- rep(TRUE, ncol(ridge_beta))
    if (!is.null(initial_fit$converged) &&
        length(initial_fit$converged) == ncol(ridge_beta)) {
      ridge_valid <- as.logical(initial_fit$converged)
    }
    ridge_valid <- ridge_valid & apply(ridge_beta, 2L, function(z) all(is.finite(z)))
    if (!any(ridge_valid)) {
      stop("No usable ridge solution was available for adaptive-LASSO weights.")
    }
    ridge_index <- which.max(replace(as.numeric(initial_fit$logLik), !ridge_valid, -Inf))
    initial_beta_std  <- as.numeric(ridge_beta[, ridge_index])
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

  fit_alasso <- fit_fastcrrp(
    penalty = "LASSO",
    penalty_factor = adaptive_weights,
    lambda = reference_lambda
  )

  # SCAD (gamma = 3.7) and MCP (gamma = 2.7): the paper's fixed concavity params.
  fit_scad <- fit_fastcrrp(penalty = "SCAD", gamma = 3.7, lambda = reference_lambda)
  fit_mcp  <- fit_fastcrrp(penalty = "MCP",  gamma = 2.7, lambda = reference_lambda)

  fits <- list(LASSO = fit_lasso, aLASSO = fit_alasso, SCAD = fit_scad, MCP = fit_mcp)

  # All comparators were handed the same reference grid explicitly; verify so
  # comparisons cannot silently use different tuning grids.
  same_lambda_grid <- vapply(
    fits,
    function(z) isTRUE(all.equal(z$lambda.path, reference_lambda, tolerance = 1e-12)),
    logical(1)
  )
  if (!all(same_lambda_grid)) {
    stop("The fitted models did not use the same lambda grid: ",
         paste(names(same_lambda_grid)[!same_lambda_grid], collapse = ", "))
  }

  selected_models <- list(
    LASSO  = select_by_bic(fit_lasso,  "LASSO"),
    aLASSO = select_by_bic(fit_alasso, "adaptive LASSO"),
    SCAD   = select_by_bic(fit_scad,   "SCAD"),
    MCP    = select_by_bic(fit_mcp,    "MCP")
  )

  other_time2 <- Sys.time()
  other_time  <- other_time2 - other_time1

  ## ----- Bundle everything needed for downstream diagnostics -----
  list(
    meta = list(nobs = nobs,
                npredictors = npredictors,
                beta1 = unname(beta1),
                beta2 = unname(beta2),
                zero_gap_target = zero_gap_target,
                timestamp = Sys.time(),
                ssl_psdh_time = ssl_psdh_time,
                other_time = other_time,
                # Per-predictor correlation stratum. Saved with the run rather
                # than recomputed downstream, so a later edit to block_props
                # cannot silently relabel results already on disk. Joined by
                # name in parse_batch_run.R, never by position.
                block_map = design$meta[, c("name", "block", "rho",
                                            "designed", "active1", "active2")],
                block_size = design$block_size,
                block_rho = design$block_rho,
                # Fingerprint of the column-to-block map. Runs pooled into one
                # table must share this string, or their blocks mean different
                # things and per-block summaries are meaningless.
                block_signature = design$signature,
                design_settings = design$settings,
                cens_rate = cens_rate,
                event_mix = event_mix),
    sim_data  = sim_data,
    x         = x,
    y         = y,
    mod       = mod,
    autotune  = autotune,
    best      = list(s0 = best_s0, s1 = best_s1),
    final_mod = final_mod,
    fg_true   = fg_true,
    # BIC-selected comparators (LASSO / aLASSO / SCAD / MCP). Each element is the
    # select_by_bic() summary, incl. beta_raw (original units) and beta_std.
    selected_models = selected_models
  )
}


#################### Sweep ################
for (npredictors in npredictors_grid) {

  # Inner loop now iterates over your specific run sequence
  for (run in run_start:run_end) {

    tag <- sprintf("n%d_p%d_run%d", nobs, npredictors, run)
    message(sprintf("[%s] %s starting...", format(Sys.time(), "%H:%M:%S"), tag))

    result <- tryCatch(
      run_once(nobs, npredictors, beta1_active, beta2_active,
               active_block, block_props, block_rho, latent_type,
               latent_q, cens_rate, u_min, u_max, zero_gap_target),
      error = function(e) {
        message(sprintf("  !! %s FAILED: %s", tag, conditionMessage(e)))
        list(meta = list(nobs = nobs, npredictors = npredictors, run = run,
                         timestamp = Sys.time(), error = conditionMessage(e)))
      }
    )

    out_file <- file.path(out_dir, paste0(tag, ".Rdata"))
    save(result, file = out_file, compress = "xz")
    message(sprintf("  -> saved %s", out_file))
  }
}


message("Batch complete.")
