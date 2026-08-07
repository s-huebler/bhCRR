#################### CHPC Batch Runner ################
# CHPC/SLURM-portable adaptation of Sim/code/batch_run_comparison.R.
#
# All configuration comes from environment variables (set by config.sh and
# passed through by run_array.slurm). Nothing is hardcoded to the Mac.
#
#   REPO_ROOT        path to the repo checkout on CHPC
#   SCRATCH_RUN_DIR  shared-scratch dir where per-run RData is written
#   NOBS             sample size
#   NPREDICTORS      comma list, e.g. "25,50,206"
#   RUN_SUBSTART     first run index this array task handles (from run_array.slurm)
#   RUN_SUBEND       last  run index this array task handles
#   ZERO_GAP_TARGET  clinically-relevant minimum treatment effect
#   BETA1_ACTIVE     comma list of active cause-1 coefficients
#   BETA2_ACTIVE     comma list of active cause-2 coefficients
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

repo_root       <- get_env("REPO_ROOT")
out_dir         <- get_env("SCRATCH_RUN_DIR")
nobs            <- as.integer(get_env("NOBS", "200"))
npredictors_grid<- int_list(get_env("NPREDICTORS", "206"))
run_start       <- as.integer(get_env("RUN_SUBSTART", get_env("RUN_START", "1")))
run_end         <- as.integer(get_env("RUN_SUBEND",   get_env("RUN_END",   "10")))
zero_gap_target <- as.numeric(get_env("ZERO_GAP_TARGET", "0.1"))
beta1_active    <- num_list(get_env("BETA1_ACTIVE", "0.40,-0.50,0.60,0.75,-0.80"))
beta2_active    <- num_list(get_env("BETA2_ACTIVE", "0,0.3,0,0,-0.2"))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message(sprintf("run_batch.R  repo=%s", repo_root))
message(sprintf("  out_dir=%s", out_dir))
message(sprintf("  nobs=%d  npredictors=%s  runs=%d..%d  zero_gap=%.3g",
                nobs, paste(npredictors_grid, collapse = ","),
                run_start, run_end, zero_gap_target))


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

# Start from a clean, non-deterministic RNG state.
if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)


#################### One Run ################
# Verbatim from Sim/code/batch_run_comparison.R -- generates fresh data, fits
# the initial model, autotunes, refits, and returns a bundled result list.
run_once <- function(nobs, npredictors, beta1_active, beta2_active,
                     zero_gap_target) {

  ## ----- Data generation -----
  beta1 <- c(beta1_active, rep(0, npredictors - length(beta1_active)))
  beta2 <- c(beta2_active, rep(0, npredictors - length(beta2_active)))

  sim <- simulateTwoCauseFineGrayModel(nobs = nobs,
                                       beta1 = beta1,
                                       beta2 = beta2,
                                       X = NULL,
                                       u.min = 100,
                                       u.max = 100,
                                       p = 0.5,
                                       returnX = TRUE)

  sim_data <- cbind(data.frame(ID = 1:nobs,
                               TTE = sim$ftime,
                               Status = as.integer(sim$fstatus)),
                    sim$X)

  names(sim_data) <- c("ID", "TTE", "Status", paste0("X_", 1:(ncol(sim_data) - 3)))

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

  ## ----- Other comparators (LASSO / aLASSO / SCAD / MCP) -----
  predictor_names <- paste0("X_", seq_len(npredictors))
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

  select_by_wolbers <- function(fit, method_label, eval_quantile = 0.5){

    beta_path <- as.matrix(stats::coef(fit))
    if (nrow(beta_path) != npredictors && ncol(beta_path) == npredictors) {
      beta_path <- t(beta_path)
    }

    df_path  <- colSums(beta_path != 0)

    eval_time <- as.numeric(stats::quantile(fit_data$TTE[fit_data$Status == 1], eval_quantile))

    X <- as.matrix(fit_data[, predictor_names])
    y_true_mat <- as.matrix(fit_data[, c("TTE", "Status")])

    n_lambda <- ncol(beta_path)
    wolbers_path <- rep(NA, n_lambda)

    dummy_data <- fit_data[, c("TTE", "Status", predictor_names)]
    dummy_data$Event <- as.integer(dummy_data$Status == 1)

    dummy_formula <- stats::as.formula(paste("survival::Surv(TTE, Event) ~", paste(predictor_names, collapse = " + ")))

    for(k in 1:n_lambda) {
      coefs <- as.numeric(beta_path[, k])
      names(coefs) <- predictor_names

      if (any(is.na(coefs))) next

      bh_result <- tryCatch({
        dummy_model <- survival::coxph(
          dummy_formula,
          data = dummy_data,
          init = coefs,
          iter = 0
        )
        survival::basehaz(dummy_model, centered = FALSE)
      }, error = function(e) {
        NULL
      })

      if (is.null(bh_result)) next

      idx <- findInterval(eval_time, bh_result$time)
      base_haz_t <- if (idx == 0) 0 else bh_result$hazard[idx]

      linear_predictors <- X %*% coefs
      risk_score <- as.vector(1 - exp(-base_haz_t * exp(linear_predictors)))

      wolbers_path[k] <- wolbers_c(y_true = y_true_mat, risk_score = risk_score, evaluation_time = eval_time)
    }

    valid <- is.finite(wolbers_path)
    if (!is.null(fit$converged) && length(fit$converged) == length(wolbers_path)) {
      valid <- valid & as.logical(fit$converged)
    }
    if (!any(valid)) {
      stop("No finite, converged solution was available for ", method_label, ".")
    }

    metric_for_selection <- wolbers_path
    metric_for_selection[!valid] <- -Inf
    index <- which.max(metric_for_selection)

    beta_std <- as.numeric(beta_path[, index])
    names(beta_std) <- predictor_names

    beta_raw <- beta_std / x_scale
    names(beta_raw) <- predictor_names

    selected <- beta_std != 0

    list(
      method       = method_label,
      index        = index,
      lambda       = as.numeric(fit$lambda.path[index]),
      wolbers_c    = as.numeric(wolbers_path[index]),
      df           = as.integer(df_path[index]),
      converged    = if (length(fit$converged) == length(wolbers_path)) {
        as.logical(fit$converged[index])
      } else {
        NA
      },
      beta_std     = beta_std,
      beta_raw     = beta_raw,
      selected     = selected,
      wolbers_path = wolbers_path,
      df_path      = df_path
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

  fit_scad <- fit_fastcrrp(penalty = "SCAD", gamma = 3.7, lambda = reference_lambda)
  fit_mcp  <- fit_fastcrrp(penalty = "MCP",  gamma = 2.7, lambda = reference_lambda)

  fits <- list(LASSO = fit_lasso, aLASSO = fit_alasso, SCAD = fit_scad, MCP = fit_mcp)

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
    LASSO_BIC = select_by_bic(fit_lasso, "LASSO_BIC"),
    LASSO_Wolbers = select_by_wolbers(fit_lasso, "LASSO_Wolbers"),
    aLASSO_BIC = select_by_bic(fit_alasso, "aLASSO_BIC"),
    aLASSO_Wolbers = select_by_wolbers(fit_alasso, "aLASSO_Wolbers"),
    SCAD_BIC = select_by_bic(fit_scad, "SCAD_BIC"),
    SCAD_Wolbers = select_by_wolbers(fit_scad, "SCAD_Wolbers"),
    MCP_BIC = select_by_bic(fit_mcp, "MCP_BIC"),
    MCP_Wolbers = select_by_wolbers(fit_mcp, "MCP_Wolbers")
  )

  other_time2 <- Sys.time()
  other_time  <- other_time2 - other_time1

  list(
    meta = list(nobs = nobs,
                npredictors = npredictors,
                beta1 = beta1,
                beta2 = beta2,
                zero_gap_target = zero_gap_target,
                timestamp = Sys.time(),
                ssl_psdh_time = ssl_psdh_time,
                other_time = other_time),
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
  for (run in run_start:run_end) {

    tag <- sprintf("n%d_p%d_run%d", nobs, npredictors, run)
    message(sprintf("[%s] %s starting...", format(Sys.time(), "%H:%M:%S"), tag))

    result <- tryCatch(
      run_once(nobs, npredictors, beta1_active, beta2_active, zero_gap_target),
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

message("Batch (task) complete.")
