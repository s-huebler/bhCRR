#################### Set Up ################
library(tidyverse)
library(fastcmprsk, lib.loc = "/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library")
library(survival)

 source("/Users/sophiehuebler/Documents/bhCRR/R/helpers.r")
 source("/Users/sophiehuebler/Documents/bhCRR/R/utils.r")
 source("/Users/sophiehuebler/Documents/bhCRR/R/dedupe_warnings.r")

 source("/Users/sophiehuebler/Documents/bhCRR/R/update_betas.r")
 source("/Users/sophiehuebler/Documents/bhCRR/R/expected_inclusion_probs.r")
 source("/Users/sophiehuebler/Documents/bhCRR/R/expected_penalty_weights.r")
 source("/Users/sophiehuebler/Documents/bhCRR/R/update_mixture_prob.r")

 Rcpp::sourceCpp("/Users/sophiehuebler/Documents/bhCRR/src/cv_fastcrrp.cpp")
 Rcpp::sourceCpp("/Users/sophiehuebler/Documents/bhCRR/src/RcppExports.cpp")
 source("/Users/sophiehuebler/Documents/bhCRR/R/cv_fastCrrp_cpp.r")
 source("/Users/sophiehuebler/Documents/bhCRR/R/fit_ssl_psdh.r")

 source("/Users/sophiehuebler/Documents/bhCRR/R/generate_foldid.r")
 source("/Users/sophiehuebler/Documents/bhCRR/R/predict_from_ssl_psdh.r")
 source("/Users/sophiehuebler/Documents/bhCRR/R/wolbers_c.r")
 source("/Users/sophiehuebler/Documents/bhCRR/R/cv_ssl_psdh.r")
 source("/Users/sophiehuebler/Documents/bhCRR/R/tune_ssl_psdh.r")
 source("/Users/sophiehuebler/Documents/bhCRR/R/threshold.R")
 source("/Users/sophiehuebler/Documents/bhCRR/R/tuning_diagnostics.r")



#################### Data Generation ################

# To change per scenario
nobs = 200
npredictors = 206
beta1_active = c(0.40, - 0.50, 0.60, 0.75, - 0.80)
beta2_active = c(0, 0.3, 0, 0, -0.2)


# Calculated (fixed)
beta1 <- c(beta1_active, rep(0, npredictors-length(beta1_active)))
beta2 <- c(beta2_active, rep(0, npredictors-length(beta2_active)))

# Simulated
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
predictor_names <- paste0("X_", seq_len(npredictors))
names(sim_data) <- c("ID", "TTE", "Status", predictor_names)

names(sim_data)<- c("ID", "TTE", "Status", paste0("X_", 1:(ncol(sim_data)-3)))

x <- as.matrix(sim_data %>%select(starts_with("X")))
y <- as.matrix(sim_data %>%
                 mutate(Status = as.numeric(Status))%>%
                 select(TTE, Status))

storage.mode(x) <- "double"

if (!all(sim_data$Status %in% c(0L, 1L, 2L))) {
  stop("Status must use 0 = censored, 1 = event of interest, 2 = competing event.")
}
if (sum(sim_data$Status == 1L) == 0L) {
  stop("The simulated data contain no cause-1 events; change the seed or scenario.")
}

# Standardize ONCE and feed this same matrix to every method. This avoids
# method-to-method differences from hidden preprocessing, and it puts the
# adaptive weights on the same coefficient scale as the penalized fits.
x_center <- colMeans(x)
x_scale <- apply(x, 2L, stats::sd)

bad_scale <- !is.finite(x_scale) | x_scale <= sqrt(.Machine$double.eps)
if (any(bad_scale)) {
  stop(
    "Predictors with zero or invalid standard deviation: ",
    paste(predictor_names[bad_scale], collapse = ", ")
  )
}

x_std <- sweep(x, 2L, x_center, FUN = "-")
x_std <- sweep(x_std, 2L, x_scale, FUN = "/")
colnames(x_std) <- predictor_names

fit_data <- data.frame(
  TTE = sim_data$TTE,
  Status = sim_data$Status,
  x_std,
  check.names = FALSE
)

fg_formula <- stats::as.formula(
  paste0(
    "Crisk(TTE, Status, cencode = 0, failcode = 1) ~ ",
    paste(predictor_names, collapse = " + ")
  )
)


#################### FG Model ################
active_idx <- which(beta1 != 0)
active_x <- x[, active_idx, drop = FALSE]

fg_true <- fastcmprsk::fastCrr(fastcmprsk::Crisk(sim$ftime, as.integer(sim$fstatus)) ~ active_x,
                               variance = FALSE)




#################### SSL Model ################

ssl_psdh_time1 <- Sys.time()

#Initial model fit
mod <- fit_ssl_psdh(x, y,
                    ss=c(0.04, 0.6),
                    initial_sparsity = 0.05,
                    theta_a = 1,
                    theta_b = 1,
                    maxit = 50,
                    epsilon=1e-04,
                    init = NULL,
                    init_lam_path = 10^seq(log10(0.1),
                                                log10(0.001),
                                                length = 10),
                    inner_maxit_start = 1000)

#Tuning
#
# The manual pre-flight / post-flight checks are now automated in
# R/tuning_diagnostics.r (bhcrr_autotune). It: estimates theta from the untuned
# fit, clamps the reasonable s1 range into the feasible band, solves for the s0
# region around the clinical zero-gap target, recommends and validates a grid,
# runs tune_ssl_psdh(), and auto-widens + re-tunes once if the optimum lands on
# a grid edge. Everything runs unattended and reports flags for review.

zero_gap_target <- 0.1   # clinically relevant minimum treatment effect

autotune <- bhcrr_autotune(mod,
                           beta_min   = zero_gap_target,
                           beta_floor = 0.01,
                           nfolds     = 10,
                           ncv        = 2,
                           foldid     = NULL,
                           max_widen  = 1)

# Pre-flight diagnostics, the validation checks, and the chosen pair.
# print(autotune$preflight)
# print(autotune)
# print(autotune$validation$checks)
# plot_score_heatmap(autotune$tuning)

# Tuned scale parameters (spike s0, slab s1) selected by the autotuner.
best_s0 <- autotune$best$s0
best_s1 <- autotune$best$s1



final_mod <-fit_ssl_psdh(x, y,
                         ss=c(best_s0, best_s1),
                         initial_sparsity = 0.05,
                         theta_a = 1,
                         theta_b = 1,
                         maxit = 100,
                         epsilon=1e-04,
                         init = NULL,
                         init_lam_path = 10^seq(log10(0.1),
                                                log10(0.001),
                                                length = 10),
                         inner_maxit_start = 1000)



ssl_psdh_time2 <- Sys.time()
ssl_psdh_time = ssl_psdh_time2 - ssl_psdh_time1


#################### Other  Models Setup ################

nlambda <- 100L
lambda_min_ratio <- 0.001
algorithm_tolerance <- 1e-7
max_iterations <- 5000L

# fastcmprsk releases differ slightly in whether fastCrrp() exposes
# getBreslowJumps. Add it only when the installed release supports it.
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

  if (!is.null(gamma)) {
    args$gamma <- gamma
  }

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

  # fastCrrp normally returns predictors by rows and lambda values by columns.
  # Transpose defensively if an installed release returns the opposite layout.
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

  # This is the paper's simplified df choice: model size rather than the
  # effective-df trace. It is also what AIC.fcrrp uses internally.
  df_path <- colSums(beta_path != 0)
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

  # Since x_std = (x - center) / scale, divide by scale to return coefficients
  # to the original predictor units. Centering is absorbed by the baseline
  # subdistribution hazard.
  beta_raw <- beta_std / x_scale
  names(beta_raw) <- predictor_names

  selected <- beta_std != 0

  list(
    method = method_label,
    index = index,
    lambda = as.numeric(fit$lambda.path[index]),
    bic = as.numeric(bic_path[index]),
    df = as.integer(df_path[index]),
    converged = if (length(fit$converged) == length(bic_path)) {
      as.logical(fit$converged[index])
    } else {
      NA
    },
    beta_std = beta_std,
    beta_raw = beta_raw,
    selected = selected,
    bic_path = bic_path,
    df_path = df_path
  )
}


select_by_wolbers <- function(fit, method_label, eval_quantile = 0.5){

  beta_path <- as.matrix(stats::coef(fit))
  if (nrow(beta_path) != npredictors && ncol(beta_path) == npredictors) {
    beta_path <- t(beta_path)
  }

  df_path  <- colSums(beta_path != 0)

  # Calculate evaluation time as median event time for cause 1
  eval_time <- as.numeric(stats::quantile(fit_data$TTE[fit_data$Status == 1], eval_quantile))

  # Prepare predictor matrix and true responses for all patients
  X <- as.matrix(fit_data[, predictor_names])
  y_true_mat <- as.matrix(fit_data[, c("TTE", "Status")])

  n_lambda <- ncol(beta_path)
  wolbers_path <- rep(NA, n_lambda)

  # Setup dummy data once to save time
  dummy_data <- fit_data[, c("TTE", "Status", predictor_names)]
  dummy_data$Event <- as.integer(dummy_data$Status == 1)

  # Build the formula explicitly to avoid the EncodeVars() warning
    dummy_formula <- stats::as.formula(paste("survival::Surv(TTE, Event) ~", paste(predictor_names, collapse = " + ")))

    for(k in 1:n_lambda) {
      coefs <- as.numeric(beta_path[, k])
      names(coefs) <- predictor_names

      # Skip if coefficients are NA
      if (any(is.na(coefs))) next

      # Fit dummy model to extract baseline hazard using predefined coefficients
      # Wrapped in tryCatch: extreme coefficients at small lambdas can cause
      # overflow in exp(X*beta) inside coxph, throwing a .Machine$double.xmax error.
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

    # Calculate absolute risk
    linear_predictors <- X %*% coefs
    risk_score <- as.vector(1 - exp(-base_haz_t * exp(linear_predictors)))

    # Calculate IPCW Wolbers C-index
    wolbers_path[k] <- wolbers_c(y_true = y_true_mat, risk_score = risk_score, evaluation_time = eval_time)
  }

  valid <- is.finite(wolbers_path)
  if (!is.null(fit$converged) && length(fit$converged) == length(wolbers_path)) {
    valid <- valid & as.logical(fit$converged)
  }
  if (!any(valid)) {
    stop("No finite, converged solution was available for ", method_label, ".")
  }

  # Maximize the Wolbers C-index
  metric_for_selection <- wolbers_path
  metric_for_selection[!valid] <- -Inf
  index <- which.max(metric_for_selection)

  beta_std <- as.numeric(beta_path[, index])
  names(beta_std) <- predictor_names

  # x_std = (x - center) / scale, so divide by scale to return coefficients to
  # original predictor units. Centering is absorbed by the baseline hazard.
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

#################### Shared Lambda Grid ################

# fastCrrp()'s built-in lambda grid is generated from a null-score projection,
# max(t(w0 * z) %*% X) / n. When the data contain no censored observations (as
# here: Status is only 1/2), the internal IPCW/score path degenerates and that
# projection collapses to floating-point roundoff (~1e-17), yielding an all-but-
# unpenalized grid on which every comparator selects all coefficients. We instead
# compute lambda.max directly from the KKT condition of the null Fine-Gray model
# and build one clean geometric grid shared by all comparators. The guard makes a
# future collapse fail loudly instead of silently propagating a degenerate path.
fg_lambda_path <- function(x_std, ftime, status,
                           failcode = 1L, cencode = 0L,
                           nlambda = 100L, lambda_min_ratio = 0.001) {
  n <- length(ftime)
  p <- ncol(x_std)

  # Censoring (Kaplan-Meier) survival G(t) for the IPCW weights. With no
  # censored subjects G(t) == 1 and the weights reduce to 1 (exact here).
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

  # Gradient of the null Fine-Gray log-pseudo-likelihood wrt each beta_j:
  #   sum_{i: event} ( x_ij - weighted mean of x_.j over the FG risk set of i ).
  grad <- numeric(p)
  for (i in ev) {
    ti          <- ftime[i]
    at_risk     <- ftime >= ti                                   # weight 1
    comp_before <- ftime < ti & status != failcode & status != cencode
    w           <- numeric(n)
    w[at_risk]  <- 1
    if (any(comp_before)) {
      w[comp_before] <- Gfun(ti) / Gfun(ftime[comp_before])      # IPCW weight
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

#################### LASSO Model ################

lasso_time1 <- Sys.time()

# 4a. LASSO. Build the robust shared lambda grid, then fit on it. Every other
# comparator below already reuses reference_lambda, so all four share this grid.
reference_lambda <- fg_lambda_path(
  x_std, fit_data$TTE, fit_data$Status,
  failcode = 1L, cencode = 0L,
  nlambda = nlambda, lambda_min_ratio = lambda_min_ratio
)
fit_lasso <- fit_fastcrrp(penalty = "LASSO", lambda = reference_lambda)

lasso_time2 <- Sys.time()
lasso_time = lasso_time2 - lasso_time1



d  <- fit_data
# recompute exactly what fastCrrp does internally:
sc <- t(reference_lambda)  # placeholder; instead inspect directly:
range(fit_lasso$lambda.path)                 # ~1e-17 .. 1e-20  -> collapsed
c(sd_TTE = sd(d$TTE),
  n_unique_time = length(unique(d$TTE)),
  n_unique_evt_time = length(unique(d$TTE[d$Status==1])),
  table_status = NA); print(table(d$Status))

#################### aLASSO Model ################
alasso_time1 <- Sys.time()

# 4b. Adaptive LASSO
# The paper uses theta_j = 1 / abs(beta_hat_j), where beta_hat is the
# unpenalized Fine-Gray estimate. With p < n, use fastCrr(). If that fit is
# unavailable (e.g., in a high-dimensional extension), use a ridge endpoint.
initial_fit <- NULL
initial_beta_std <- NULL
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
      initial_beta_std <- candidate
      initial_estimator <- "unpenalized fastCrr"
    }
  }
}

if (is.null(initial_beta_std)) {
  warning(
    "Using a ridge-regularized initial estimate for adaptive-LASSO weights."
  )
  initial_fit <- fit_fastcrrp(
    penalty = "RIDGE",
    lambda = reference_lambda
  )
  ridge_beta <- as.matrix(stats::coef(initial_fit))
  ridge_valid <- rep(TRUE, ncol(ridge_beta))
  if (!is.null(initial_fit$converged) &&
      length(initial_fit$converged) == ncol(ridge_beta)) {
    ridge_valid <- as.logical(initial_fit$converged)
  }
  ridge_valid <- ridge_valid & apply(ridge_beta, 2L, function(z) all(is.finite(z)))
  if (!any(ridge_valid)) {
    stop("No usable ridge solution was available for adaptive-LASSO weights.")
  }
  # Among converged ridge solutions, use the one with the largest log-likelihood
  # (typically the least-penalized endpoint).
  ridge_index <- which.max(replace(
    as.numeric(initial_fit$logLik),
    !ridge_valid,
    -Inf
  ))
  initial_beta_std <- as.numeric(ridge_beta[, ridge_index])
  initial_estimator <- "ridge fastCrrp endpoint"
}

names(initial_beta_std) <- predictor_names

# Numerical stabilization of 1 / abs(beta_hat): floor avoids infinite weights,
# cap avoids extreme conditioning. Dividing by the smallest weight preserves
# all relative weights while making every penalty factor at least one; this helps
# the shared lambda grid start at a sparse/all-zero adaptive-LASSO solution.
adaptive_power <- 1
aLASSO_floor <- 1e-4
aLASSO_weight_cap <- 1e4
adaptive_weights <- 1 / pmax(abs(initial_beta_std), aLASSO_floor)^adaptive_power
adaptive_weights <- pmin(adaptive_weights, aLASSO_weight_cap)
adaptive_weights <- adaptive_weights / min(adaptive_weights)
names(adaptive_weights) <- predictor_names

fit_alasso <- fit_fastcrrp(
  penalty = "LASSO",
  penalty_factor = adaptive_weights,
  lambda = reference_lambda
)

alasso_time2 <- Sys.time()
alasso_time = alasso_time2 - alasso_time1

#################### SCAD Model ################

scad_time1 <- Sys.time()
# 4c. SCAD. Explicit gamma = 3.7 is the paper's fixed concavity parameter.
fit_scad <- fit_fastcrrp(
  penalty = "SCAD",
  gamma = 3.7,
  lambda = reference_lambda
)

scad_time2 <- Sys.time()
scad_time = scad_time2 - scad_time1

#################### MCP Model ################

mcp_time1 <- Sys.time()

# 4d. MCP. Explicit gamma = 2.7 is the paper's fixed concavity parameter.
fit_mcp <- fit_fastcrrp(
  penalty = "MCP",
  gamma = 2.7,
  lambda = reference_lambda
)

mcp_time2 <- Sys.time()
mcp_time = mcp_time2 - mcp_time1

#################### Fits Model ################

fits <- list(
  LASSO = fit_lasso,
  aLASSO = fit_alasso,
  SCAD = fit_scad,
  MCP = fit_mcp
)

# The LASSO-generated reference grid was passed explicitly to the other fits.
# Verify equality so comparisons cannot silently use different tuning grids.
same_lambda_grid <- vapply(
  fits,
  function(z) isTRUE(all.equal(z$lambda.path, reference_lambda, tolerance = 1e-12)),
  logical(1)
)
if (!all(same_lambda_grid)) {
  stop(
    "The fitted models did not use the same lambda grid: ",
    paste(names(same_lambda_grid)[!same_lambda_grid], collapse = ", ")
  )
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


mod_comparison <- bind_cols(
  final_mod$coefficients,
  beta1,
  c(fg_true$coef, rep(NA, npredictors-length(beta1_active))),
  selected_models$LASSO_BIC$beta_raw,
  selected_models$LASSO_Wolbers$beta_raw,
  selected_models$aLASSO_BIC$beta_raw,
  selected_models$aLASSO_Wolbers$beta_raw,
  selected_models$SCAD_BIC$beta_raw,
  selected_models$SCAD_Wolbers$beta_raw,
  selected_models$MCP_BIC$beta_raw,
  selected_models$MCP_Wolbers$beta_raw
)

names(mod_comparison)<- c("predictor", "SSL", "True", "FG", "LASSO_BIC", "LASSO_Wolbers", "aLASSO_BIC", "aLASSO_Wolbers", "SCAD_BIC", "SCAD_Wolbers", "MCP_BIC", "MCP_Wolbers")

mod_comparison <- mod_comparison %>%
  select(c("predictor", "True", "FG", "SSL", "LASSO_BIC", "LASSO_Wolbers", "aLASSO_BIC", "aLASSO_Wolbers", "SCAD_BIC", "SCAD_Wolbers", "MCP_BIC", "MCP_Wolbers"))

mod_comparison_long <- mod_comparison %>%
  pivot_longer(c("SSL", "LASSO_BIC", "LASSO_Wolbers", "aLASSO_BIC", "aLASSO_Wolbers", "SCAD_BIC", "SCAD_Wolbers", "MCP_BIC", "MCP_Wolbers"),
               names_to = "Model",
               values_to = "Estimate")%>%
  mutate(Bias = Estimate - True)

# Which lambdas actually converged, and is BIC/Wolbers being forced to a grid edge?
diag_one <- function(fit, sel, nm) {
  # Dynamically handle whether this model used Wolbers (max) or BIC (min)
  if ("wolbers_path" %in% names(sel)) {
    metric_path <- sel$wolbers_path
    idx <- which.max(replace(metric_path, !is.finite(metric_path), -Inf))
  } else {
    metric_path <- sel$bic_path
    idx <- which.min(replace(metric_path, !is.finite(metric_path), Inf))
  }

  data.frame(
    method       = nm,
    n_lambda     = length(metric_path),
    n_converged  = sum(as.logical(fit$converged)),
    metric_idx   = idx,
    at_grid_edge = idx %in% c(1L, length(metric_path)),
    df_at_sel    = sel$df,
    df_range     = paste0(min(sel$df_path), "-", max(sel$df_path)),
    lambda_sel   = signif(sel$lambda, 3)
  )
}

# Pass the correct fit model corresponding to each selection metric
fit_list <- list(
  fit_lasso, fit_lasso,
  fit_alasso, fit_alasso,
  fit_scad, fit_scad,
  fit_mcp, fit_mcp
)

do.call(rbind, Map(diag_one,
                   fit_list,
                   selected_models,
                   names(selected_models)))

# aLASSO grid mis-scaling: how many leading all-zero (null) columns each path has.
sapply(selected_models, function(s) which(s$df_path > 0)[1] - 1L)

