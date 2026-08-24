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

 # Block-correlated design matrix generator (Binder & Schumacher 2008 sec 3.1;
 # Binder et al. 2009 sec 3.2.1). Kept in its own file so batch_run_complex.R
 # uses the identical generator instead of a drifting copy.
 source("/Users/sophiehuebler/Documents/bhCRR/Sim/code/generate_block_X.R")



#################### Data Generation ################

# To change per scenario
nobs = 200
npredictors = 25

# Effects of the designed active predictors. beta1_active, beta2_active and
# active_block are aligned elementwise: entry k is one predictor, with its
# cause-1 effect, its cause-2 effect, and the block it is planted in. An entry
# with one beta set to 0 is still a designed active, which is how Binder et al.
# (2009) place covariates acting on a single cause.
beta1_active = c(0.40, - 0.50, 0.60, 0.75, - 0.80)
beta2_active = c(0, 0.3, 0, 0, -0.2)
active_block = c("strong", "strong", "weak", "weak", "indep")

# Predictor correlation blocks. block_props are relative shares of npredictors
# (apportioned by largest remainder so they sum exactly); block_rho is the
# within-block correlation, which is exchangeable inside a block and exactly
# zero between blocks. Names must match across all three vectors.
block_props = c(indep = 0.50, weak = 0.25, strong = 0.25)
block_rho   = c(indep = 0.00, weak = 0.30, strong = 0.60)

# Latent driving the within-block correlation. "gaussian" hits block_rho exactly
# in every sample. "bernoulli" is the papers' own construction, giving bimodal
# microarray-like marginals, and then latent_q sets the shift probability per
# block. To replicate Binder et al. (2009) directly, use latent = "bernoulli"
# with block_rho = c(indep = 0, weak = 0.05, strong = 0.5) and
# latent_q = c(indep = 0.5, weak = 0.7, strong = 0.5).
latent_type = "gaussian"
latent_q    = c(indep = 0.5, weak = 0.4, strong = 0.5)

# Random censoring. Censoring times are drawn as C_i ~ Exponential(cens_rate),
# independently of the predictors and of the event process. A larger rate means
# censoring happens earlier, so a larger rate gives a higher censored fraction.
# The mean censoring time is 1 / cens_rate, so pick the rate relative to the
# scale of the event times (roughly, rate ~ 1 / (typical follow-up time)).
# Set cens_rate = 0 to turn random censoring off entirely.
cens_rate = 0.05

# Administrative censoring bounds handed to simulateTwoCauseFineGrayModel().
# Keeping u_min = u_max = 100 leaves the generator's own uniform censoring
# inactive (no event time reaches 100), so cens_rate is the only censoring
# mechanism. Lower these if administrative censoring is wanted as well.
u_min = 100
u_max = 100


# Calculated (fixed)
# generate_block_X() builds the block-correlated design matrix, places each
# designed active in the block named for it, and scatters the betas to the
# matching columns. With layout = "active_first" the actives occupy columns
# 1..length(beta1_active) in the order written above, so beta1[k] is the effect
# of beta1_active[k]. design$meta carries the per-predictor block labels.
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

# Realized block structure for this replicate. Worth eyeballing when a new
# block specification is used for the first time.
# print(design$block_size)
# print(block_corr_summary(X_design, design$block_id, design$block_rho))

# Simulated
sim <- simulateTwoCauseFineGrayModel(nobs = nobs,
                                     beta1 = beta1,
                                     beta2 = beta2,
                                     X = X_design,
                                     u.min = u_min,
                                     u.max = u_max,
                                     p = 0.5,
                                     returnX = TRUE)

# The event times are only interpretable if the generator used the design matrix
# we handed it rather than one of its own. This catches a silent substitution.
# It cannot catch internal rescaling of X prior to generating the event times,
# which is worth confirming once against the fastcmprsk source.
if (!isTRUE(all.equal(unname(as.matrix(sim$X)), unname(X_design),
                      tolerance = 1e-12))) {
  stop("simulateTwoCauseFineGrayModel() did not return the supplied design ",
       "matrix unchanged; the event times were not generated from X_design.")
}

# Apply the exponential random censoring on top of the generated event times.
# Subject i is censored when C_i lands before the event time already simulated
# for them, in which case the observed time becomes C_i and Status becomes 0.
# Overwriting sim$ftime / sim$fstatus here (rather than only sim_data) keeps
# every downstream consumer of `sim`, including the oracle Fine-Gray fit, on the
# same censored data.
if (!is.finite(cens_rate) || cens_rate < 0) {
  stop("cens_rate must be a finite, non-negative Exponential rate.")
}

if (cens_rate > 0) {
  cens_time <- stats::rexp(nobs, rate = cens_rate)
  is_censored <- cens_time < sim$ftime
  sim$ftime[is_censored] <- cens_time[is_censored]
  sim$fstatus[is_censored] <- 0
}

# Realized outcome mix for this replicate. Worth logging per scenario: the
# censored fraction is a consequence of cens_rate and the event-time scale, not
# something set directly.
cens_prop <- mean(sim$fstatus == 0)
event_mix <- c(censored = cens_prop,
               cause1 = mean(sim$fstatus == 1),
               cause2 = mean(sim$fstatus == 2))

sim_data <- cbind(data.frame(ID = 1:nobs,
                             TTE = sim$ftime,
                             Status = as.integer(sim$fstatus)),
                  X_design)
predictor_names <- design$meta$name
names(sim_data) <- c("ID", "TTE", "Status", predictor_names)

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
if (cens_rate > 0 && cens_prop > 0.9) {
  warning("Over 90% of observations are censored at cens_rate = ", cens_rate,
          "; the scenario may be uninformative.")
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

#################### Shared Lambda Grid ################

# fastCrrp()'s built-in lambda grid is generated from a null-score projection,
# max(t(w0 * z) %*% X) / n. When the data contain no censored observations (the
# cens_rate = 0 case, where Status is only 1/2), the internal IPCW/score path
# degenerates and that projection collapses to roundoff (~1e-17), yielding an all-but-
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
  LASSO = select_by_bic(fit_lasso, "LASSO"),
  aLASSO = select_by_bic(fit_alasso, "adaptive LASSO"),
  SCAD = select_by_bic(fit_scad, "SCAD"),
  MCP = select_by_bic(fit_mcp, "MCP")
)



# The oracle Fine-Gray fit only estimates the cause-1 active columns, so scatter
# its coefficients back to their positions by index. Building this column as
# c(fg_true$coef, rep(NA, ...)) assumed the cause-1 actives were the leading
# columns of x, which stops being true as soon as a beta1_active entry is 0.
fg_col <- rep(NA_real_, npredictors)
fg_col[active_idx] <- as.numeric(fg_true$coef)

mod_comparison <- bind_cols(final_mod$coefficients,beta1, fg_col, selected_models$LASSO$beta_raw, selected_models$aLASSO$beta_raw, selected_models$SCAD$beta_raw, selected_models$MCP$beta_raw, )
names(mod_comparison)<- c("predictor", "SSL", "True", "FG", "LASSO_BIC", "aLASSO_BIC", "SCAD_BIC", "MCP_BIC")
mod_comparison <- mod_comparison %>%
  select(c("predictor", "True", "FG", "SSL",  "LASSO_BIC", "aLASSO_BIC", "SCAD_BIC", "MCP_BIC")) %>%
  # Correlation stratum for each predictor, so bias and selection can be read
  # by block. design$meta is in column order, matching mod_comparison's rows.
  mutate(Block = design$meta$block, .after = predictor)

mod_comparison_long <- mod_comparison %>%
  pivot_longer(c("SSL",  "LASSO_BIC", "aLASSO_BIC", "SCAD_BIC", "MCP_BIC"),
               names_to = "Model",
               values_to = "Estimate")%>%
  mutate(Bias = Estimate - True)

plot_model_comparison(mod_comparison_long)

