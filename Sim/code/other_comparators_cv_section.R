  ## ----- Other comparators (LASSO / aLASSO / SCAD / MCP) -----
  # Standardize ONCE and feed this same matrix to every comparator, mirroring
  # single_run.R. This puts the adaptive weights on the same coefficient scale
  # as the penalized fits and avoids method-to-method preprocessing differences.
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

  # Match the SSL-PSDH tuning structure: repeated k-fold CV, a single evaluation
  # horizon, and pooled out-of-fold Wolbers C-index within each CV repetition.
  comparator_nfolds       <- 10L
  comparator_ncv          <- 2L
  comparator_eval_quantile <- 0.5

  # fastcmprsk releases differ in whether fastCrrp() exposes getBreslowJumps.
  # `data` is an argument so the same fitting routine can be reused in each
  # training fold without changing the full-data fit.
  fit_fastcrrp <- function(penalty, penalty_factor = NULL, gamma = NULL,
                           lambda = NULL, data = fit_data,
                           get_breslow_jumps = FALSE) {
    args <- list(
      formula = fg_formula,
      data = data,
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
      args$getBreslowJumps <- get_breslow_jumps
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

  fit_initial_finegray <- function(data = fit_data) {
    args <- list(
      formula = fg_formula,
      data = data,
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

  adaptive_power    <- 1
  aLASSO_floor      <- 1e-4
  aLASSO_weight_cap <- 1e4

  # Derive adaptive-LASSO weights from the data supplied to this helper. During
  # CV this means the weights are estimated from the training fold only, which
  # avoids outcome leakage from the held-out fold. The full-data final fit still
  # uses weights estimated from the full data, as before.
  get_adaptive_lasso_weights <- function(data, lambda, warn_on_ridge = FALSE) {
    n_train          <- nrow(data)
    initial_fit      <- NULL
    initial_beta_std <- NULL
    estimator        <- NULL

    if (npredictors < n_train) {
      initial_fit <- tryCatch(
        fit_initial_finegray(data = data),
        error = function(e) {
          if (warn_on_ridge) {
            warning("Unpenalized fastCrr() failed: ", conditionMessage(e))
          }
          NULL
        }
      )

      if (!is.null(initial_fit)) {
        candidate <- as.numeric(stats::coef(initial_fit))
        initial_converged <- TRUE
        if (!is.null(initial_fit$converged)) {
          initial_status <- as.logical(initial_fit$converged)
          initial_converged <- all(!is.na(initial_status) & initial_status)
        }
        if (initial_converged &&
            length(candidate) == npredictors &&
            all(is.finite(candidate))) {
          initial_beta_std <- candidate
          estimator <- "unpenalized fastCrr"
        }
      }
    }

    if (is.null(initial_beta_std)) {
      if (warn_on_ridge) {
        warning("Using a ridge-regularized initial estimate for adaptive-LASSO weights.")
      }

      initial_fit <- fit_fastcrrp(
        penalty = "RIDGE",
        lambda = lambda,
        data = data
      )

      ridge_beta <- as.matrix(stats::coef(initial_fit))
      ridge_loglik <- as.numeric(initial_fit$logLik)
      if (length(ridge_loglik) != ncol(ridge_beta)) {
        stop("Ridge log-likelihood and coefficient paths have incompatible lengths.")
      }

      ridge_valid <- apply(ridge_beta, 2L, function(z) all(is.finite(z))) &
        is.finite(ridge_loglik)
      if (!is.null(initial_fit$converged) &&
          length(initial_fit$converged) == ncol(ridge_beta)) {
        ridge_converged <- as.logical(initial_fit$converged)
        ridge_converged[is.na(ridge_converged)] <- FALSE
        ridge_valid <- ridge_valid & ridge_converged
      }
      if (!any(ridge_valid)) {
        stop("No usable ridge solution was available for adaptive-LASSO weights.")
      }

      ridge_for_selection <- ridge_loglik
      ridge_for_selection[!ridge_valid] <- -Inf
      ridge_index <- which.max(ridge_for_selection)
      initial_beta_std <- as.numeric(ridge_beta[, ridge_index])
      estimator <- "ridge fastCrrp endpoint"
    }

    names(initial_beta_std) <- predictor_names

    weights <- 1 / pmax(abs(initial_beta_std), aLASSO_floor)^adaptive_power
    weights <- pmin(weights, aLASSO_weight_cap)
    weights <- weights / min(weights)
    names(weights) <- predictor_names

    list(
      initial_fit = initial_fit,
      initial_beta_std = initial_beta_std,
      initial_estimator = estimator,
      adaptive_weights = weights
    )
  }

  # Evaluate the full risk path with the same Wolbers metric used by
  # cv_ssl_psdh(). The compiled kernel is used when available. Passing an
  # identity design matrix makes its matrix product equal to the already
  # computed out-of-fold risk matrix. Any column containing missing predictions
  # (or any compiled-kernel failure) falls back to wolbers_c() directly.
  score_wolbers_path <- function(y_true, risk_matrix, evaluation_time) {
    y_true <- as.matrix(y_true)
    risk_matrix <- as.matrix(risk_matrix)

    if (ncol(y_true) < 2L) {
      stop("y_true must contain time and status columns.")
    }
    if (nrow(risk_matrix) != nrow(y_true)) {
      stop("Risk matrix and outcome matrix have incompatible row counts.")
    }
    if (length(evaluation_time) != 1L || !is.finite(evaluation_time)) {
      stop("evaluation_time must be one finite value.")
    }

    n_lambda <- ncol(risk_matrix)
    scores <- rep(NA_real_, n_lambda)
    n_available <- colSums(!is.na(risk_matrix))
    complete <- n_available == nrow(risk_matrix)

    if (any(complete) && exists("cpp_cv_cindex", mode = "function")) {
      compiled_scores <- tryCatch(
        cpp_cv_cindex(
          x_test = diag(nrow(risk_matrix)),
          coef = risk_matrix[, complete, drop = FALSE],
          time_test = as.numeric(y_true[, 1]),
          status_test = as.numeric(y_true[, 2]),
          tuning = "wolbers",
          evaluation_time = evaluation_time,
          reverse = FALSE
        ),
        error = function(e) NULL
      )

      if (!is.null(compiled_scores) &&
          length(compiled_scores) == sum(complete)) {
        scores[complete] <- as.numeric(compiled_scores)
      }
    }

    fallback <- which(is.na(scores) & n_available > 0L)
    for (lambda_index in fallback) {
      scores[lambda_index] <- tryCatch(
        wolbers_c(
          y_true = y_true,
          risk_score = risk_matrix[, lambda_index],
          evaluation_time = evaluation_time
        ),
        error = function(e) NA_real_
      )
    }

    scores
  }

  # Predict the cause-1 cumulative incidence at one time point for every
  # lambda in a fitted fastCrrp path. This mirrors predict_from_ssl_psdh(), but
  # handles fastCrrp's one-baseline-hazard-column-per-lambda representation.
  predict_fastcrrp_path <- function(fit, newx, prediction_time, valid_lambda) {
    newx <- as.matrix(newx)
    beta_path <- as.matrix(stats::coef(fit))
    n_lambda <- ncol(beta_path)

    if (ncol(newx) != nrow(beta_path)) {
      stop("newx and the fitted coefficient path have incompatible dimensions.")
    }
    if (length(valid_lambda) != n_lambda) {
      stop("valid_lambda and the fitted coefficient path have incompatible lengths.")
    }

    jumps <- fit$breslowJump
    if (is.null(jumps) || isFALSE(jumps)) {
      stop("fastCrrp did not return Breslow baseline-hazard jumps.")
    }

    jumps <- as.matrix(jumps)
    if (ncol(jumps) == n_lambda + 1L) {
      failure_times <- as.numeric(jumps[, 1])
      jump_path <- jumps[, -1, drop = FALSE]
    } else if (ncol(jumps) == n_lambda &&
               length(fit$uftime) == nrow(jumps)) {
      failure_times <- as.numeric(fit$uftime)
      jump_path <- jumps
    } else {
      stop("Unexpected Breslow-jump dimensions returned by fastCrrp().")
    }

    before_horizon <- which(failure_times <= prediction_time)
    cumulative_baseline <- if (length(before_horizon) == 0L) {
      rep(0, n_lambda)
    } else {
      colSums(jump_path[before_horizon, , drop = FALSE])
    }

    valid_lambda <- valid_lambda &
      is.finite(cumulative_baseline) & cumulative_baseline >= 0

    predicted_risk <- matrix(
      NA_real_,
      nrow = nrow(newx),
      ncol = n_lambda
    )

    for (lambda_index in which(valid_lambda)) {
      if (cumulative_baseline[lambda_index] == 0) {
        predicted_risk[, lambda_index] <- 0
      } else {
        lp <- as.vector(newx %*% beta_path[, lambda_index])
        subject_cumulative_hazard <- exp(
          lp + log(cumulative_baseline[lambda_index])
        )
        predicted_risk[, lambda_index] <-
          -expm1(-subject_cumulative_hazard)
      }
    }

    list(
      risk = predicted_risk,
      valid_lambda = valid_lambda,
      cumulative_baseline = cumulative_baseline
    )
  }

  finite_mean <- function(z) {
    z <- z[is.finite(z)]
    if (length(z) == 0L) NA_real_ else mean(z)
  }

  finite_sd <- function(z) {
    z <- z[is.finite(z)]
    if (length(z) <= 1L) NA_real_ else stats::sd(z)
  }

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

  # Cross-validate one comparator over the complete shared lambda path. For each
  # repetition, every subject receives an out-of-fold predicted cumulative
  # incidence at every lambda before the pooled Wolbers C-index is calculated.
  cv_fastcrrp_wolbers <- function(penalty, method_label,
                                  penalty_factor = NULL, gamma = NULL,
                                  adaptive = FALSE) {
    if (adaptive && !is.null(penalty_factor)) {
      stop("Use either adaptive = TRUE or a fixed penalty_factor, not both.")
    }

    n_lambda <- length(reference_lambda)
    pooled_scores <- matrix(
      NA_real_,
      nrow = comparator_ncv,
      ncol = n_lambda,
      dimnames = list(
        paste0("rep_", seq_len(comparator_ncv)),
        paste0("lambda_", seq_len(n_lambda))
      )
    )
    fold_scores <- array(
      NA_real_,
      dim = c(comparator_ncv, comparator_nfolds, n_lambda),
      dimnames = list(
        paste0("rep_", seq_len(comparator_ncv)),
        paste0("fold_", seq_len(comparator_nfolds)),
        paste0("lambda_", seq_len(n_lambda))
      )
    )
    unavailable <- array(
      FALSE,
      dim = c(comparator_ncv, comparator_nfolds, n_lambda),
      dimnames = dimnames(fold_scores)
    )
    n_predicted <- matrix(
      0L,
      nrow = comparator_ncv,
      ncol = n_lambda,
      dimnames = dimnames(pooled_scores)
    )
    failed_folds <- vector("list", comparator_ncv)

    for (cv_rep in seq_len(comparator_ncv)) {
      risk_all <- matrix(NA_real_, nrow = nobs, ncol = n_lambda)
      failed <- integer(0)

      for (fold in seq_len(comparator_nfolds)) {
        omit_indices <- comparator_fold_indices[[cv_rep]][[fold]]
        train_data <- fit_data[-omit_indices, , drop = FALSE]
        x_test <- as.matrix(
          fit_data[omit_indices, predictor_names, drop = FALSE]
        )

        this_penalty_factor <- penalty_factor
        if (adaptive) {
          adaptive_info <- tryCatch(
            suppressWarnings(
              get_adaptive_lasso_weights(
                data = train_data,
                lambda = reference_lambda,
                warn_on_ridge = FALSE
              )
            ),
            error = function(e) e
          )
          if (inherits(adaptive_info, "error")) {
            message(sprintf(
              "[comparator CV] %s adaptive weights failed in fold %d (rep %d): %s",
              method_label, fold, cv_rep, conditionMessage(adaptive_info)
            ))
            unavailable[cv_rep, fold, ] <- TRUE
            failed <- c(failed, fold)
            next
          }
          this_penalty_factor <- adaptive_info$adaptive_weights
        }

        fit_fold <- tryCatch(
          suppressWarnings(
            fit_fastcrrp(
              penalty = penalty,
              penalty_factor = this_penalty_factor,
              gamma = gamma,
              lambda = reference_lambda,
              data = train_data,
              get_breslow_jumps = TRUE
            )
          ),
          error = function(e) e
        )

        if (inherits(fit_fold, "error")) {
          message(sprintf(
            "[comparator CV] %s fit failed in fold %d (rep %d): %s",
            method_label, fold, cv_rep, conditionMessage(fit_fold)
          ))
          unavailable[cv_rep, fold, ] <- TRUE
          failed <- c(failed, fold)
          next
        }

        if (!isTRUE(all.equal(
          as.numeric(fit_fold$lambda.path),
          as.numeric(reference_lambda),
          tolerance = 1e-12
        ))) {
          message(sprintf(
            "[comparator CV] %s returned a different lambda path in fold %d (rep %d).",
            method_label, fold, cv_rep
          ))
          unavailable[cv_rep, fold, ] <- TRUE
          failed <- c(failed, fold)
          next
        }

        beta_path <- as.matrix(stats::coef(fit_fold))
        valid_lambda <- apply(beta_path, 2L, function(z) all(is.finite(z)))
        if (!is.null(fit_fold$converged) &&
            length(fit_fold$converged) == n_lambda) {
          fold_converged <- as.logical(fit_fold$converged)
          fold_converged[is.na(fold_converged)] <- FALSE
          valid_lambda <- valid_lambda & fold_converged
        }
        prediction <- tryCatch(
          predict_fastcrrp_path(
            fit = fit_fold,
            newx = x_test,
            prediction_time = comparator_eval_time,
            valid_lambda = valid_lambda
          ),
          error = function(e) e
        )
        if (inherits(prediction, "error")) {
          message(sprintf(
            "[comparator CV] %s prediction failed in fold %d (rep %d): %s",
            method_label, fold, cv_rep, conditionMessage(prediction)
          ))
          unavailable[cv_rep, fold, ] <- TRUE
          failed <- c(failed, fold)
          next
        }

        unavailable[cv_rep, fold, ] <- !prediction$valid_lambda
        risk_fold <- prediction$risk
        if (!any(prediction$valid_lambda)) failed <- c(failed, fold)

        risk_all[omit_indices, ] <- risk_fold
        y_test <- comparator_y[omit_indices, , drop = FALSE]
        fold_scores[cv_rep, fold, ] <- score_wolbers_path(
          y_true = y_test,
          risk_matrix = risk_fold,
          evaluation_time = comparator_eval_time
        )
      }

      pooled_scores[cv_rep, ] <- score_wolbers_path(
        y_true = comparator_y,
        risk_matrix = risk_all,
        evaluation_time = comparator_eval_time
      )
      n_predicted[cv_rep, ] <- colSums(!is.na(risk_all))
      failed_folds[[cv_rep]] <- sort(unique(failed))
    }

    score_mean <- apply(pooled_scores, 2L, finite_mean)
    score_sd <- if (comparator_ncv == 1L) {
      rep(NA_real_, n_lambda)
    } else {
      apply(pooled_scores, 2L, finite_sd)
    }

    fold_sd_by_rep <- matrix(
      NA_real_,
      nrow = comparator_ncv,
      ncol = n_lambda
    )
    for (cv_rep in seq_len(comparator_ncv)) {
      for (lambda_index in seq_len(n_lambda)) {
        fold_sd_by_rep[cv_rep, lambda_index] <- finite_sd(
          fold_scores[cv_rep, , lambda_index]
        )
      }
    }
    fold_sd <- apply(fold_sd_by_rep, 2L, finite_mean)

    list(
      method = method_label,
      lambda = reference_lambda,
      score_mean = score_mean,
      score_sd = score_sd,
      fold_sd = fold_sd,
      pooled_scores = pooled_scores,
      fold_scores = fold_scores,
      unavailable = unavailable,
      n_predicted = n_predicted,
      failed_folds = failed_folds,
      n_failed = sum(lengths(failed_folds)),
      foldid = comparator_foldid,
      eval_time = comparator_eval_time,
      eval_quantile = comparator_eval_quantile,
      nfolds = comparator_nfolds,
      ncv = comparator_ncv
    )
  }

  # Select the full-data coefficient vector at the lambda maximizing the mean
  # pooled cross-validated Wolbers C-index. BIC and model-size paths are retained
  # only as diagnostics for backward compatibility; they do not drive selection.
  select_by_cv_wolbers <- function(fit, cv_result, method_label) {
    beta_path <- as.matrix(stats::coef(fit))
    if (nrow(beta_path) != npredictors && ncol(beta_path) == npredictors) {
      beta_path <- t(beta_path)
    }

    lambda_path <- as.numeric(fit$lambda.path)
    if (length(lambda_path) != ncol(beta_path) ||
        !isTRUE(all.equal(
          lambda_path,
          as.numeric(cv_result$lambda),
          tolerance = 1e-12
        ))) {
      stop("CV and full-data lambda paths are incompatible for ", method_label, ".")
    }

    cv_path <- as.numeric(cv_result$score_mean)
    if (length(cv_path) != ncol(beta_path)) {
      stop("CV score and coefficient paths have incompatible lengths for ",
           method_label, ".")
    }

    full_valid <- apply(beta_path, 2L, function(z) all(is.finite(z)))
    if (!is.null(fit$converged) && length(fit$converged) == length(cv_path)) {
      full_converged <- as.logical(fit$converged)
      full_converged[is.na(full_converged)] <- FALSE
      full_valid <- full_valid & full_converged
    }

    valid <- is.finite(cv_path) & full_valid
    if (!any(valid)) {
      stop("No finite cross-validated, converged solution was available for ",
           method_label, ".")
    }

    cv_for_selection <- cv_path
    cv_for_selection[!valid] <- -Inf
    index <- which.max(cv_for_selection)

    beta_std <- as.numeric(beta_path[, index])
    names(beta_std) <- predictor_names

    # x_std = (x - center) / scale, so divide by scale to return coefficients to
    # original predictor units. Centering is absorbed by the baseline hazard.
    beta_raw <- beta_std / x_scale
    names(beta_raw) <- predictor_names

    selected <- beta_std != 0

    # Retain the former BIC/df diagnostics at the CV-selected lambda so existing
    # downstream code that reads these fields does not break.
    df_path <- colSums(beta_path != 0)
    loglik_path <- as.numeric(fit$logLik)
    if (length(loglik_path) != ncol(beta_path)) {
      stop("Log-likelihood and coefficient paths have incompatible lengths.")
    }
    bic_path <- -2 * loglik_path + log(nobs) * df_path

    list(
      method            = method_label,
      index             = index,
      lambda            = lambda_path[index],
      lambda_path       = lambda_path,
      cv_c_index        = cv_path[index],
      cv_c_index_sd     = as.numeric(cv_result$score_sd[index]),
      cv_fold_sd        = as.numeric(cv_result$fold_sd[index]),
      bic               = as.numeric(bic_path[index]),
      df                = as.integer(df_path[index]),
      converged         = if (length(fit$converged) == length(cv_path)) {
        as.logical(fit$converged[index])
      } else {
        NA
      },
      beta_std          = beta_std,
      beta_raw          = beta_raw,
      selected          = selected,
      cv_c_index_path   = cv_path,
      cv_c_index_sd_path = as.numeric(cv_result$score_sd),
      cv_fold_sd_path   = as.numeric(cv_result$fold_sd),
      cv_raw_matrix     = cv_result$pooled_scores,
      cv_fold_scores    = cv_result$fold_scores,
      cv_unavailable    = cv_result$unavailable,
      cv_n_predicted    = cv_result$n_predicted,
      cv_failed_folds   = cv_result$failed_folds,
      cv_n_failed       = cv_result$n_failed,
      cv_foldid         = cv_result$foldid,
      cv_eval_time      = cv_result$eval_time,
      cv_eval_quantile  = cv_result$eval_quantile,
      cv_nfolds         = cv_result$nfolds,
      cv_ncv            = cv_result$ncv,
      selection_metric  = "mean pooled k-fold CV Wolbers C-index",
      bic_path          = bic_path,
      df_path           = df_path
    )
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
  adaptive_full <- get_adaptive_lasso_weights(
    data = fit_data,
    lambda = reference_lambda,
    warn_on_ridge = TRUE
  )
  initial_fit       <- adaptive_full$initial_fit
  initial_beta_std  <- adaptive_full$initial_beta_std
  initial_estimator <- adaptive_full$initial_estimator
  adaptive_weights  <- adaptive_full$adaptive_weights

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

  # Generate the folds once and reuse them across all methods and all lambdas,
  # matching the SSL-PSDH tuning strategy.
  comparator_folds <- generate_foldid(
    nobs = nobs,
    nfolds = comparator_nfolds,
    foldid = NULL,
    ncv = comparator_ncv
  )
  comparator_foldid <- comparator_folds$foldid
  comparator_nfolds <- comparator_folds$nfolds
  comparator_ncv <- comparator_folds$ncv

  comparator_y <- as.matrix(fit_data[, c("TTE", "Status")])
  cause1_times <- comparator_y[comparator_y[, 2] == 1, 1]
  if (length(cause1_times) == 0L) {
    stop("No cause-1 events were available for comparator cross-validation.")
  }
  comparator_eval_time <- as.numeric(
    stats::quantile(cause1_times, comparator_eval_quantile)
  )

  comparator_fold_indices <- lapply(
    seq_len(comparator_ncv),
    function(cv_rep) {
      lapply(
        seq_len(comparator_nfolds),
        function(fold) which(comparator_foldid[, cv_rep] == fold)
      )
    }
  )

  cv_lasso <- cv_fastcrrp_wolbers(
    penalty = "LASSO",
    method_label = "LASSO"
  )
  cv_alasso <- cv_fastcrrp_wolbers(
    penalty = "LASSO",
    method_label = "adaptive LASSO",
    adaptive = TRUE
  )
  cv_scad <- cv_fastcrrp_wolbers(
    penalty = "SCAD",
    method_label = "SCAD",
    gamma = 3.7
  )
  cv_mcp <- cv_fastcrrp_wolbers(
    penalty = "MCP",
    method_label = "MCP",
    gamma = 2.7
  )

  selected_models <- list(
    LASSO  = select_by_cv_wolbers(fit_lasso,  cv_lasso,  "LASSO"),
    aLASSO = select_by_cv_wolbers(fit_alasso, cv_alasso, "adaptive LASSO"),
    SCAD   = select_by_cv_wolbers(fit_scad,   cv_scad,   "SCAD"),
    MCP    = select_by_cv_wolbers(fit_mcp,    cv_mcp,    "MCP")
  )

  other_time2 <- Sys.time()
  other_time  <- other_time2 - other_time1

