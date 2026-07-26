#' Cross-validated Tuning for fastCrrp (RcppArmadillo backend)
#'
#' A drop-in alternative to [cv_fastCrrp()] that delegates the cross-validation
#' metric computation to compiled RcppArmadillo kernels. Model fitting is still
#' performed by [fastcmprsk::fastCrrp()] (an external, already-compiled
#' routine); the speedup comes from moving the per-lambda prediction and the
#' O(n^2) concordance computations into C++ (see \code{src/cv_fastcrrp.cpp}).
#'
#' The fold assignment, fitting calls, and returned object are identical in
#' structure to [cv_fastCrrp()], so results should match up to the C-index
#' implementation details documented below.
#'
#' @param x Numeric predictor matrix (n x p).
#' @param time Numeric vector of event/censoring times.
#' @param status Numeric vector of status codes (0 = censored, 1 = cause-1
#'   event, 2 = competing event).
#' @param k Integer number of cross-validation folds.
#' @param penalty Penalty passed to [fastcmprsk::fastCrrp()].
#' @param lambda_path Optional numeric vector of lambda values; if \code{NULL}
#'   the path is chosen by [fastcmprsk::fastCrrp()].
#' @param tuning Either \code{"normal"} (Harrell C-index, matching
#'   \code{survival::concordance}) or \code{"wolbers"} (IPCW competing-risks
#'   C-index, matching [wolbers_c()]).
#' @param eval_quantile Quantile of the training times used as the evaluation
#'   horizon for \code{tuning = "wolbers"}.
#' @param reverse Logical orientation flag for the \code{"normal"} Harrell
#'   C-index. \code{FALSE} (default) matches the \code{survival::concordance}
#'   default. Ignored when \code{tuning = "wolbers"}.
#'
#' @returns A list with elements \code{lambda}, \code{cv_c_index},
#'   \code{lambda_min}, \code{full_model}, and \code{cv_raw_matrix}, matching
#'   [cv_fastCrrp()].
#'
#' @seealso [cv_fastCrrp()], [wolbers_c()]
#'
#' @importFrom fastcmprsk fastCrrp Crisk
#' @importFrom stats quantile
#' @export
#'
#' @examples
#' \dontrun{
#' cv_fastCrrp_cpp(x, time, status, k = 5, penalty = "LASSO",
#'                 tuning = "wolbers", eval_quantile = 0.5)
#' }
cv_fastCrrp_cpp <- function(x, time, status, k = 5, penalty = "LASSO",
                            lambda_path = NULL, tuning = "normal",
                            eval_quantile = 0.5, reverse = FALSE) {

  if (!tuning %in% c("normal", "wolbers")) {
    stop("`tuning` must be 'normal' or 'wolbers'.")
  }

  x      <- as.matrix(x)
  time   <- as.numeric(time)
  status <- as.numeric(status)
  n      <- nrow(x)

  # Assign CV folds reproducibly without leaking the seed into the caller's
  # RNG stream: save the current global seed, set our fixed seed for the fold
  # assignment only, then restore the caller's state.
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv))
    get(".Random.seed", envir = .GlobalEnv) else NULL
  set.seed(9134)
  fold_ids <- sample(rep(1:k, length.out = n))
  if (!is.null(old_seed)) {
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  } else {
    rm(".Random.seed", envir = .GlobalEnv)
  }

  # 1. Fit the full model
  fit_full <- fastCrrp(Crisk(time, status, cencode = 0, failcode = 1) ~ x,
                       penalty = penalty,
                       lambda  = lambda_path)

  lambdas    <- fit_full$lambda
  n_lambda   <- length(lambdas)
  cv_metrics <- matrix(NA_real_, nrow = k, ncol = n_lambda)

  for (i in 1:k) {
    test_idx <- which(fold_ids == i)

    x_train      <- x[-test_idx, , drop = FALSE]
    time_train   <- time[-test_idx]
    status_train <- status[-test_idx]

    x_test       <- x[test_idx, , drop = FALSE]
    time_test    <- time[test_idx]
    status_test  <- status[test_idx]

    # 2. Fit training folds (external compiled routine)
    fit_train <- fastCrrp(Crisk(time_train, status_train, cencode = 0, failcode = 1) ~ x_train,
                          penalty = penalty,
                          lambda  = lambdas)

    coef_mat <- as.matrix(fit_train$coef)  # (p x n_lambda)

    # 3. Predictions + concordance for every lambda, computed in C++.
    eval_time <- if (tuning == "wolbers") {
      as.numeric(quantile(time_train[status_train == 1], eval_quantile))
    } else {
      0
    }

    cv_metrics[i, ] <- cpp_cv_cindex(
      x_test          = x_test,
      coef            = coef_mat,
      time_test       = time_test,
      status_test     = status_test,
      tuning          = tuning,
      evaluation_time = eval_time,
      reverse         = reverse
    )
  }

  mean_cv_metrics <- colMeans(cv_metrics, na.rm = TRUE)
  opt_idx    <- which.max(mean_cv_metrics)
  lambda_min <- lambdas[opt_idx[1]]

  return(list(
    lambda        = lambdas,
    cv_c_index    = mean_cv_metrics,
    lambda_min    = lambda_min,
    full_model    = fit_full,
    cv_raw_matrix = cv_metrics
  ))
}
