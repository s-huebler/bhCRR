#' Cross-validated Tuning for fastCrrp
#'
#' @param x
#' @param time
#' @param status
#' @param k
#' @param penalty
#' @param lambda_path
#' @param tuning
#' @param eval_quantile
#'
#' @returns
#' @export
#'
#' @examples
cv_fastCrrp <- function(x, time, status, k = 5, penalty = "LASSO",
                        lambda_path = NULL, tuning = "normal",
                        eval_quantile = 0.5) {

  n <- nrow(x)
  set.seed(42)
  fold_ids <- sample(rep(1:k, length.out = n))

  # 1. Fit the full model
  fit_full <- fastCrrp(Crisk(time, status, cencode = 0, failcode = 1) ~ x,
                       penalty = penalty,
                       lambda  = lambda_path)

  lambdas    <- fit_full$lambda
  n_lambda   <- length(lambdas)
  cv_metrics <- matrix(NA, nrow = k, ncol = n_lambda)

  for (i in 1:k) {
    test_idx <- which(fold_ids == i)

    x_train      <- x[-test_idx, , drop = FALSE]
    time_train   <- time[-test_idx]
    status_train <- status[-test_idx]

    x_test       <- x[test_idx, , drop = FALSE]
    time_test    <- time[test_idx]
    status_test  <- status[test_idx]

    # 2. Fit training folds
    fit_train <- fastCrrp(Crisk(time_train, status_train, cencode = 0, failcode = 1) ~ x_train,
                          penalty = penalty,
                          lambda  = lambdas)

    # 3. MANUALLY CALCULATE PREDICTIONS (X %*% beta)
    # x_test is (n_test x p), fit_train$coef is (p x n_lambda)
    # The result is an (n_test x n_lambda) matrix of risk scores
    pred_risk <- as.matrix(x_test) %*% as.matrix(fit_train$coef)

    if (tuning == "normal") {
      for (j in 1:n_lambda) {
        risk_scores <- pred_risk[, j]
        cv_metrics[i, j] <- tryCatch(
          survival::concordance(
            Surv(time_test, status_test == 1) ~ risk_scores
          )$concordance,
          error = function(e) NA_real_
        )
      }
    } else if (tuning == "wolbers") {
      eval_time  <- as.numeric(quantile(time_train[status_train == 1], eval_quantile))
      y_test_mat <- cbind(time_test, status_test)
      for (j in 1:n_lambda) {
        risk_scores <- pred_risk[, j]
        cv_metrics[i, j] <- tryCatch(
          measure_ssl_psdh(y_true          = y_test_mat,
                           risk_score      = risk_scores,
                           evaluation_time = eval_time),
          error = function(e) NA_real_
        )
      }
    } else {
      stop("`tuning` must be 'normal' or 'wolbers'.")
    }
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
