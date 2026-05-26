#' Cross-validate a spike-and-slab PSDH model
#'
#' Evaluates a specific spike-and-slab scale pair \code{(s0, s1)} via
#' \eqn{k}-fold cross-validation, using the IPCW concordance index
#' (from \code{\link{measure_ssl_psdh}}) as the performance metric.
#' Optionally repeats cross-validation \code{ncv} times with different
#' random splits and averages the results.
#'
#' @param object A fitted model object returned by \code{\link{fit_ssl_psdh}},
#'   used to extract \code{$x} (feature matrix) and \code{$y} (outcome
#'   matrix) for re-fitting on each training fold.
#' @param foldid Integer matrix of dimensions \eqn{n \times \code{ncv}}.
#'   Each column assigns observations to one of \code{nfolds} folds for one
#'   CV repetition; typically produced by \code{\link{generate_foldid}}.
#' @param s0 Numeric scalar. Spike scale parameter (small value enforcing
#'   shrinkage toward zero for inactive features).
#' @param s1 Numeric scalar. Slab scale parameter (larger value permitting
#'   non-zero effects for active features). Must satisfy \code{s1 > s0}.
#' @param ncv Integer. Number of independent cross-validation repetitions
#'   over which results are averaged. Default \code{1}.
#' @param eval_quantile Numeric in \code{(0, 1)}. Quantile of observed
#'   event times used as the evaluation horizon for the C-index.
#'   Default \code{0.5} (median event time).
#' @param initial_sparsity Numeric in \code{(0, 1)}. Starting value for the
#'   global mixture probability (prior proportion of active features).
#'   Default \code{0.05}.
#'
#' @returns A list with the following elements:
#'   \describe{
#'     \item{\code{measures}}{A \eqn{2 \times 1} matrix with rows
#'       \code{"mean"} and \code{"sd"} giving the mean and standard
#'       deviation of the pooled cross-validated C-index across the
#'       \code{ncv} repetitions.  SD is \code{NA} when \code{ncv = 1}.}
#'     \item{\code{fold_scores}}{An \eqn{\code{ncv} \times \code{nfolds}}
#'       matrix of per-fold C-indices (one row per CV repetition).
#'       Entries may be \code{NA} for folds in which the fit failed or
#'       the held-out subset had no cause-1 events.}
#'     \item{\code{fold_sd}}{Numeric scalar.  Mean over \code{ncv}
#'       repetitions of the within-repetition SD of the per-fold
#'       C-indices.  Provides a measure of fold-to-fold variability
#'       even when \code{ncv = 1}.}
#'     \item{\code{failed_folds}}{Length-\code{ncv} list of integer
#'       vectors, each giving the fold indices in which
#'       \code{fit_ssl_psdh} failed (empty if all folds succeeded).}
#'     \item{\code{n_failed}}{Total number of (rep, fold) pairs in
#'       which \code{fit_ssl_psdh} failed.}
#'   }
#'
#' @importFrom stats quantile sd
#'
#' @seealso \code{\link{fit_ssl_psdh}}, \code{\link{tune_ssl_psdh}},
#'   \code{\link{measure_ssl_psdh}}, \code{\link{generate_foldid}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fit  <- fit_ssl_psdh(x, y)
#' fols <- generate_foldid(nobs = nrow(x), nfolds = 5)
#' cv_ssl_psdh(fit, foldid = fols$foldid, s0 = 0.04, s1 = 0.5)
#' }
cv_ssl_psdh <- function(object, foldid, s0, s1, ncv=1, eval_quantile = 0.5, initial_sparsity) {
  # Extract data
  y <- object$y
  x <- object$x
  n <- NROW(y)
  nfolds <- max(foldid)

  eval_time <- as.numeric(quantile(y[,1], eval_quantile))

  pooled_measures <- rep(NA_real_, ncv)
  fold_measures_list <- vector("list", ncv)
  failed_folds_list <- vector("list", ncv)

  for (k in 1:ncv) {
    lp_all       <- rep(NA_real_, n)
    fold_scores  <- rep(NA_real_, nfolds)
    failed       <- integer(0)

    for (i in 1:nfolds) {
      # Identify hold-out set
      omit_indices <- which(foldid[, k] == i)

      # Create training sets
      y_train <- as.matrix(y[-omit_indices, ])
      x_train <- as.matrix(x[-omit_indices, , drop = FALSE])

      # RE-FIT (errors are captured rather than silenced so we can report them)
      suppressWarnings({
        fit <- try(fit_ssl_psdh(x = x_train, y = y_train, ss = c(s0, s1),
                                initial_sparsity = initial_sparsity,
                                maxit = 50,
                                epsilon = 1e-04), silent = TRUE)
      })

      if (inherits(fit, "try-error")) {
        err_msg <- attr(fit, "condition")$message
        message(sprintf(
          "[cv_ssl_psdh] fit failed in fold %d (rep %d) at s0=%g, s1=%g: %s",
          i, k, s0, s1, err_msg))
        failed <- c(failed, i)
        next  # try the remaining folds rather than aborting the repetition
      }

      # PREDICT on hold-out set
      x_test  <- x[omit_indices, , drop = FALSE]
      lp_fold <- predict_from_ssl_psdh(fit, newx = x_test, eval_time)
      lp_all[omit_indices] <- lp_fold

      # Per-fold C-index, evaluated on this fold's held-out subjects only.
      # NA is possible if the held-out subset has no cause-1 events at or
      # before eval_time, or if the IPCW denominator is zero.
      y_test <- as.matrix(y[omit_indices, ])
      fold_scores[i] <- tryCatch(
        measure_ssl_psdh(y_test, lp_fold, eval_time),
        error = function(e) NA_real_
      )
    }

    fold_measures_list[[k]] <- fold_scores
    failed_folds_list[[k]]  <- failed

    # Pooled cross-validated C-index across folds for this repetition
    pooled_measures[k] <- measure_ssl_psdh(y, lp_all, eval_time)
  }

  # Aggregate across NCV repetitions
  final_measures <- matrix(
    c(mean(pooled_measures, na.rm = TRUE),
      if (ncv == 1) NA_real_ else sd(pooled_measures, na.rm = TRUE)),
    nrow = 2,
    dimnames = list(c("mean", "sd"), "score")
  )

  # Within-repetition SD of per-fold C-indices, averaged across repetitions
  fold_mat <- do.call(rbind, fold_measures_list)  # ncv x nfolds
  fold_sd_per_rep <- apply(fold_mat, 1, sd, na.rm = TRUE)
  fold_sd_mean    <- mean(fold_sd_per_rep, na.rm = TRUE)

  list(measures      = final_measures,
       fold_scores   = fold_mat,
       fold_sd       = fold_sd_mean,
       failed_folds  = failed_folds_list,
       n_failed      = sum(lengths(failed_folds_list)))
}
