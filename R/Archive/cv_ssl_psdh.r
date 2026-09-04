#' Cross-validate a spike-and-slab PSDH model
#'
#' Evaluates a specific spike-and-slab scale pair \code{(s0, s1)} via
#' \eqn{k}-fold cross-validation, using the IPCW concordance index
#' (from \code{\link{wolbers_c}}) as the performance metric.
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
#' @param init Optional warm-start structure. A list of length \code{ncv},
#'   each element a list of length \code{nfolds} holding a length-\eqn{p}
#'   numeric vector of starting coefficients for that (repetition, fold)
#'   training fit; \code{NULL} entries (or \code{init = NULL}, the default)
#'   fall back to the per-fit LASSO-CV cold start. Element \code{init[[k]][[i]]}
#'   is passed as \code{init} to \code{\link{fit_ssl_psdh}} for fold \code{i}
#'   of repetition \code{k}. The matching \code{$fold_coefs} field returned by
#'   this function can be fed straight back in as \code{init} for a subsequent
#'   scale pair (this is how \code{\link{tune_ssl_psdh}} warm-starts).
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
#'     \item{\code{fold_coefs}}{A list of length \code{ncv}, each element a
#'       list of length \code{nfolds} giving the fitted coefficient vector
#'       for that (repetition, fold) training fit (\code{NULL} where the fit
#'       failed). Suitable to pass back as \code{init} for a warm start.}
#'   }
#'
#' @importFrom stats quantile sd
#'
#' @seealso \code{\link{fit_ssl_psdh}}, \code{\link{tune_ssl_psdh}},
#'   \code{\link{wolbers_c}}, \code{\link{generate_foldid}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fit  <- fit_ssl_psdh(x, y)
#' fols <- generate_foldid(nobs = nrow(x), nfolds = 5)
#' cv_ssl_psdh(fit, foldid = fols$foldid, s0 = 0.04, s1 = 0.5)
#' }
cv_ssl_psdh <- function(object, foldid, s0, s1, ncv=1, eval_quantile = 0.5,
                        init = NULL, ...) {
  # Extract data
  y <- object$y
  x <- object$x
  n <- NROW(y)
  nfolds <- max(foldid)

  eval_time <- as.numeric(quantile(y[y[,2] == 1, 1], eval_quantile))

  pooled_measures <- rep(NA_real_, ncv)
  fold_measures_list <- vector("list", ncv)
  failed_folds_list <- vector("list", ncv)

  # Per-(rep, fold) fitted coefficients, returned so a caller can warm-start
  # a subsequent scale pair. Stays NULL for folds whose fit fails.
  fitted_coefs <- lapply(seq_len(ncv), function(k) vector("list", nfolds))

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

      # Warm-start coefficients for this (rep, fold), if supplied; else NULL
      # so fit_ssl_psdh derives its own LASSO-CV cold start.
      this_init <- if (is.null(init)) NULL else init[[k]][[i]]

      # RE-FIT (errors are captured rather than silenced so we can report them)
      suppressWarnings({
        fit <- try(fit_ssl_psdh(x = x_train, y = y_train, ss = c(s0, s1),
                                maxit = 50,
                                epsilon = 1e-04, init = this_init, ...),
                   silent = TRUE)
      })

      if (inherits(fit, "try-error")) {
        err_msg <- attr(fit, "condition")$message
        message(sprintf(
          "[cv_ssl_psdh] fit failed in fold %d (rep %d) at s0=%g, s1=%g: %s",
          i, k, s0, s1, err_msg))
        failed <- c(failed, i)
        next  # try the remaining folds rather than aborting the repetition
      }

      # Record fitted coefficients for a possible warm start of the next pair.
      # Source from the raw model object (a plain length-p vector) and only
      # store when well-formed: assigning NULL to a list slot would delete it
      # and shrink the per-fold list, corrupting later warm-start reads.
      fitted_coef <- as.numeric(fit$final_model_object$coef)
      if (length(fitted_coef) == ncol(x_train)) fitted_coefs[[k]][[i]] <- fitted_coef

      # PREDICT on hold-out set
      x_test  <- x[omit_indices, , drop = FALSE]
      lp_fold <- predict_from_ssl_psdh(fit, newx = x_test, eval_time)
      lp_all[omit_indices] <- lp_fold

      # Per-fold C-index, evaluated on this fold's held-out subjects only.
      # NA is possible if the held-out subset has no cause-1 events at or
      # before eval_time, or if the IPCW denominator is zero.
      y_test <- as.matrix(y[omit_indices, ])
      fold_scores[i] <- tryCatch(
        wolbers_c(y_test, lp_fold, eval_time),
        error = function(e) NA_real_
      )
    }

    fold_measures_list[[k]] <- fold_scores
    failed_folds_list[[k]]  <- failed

    # Pooled cross-validated C-index across folds for this repetition
    pooled_measures[k] <- wolbers_c(y, lp_all, eval_time)
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
       n_failed      = sum(lengths(failed_folds_list)),
       fold_coefs    = fitted_coefs)
}
