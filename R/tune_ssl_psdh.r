#' Grid-search tuning for spike-and-slab PSDH scale parameters
#'
#' Searches over all valid \code{(s0, s1)} combinations in the Cartesian
#' product of \code{s0_seq} and \code{s1_seq} (filtering out pairs where
#' \code{s1 <= s0}), evaluating each via \code{\link{cv_ssl_psdh}}.  The
#' same fold assignments are reused across all hyperparameter pairs so that
#' performance differences reflect only the scale parameters.
#'
#' @param object A fitted model object returned by \code{\link{fit_ssl_psdh}},
#'   supplying \code{$x} and \code{$y} for cross-validation re-fitting.
#' @param s0_seq Numeric vector. Candidate spike scale values to search over.
#'   All values should be small and positive (e.g. \code{seq(0.005, 0.1,
#'   length.out = 20)}).
#' @param s1_seq Numeric vector. Candidate slab scale values to search over.
#'   Values must be larger than the corresponding \code{s0} (pairs violating
#'   \code{s1 > s0} are silently dropped).
#' @param nfolds Integer. Number of cross-validation folds. Default \code{10}.
#' @param ncv Integer. Number of independent CV repetitions per
#'   hyperparameter pair. Default \code{1}.
#' @param foldid Integer matrix (\eqn{n \times \code{ncv}}) of pre-specified
#'   fold assignments.  If \code{NULL} (default) folds are generated
#'   internally by \code{\link{generate_foldid}} and shared across all
#'   hyperparameter pairs.
#'
#' @returns A data frame with one row per valid \code{(s0, s1)} pair and
#'   columns:
#'   \describe{
#'     \item{\code{s0}, \code{s1}}{The hyperparameter pair.}
#'     \item{\code{score_mean}}{Pooled cross-validated C-index averaged
#'       over \code{ncv} repetitions.}
#'     \item{\code{score_sd}}{SD of the pooled C-index across \code{ncv}
#'       repetitions; \code{NA} when \code{ncv = 1}.}
#'     \item{\code{score_fold_sd}}{Mean (over \code{ncv} repetitions) of
#'       the within-repetition SD of the per-fold C-indices.  Captures
#'       fold-to-fold variability even when \code{ncv = 1}.}
#'     \item{\code{n_failed_folds}}{Total number of (rep, fold) pairs in
#'       which \code{fit_ssl_psdh} failed for this hyperparameter pair.}
#'   }
#'
#' @seealso \code{\link{fit_ssl_psdh}}, \code{\link{cv_ssl_psdh}},
#'   \code{\link{generate_foldid}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fit   <- fit_ssl_psdh(x, y)
#' tunes <- tune_ssl_psdh(fit,
#'                        s0_seq = seq(0.005, 0.1, length.out = 10),
#'                        s1_seq = seq(0.3,   0.9, length.out = 10))
#' tunes[which.max(tunes$score_mean), ]
#' }
tune_ssl_psdh <- function(object, s0_seq, s1_seq, nfolds=10, ncv=1, foldid=NULL) {

  .dedupe_warnings({

  # 1. Generate folds once so every hyperparameter is tested on identical splits
  n <- NROW(object$y)
  if (is.null(foldid)) {
    fol <- generate_foldid(nobs=n, nfolds=nfolds, foldid=foldid, ncv=ncv)
    foldid <- fol$foldid
  }

  # 2. Generate and Filter Hyperparameter Grid
  # Create all possible combinations of s0 and s1
  param_grid <- expand.grid(s0 = s0_seq, s1 = s1_seq)

  # Apply constraint: s1 must be greater than s0
  valid_grid <- param_grid[param_grid$s1 > param_grid$s0, ]

  # Safety check
  if (nrow(valid_grid) == 0) {
    stop("No valid hyperparameter combinations found. Ensure at least one value in 's1_seq' is greater than a value in 's0_seq'.")
  }

  # 3. Iterate over the valid combinations
  results_list <- lapply(1:nrow(valid_grid), function(i) {

    current_s0 <- valid_grid$s0[i]
    current_s1 <- valid_grid$s1[i]

    cv_res <- try(cv_ssl_psdh(object,
                              foldid        = foldid,
                              s0            = current_s0,
                              s1            = current_s1,
                              ncv           = ncv,
                              eval_quantile = 0.5))

    # If cv_ssl_psdh itself errored, surface it and return an explicit NA row
    if (inherits(cv_res, "try-error")) {
      err_msg <- attr(cv_res, "condition")$message
      message(sprintf(
        "[tune_ssl_psdh] cv_ssl_psdh failed at s0=%g, s1=%g: %s",
        current_s0, current_s1, err_msg))
      return(c(s0             = current_s0,
               s1             = current_s1,
               score_mean     = NA_real_,
               score_sd       = NA_real_,
               score_fold_sd  = NA_real_,
               n_failed_folds = NA_real_))
    }

    metric_names <- if (is.null(colnames(cv_res$measures))) {
      "score"
    } else {
      colnames(cv_res$measures)
    }

    means <- as.numeric(cv_res$measures["mean", ])
    sds   <- as.numeric(cv_res$measures["sd",   ])
    names(means) <- paste0(metric_names, "_mean")
    names(sds)   <- paste0(metric_names, "_sd")

    fold_sd <- setNames(cv_res$fold_sd, paste0(metric_names, "_fold_sd"))

    c(s0 = current_s0,
      s1 = current_s1,
      means,
      sds,
      fold_sd,
      n_failed_folds = cv_res$n_failed)
  })

  # 4. Combine and return results
  as.data.frame(do.call(rbind, results_list))

  })
}
