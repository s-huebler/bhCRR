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
#' @param initial_sparsity Numeric in \code{(0, 1)}. Starting value for the
#'   global mixture probability (prior proportion of active features).
#'   Default \code{0.05}.
#' @param reuse_lambda Logical. When \code{TRUE} (default), the expensive
#'   \code{cv_fastCrrp} global-shrinkage search is run only once per \code{s0}
#'   value: it is tuned at that \code{s0}'s smallest valid \code{s1}, and the
#'   resulting per-fold lambdas are reused (via \code{fixed_global_shrinkage})
#'   for the remaining \code{s1} values at the same \code{s0}. Set to
#'   \code{FALSE} to tune the shrinkage independently for every \code{(s0, s1)}
#'   pair (the previous behaviour).
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
tune_ssl_psdh <- function(object, s0_seq, s1_seq, nfolds=10, ncv=1, foldid=NULL, initial_sparsity, reuse_lambda = TRUE) {

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

  # 3. Iterate over the valid combinations, grouped by s0.
  #
  # For each s0 (in s0_seq order) the valid s1 values are visited in
  # ascending order.  When reuse_lambda is TRUE, the first (smallest) s1
  # tunes the global shrinkage via cv_fastCrrp and the resulting per-fold
  # lambdas are cached, then reused for the remaining s1 values at that s0
  # so the expensive search runs only once per s0.

  # Helper: turn a cv_ssl_psdh result (or try-error) into a result row.
  build_row <- function(current_s0, current_s1, cv_res) {
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
  }

  # Unique s0 values, in the order they appear in s0_seq
  s0_levels <- unique(valid_grid$s0)

  results_list <- list()

  for (current_s0 in s0_levels) {
    # s1 values valid for this s0, ascending
    s1_for_s0 <- sort(valid_grid$s1[valid_grid$s0 == current_s0])

    cached_lambdas <- NULL  # per-fold (ncv x nfolds) lambdas from first s1

    for (j in seq_along(s1_for_s0)) {
      current_s1 <- s1_for_s0[j]

      # First s1 of this s0 always tunes; later s1 reuse cached lambdas
      fixed_fgs <- if (reuse_lambda && j > 1) cached_lambdas else NULL

      cv_res <- try(cv_ssl_psdh(object,
                                foldid        = foldid,
                                s0            = current_s0,
                                s1            = current_s1,
                                ncv           = ncv,
                                eval_quantile = 0.5,
                                initial_sparsity = initial_sparsity,
                                fixed_global_shrinkage = fixed_fgs))

      # Cache the per-fold lambdas tuned at the first s1 for reuse
      if (reuse_lambda && j == 1 && !inherits(cv_res, "try-error")) {
        cached_lambdas <- cv_res$fold_lambdas
      }

      results_list[[length(results_list) + 1L]] <-
        build_row(current_s0, current_s1, cv_res)
    }
  }

  # 4. Combine and return results
  as.data.frame(do.call(rbind, results_list))

  })
}
