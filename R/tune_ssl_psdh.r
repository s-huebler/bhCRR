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
#'   columns \code{s0}, \code{s1}, and one or more performance-metric
#'   columns (e.g. \code{score_mean}, \code{score_sd}) from
#'   \code{\link{cv_ssl_psdh}}.
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
  # We loop through the rows of the valid_grid
  results_list <- lapply(1:nrow(valid_grid), function(i) {

    current_s0 <- valid_grid$s0[i]
    current_s1 <- valid_grid$s1[i]

    # Call the CV function with the specific pair
    # Note: Passed 'current_s1' to the cv function


    cv_res <- try(cv_ssl_psdh(object,
                              foldid = foldid,
                              s0 = current_s0,
                              s1 = current_s1,
                              ncv = ncv,
                              eval_quantile = 0.5))

    if("try-error" %in% class(cv_res)){
     # print(paste0("Failed when s0 = ", current_s0,
    #               "; s1 = ", current_s1))


      cv_res$measures <- c(NA,NA)
    }

    # Extract mean and sd from the results
    stats <- as.vector(cv_res$measures)

    # Handle cases where cv_res$measures might have different column names
    # Assuming standard structure: columns are metrics, rows are "mean"/"sd"
    if(is.null(colnames(cv_res$measures))) {
      metric_names <- "score" # Fallback if no names
    } else {
      metric_names <- colnames(cv_res$measures)
    }

    names(stats) <- c(paste0(metric_names, "_mean"),
                      paste0(metric_names, "_sd"))

    # Return vector including the parameters used
    return(c(s0 = current_s0, s1 = current_s1, stats))
  })

  # 4. Combine and return results
  results <- as.data.frame(do.call(rbind, results_list))

  return(results)
}
