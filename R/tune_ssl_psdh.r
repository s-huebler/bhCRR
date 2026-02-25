#' tune_ssl_psdh
#'
#' @param object
#' @param s0_seq
#' @param s1_seq
#' @param nfolds
#' @param ncv
#' @param foldid
#'
#' @returns
#' @export
#'
#' @examples
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
