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
#'
#' @returns A list with one element:
#'   \describe{
#'     \item{\code{measures}}{A \eqn{2 \times 1} matrix with rows
#'       \code{"mean"} and \code{"sd"} giving the mean and standard
#'       deviation of the C-index across the \code{ncv} repetitions.}
#'   }
#'
#' @importFrom stats quantile
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
cv_ssl_psdh <- function(object, foldid, s0, s1, ncv=1, eval_quantile = 0.5) {
  # Extract data
  y <- object$y
  x <- object$x
  n <- NROW(y)
  nfolds <- max(foldid)

  eval_time <- as.numeric(quantile(y[,1], eval_quantile))
  #print(eval_time)

  measures_list <- list()

  for (k in 1:ncv) {
    lp_all <- rep(NA, n)

    for (i in 1:nfolds) {
      # Identify hold-out set
      omit_indices <- which(foldid[, k] == i)

      # create training sets
      y_train <- y[-omit_indices,]


      #print("Checking Outcome")
      #print(table(y_train[,2]))

      # if(min(y_train[,2])>0){
      #   print("No censoring in fold")
      #   break}
      # if(all(y_train[,2] != 1)){
      #   print("No events of interest in fold")
      #   break}

      y_train <- as.matrix(y_train)



      #print(y_train)
      #print(min(y_train[,2]))
      #print(i)

      x_train <- x[-omit_indices, , drop=FALSE]
      x_train <- as.matrix(x_train)




      # RE-FIT
      suppressWarnings({
        fit <- try(fit_ssl_psdh(x = x_train, y = y_train, ss= c(s0, s1),
                                initial_sparsity = 0.05,
                                maxit = 50,
                                epsilon=1e-04), silent = TRUE)
      })

      if("try-error" %in% class(fit)){
        #print("fit_ssl_psdh failed in fold")
        #print(i)

        #print("Checking Outcome")
        #print(table(y_train[,2]))

        #failed_df[[i]]<<- cbind(y_train, x_train)
        break
      }


      # PREDICT on hold-out set
      x_test <- x[omit_indices, , drop=FALSE]
      lp_all[omit_indices] <- predict_from_ssl_psdh(fit, newx = x_test, eval_time)
      }




    measures_list[[k]] <- measure_ssl_psdh(y, lp_all, eval_time)
  }

  # Aggregate results (Mean and SD across NCV repetitions)
  measures_mat <- do.call(rbind, measures_list)

  if (nrow(measures_mat) == 1) {
    final_measures <- rbind(mean = measures_mat, sd = NA)
  } else {
    final_measures <- rbind(mean = colMeans(measures_mat),
                            sd = apply(measures_mat, 2, sd))
  }

  list(measures = final_measures)
}
