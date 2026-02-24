#' cv_ssl_psdh
#'
#' @param object
#' @param foldid
#' @param s0
#' @param s1
#' @param ncv
#' @param eval_quantile
#'
#' @returns
#' @export
#'
#' @examples
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

      if(min(y_train[,2])>0){
        print("No censoring in fold")
        break}
      if(all(y_train[,2] != 1)){
        print("No events of interest in fold")
        break}

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
        print("fit_ssl_psdh failed in fold")
        print(i)

        print("Checking Outcome")
        print(table(y_train[,2]))

        failed_df[[i]]<<- cbind(y_train, x_train)
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
