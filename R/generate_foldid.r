#' generate_foldid
#'
#' @param nobs
#' @param nfolds
#' @param foldid
#' @param ncv
#'
#' @returns
#' @export
#'
#' @examples
generate_foldid <- function(nobs, nfolds=10, foldid=NULL, ncv=1) {
  if (nfolds > nobs) nfolds <- nobs
  if (nfolds == nobs) ncv <- 1

  if (is.null(foldid)) {
    foldid <- array(NA, c(nobs, ncv))
    for(j in 1:ncv) {
      # Randomly shuffle fold assignments
      foldid[, j] <- sample(rep(seq(nfolds), length=nobs))
    }
  }
  return(list(foldid = as.matrix(foldid), nfolds = nfolds, ncv = ncv))
}
