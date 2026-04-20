#' Generate cross-validation fold assignments
#'
#' Creates an \eqn{n \times \code{ncv}} integer matrix assigning each
#' observation to one of \code{nfolds} folds, optionally repeated
#' \code{ncv} times with independent random shuffles.  If \code{foldid} is
#' already supplied it is returned as-is (after coercing to matrix), allowing
#' callers to reuse the same splits across multiple hyperparameter evaluations.
#'
#' @param nobs Integer. Number of observations \eqn{n}.
#' @param nfolds Integer. Number of folds.  Automatically capped at
#'   \code{nobs}; when equal to \code{nobs} the procedure is leave-one-out
#'   and \code{ncv} is forced to \code{1}. Default \code{10}.
#' @param foldid Integer matrix (\eqn{n \times \code{ncv}}) of pre-specified
#'   fold assignments, or \code{NULL} (default) to generate them randomly.
#'   When non-\code{NULL} the other arguments are ignored.
#' @param ncv Integer. Number of independent CV repetitions (columns of the
#'   returned \code{foldid} matrix).  Forced to \code{1} for leave-one-out.
#'   Default \code{1}.
#'
#' @returns A list with three elements:
#'   \describe{
#'     \item{\code{foldid}}{Integer matrix of dimensions
#'       \eqn{n \times \code{ncv}}.  Entry \code{[i, k]} gives the fold
#'       index assigned to observation \eqn{i} in repetition \eqn{k}.}
#'     \item{\code{nfolds}}{Integer. Effective number of folds used
#'       (may differ from the argument if capped at \code{nobs}).}
#'     \item{\code{ncv}}{Integer. Effective number of repetitions used.}
#'   }
#'
#' @seealso \code{\link{cv_ssl_psdh}}, \code{\link{tune_ssl_psdh}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fols <- generate_foldid(nobs = 100, nfolds = 10, ncv = 3)
#' dim(fols$foldid)  # 100 x 3
#' }
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
