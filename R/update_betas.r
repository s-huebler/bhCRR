#' M-step coefficient update via penalised subdistribution hazard regression
#'
#' Fits a LASSO-penalised Fine-Gray subdistribution hazard model using
#' \code{fastcmprsk::fastCrrp}.  This is the M-step of the EM algorithm in
#' \code{\link{fit_ssl_psdh}}: coefficients are updated given the current
#' per-feature penalty weights derived from the spike-and-slab prior.
#'
#' @param penalty_weights Numeric vector of length \eqn{p}. Per-feature LASSO
#'   penalty weights (inverse expected prior scale \eqn{1 / E[S_j]}), as
#'   returned by \code{\link{expected_penalty_weights}}.  Passed as
#'   \code{penalty.factor} to \code{fastCrrp}.
#' @param timing_vector Numeric vector of length \eqn{n}. Observed
#'   event/censoring times.
#' @param status_vector Integer vector of length \eqn{n}. Event status codes;
#'   see \code{cencode_num} and \code{failcode_num} for the coding convention.
#' @param feature_matrix Numeric matrix of dimensions \eqn{n \times p}.
#'   Design/feature matrix (no intercept).
#' @param cencode_num Integer. Value in \code{status_vector} that denotes a
#'   censored observation. Default \code{0}.
#' @param failcode_num Integer. Value in \code{status_vector} that denotes
#'   the event of interest. Default \code{1}.
#' @param lambda Numeric. Global LASSO tuning parameter.  Default
#'   \code{mean(penalty_weights) / length(status_vector)}.
#'
#' @returns A \code{fastCrrp} model object as returned by
#'   \code{fastcmprsk::fastCrrp}, containing (among other fields)
#'   \code{$coef} (coefficient estimates), \code{$logLik} (log-likelihood at
#'   convergence), \code{$breslowJump} (baseline hazard increments), and
#'   \code{$uftime} (unique failure times).
#'
#' @seealso \code{\link{fit_ssl_psdh}}, \code{\link{expected_penalty_weights}}
#'
#' @importFrom fastcmprsk fastCrrp
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pw  <- rep(1, ncol(x))
#' mod <- update_betas(pw, y[, 1], y[, 2], x)
#' mod$coef
#' }
update_betas <- function(
  penalty_weights,
  timing_vector,
  status_vector,
  feature_matrix,
  cencode_num = 0,
  failcode_num = 1,
  lambda = mean(penalty_weights) / length(status_vector)
) {
  fastcmprsk::fastCrrp(
    Crisk(
      timing_vector,
      status_vector,
      cencode = cencode_num,
      failcode = failcode_num
    ) ~ feature_matrix,
    penalty.factor = penalty_weights,
    penalty = "LASSO",
    lambda = lambda
  )
}
