#' CIF computation from a fitted ssl_psdh model
#'
#' Internal workhorse called by \code{\link{predict.ssl_psdh}}.  Forms the
#' linear predictor from \code{object$final_model_object$coef}, accumulates
#' the Breslow baseline cumulative hazard up to \code{prediction_time}, and
#' returns subject-specific absolute risk via
#' \eqn{1 - \exp(-\hat{\Lambda}_0(t) \exp(\hat{\eta}_i))}.
#'
#' @param object A fitted \code{"ssl_psdh"} object.
#' @param newx Numeric matrix, \eqn{m \times p}.
#' @param prediction_time Numeric scalar.
#'
#' @return Numeric vector of length \eqn{m}.
#' @keywords internal
.predict_cif <- function(object, newx, prediction_time) {

  # 1) Calculate linear predictor
  beta <- as.vector(object$final_model_object$coef)
  lp <- as.vector(newx %*% beta)

  # 2) Extract baseline cumulative hazard
  # from fastCrrp: breslowJump (hazard increments) and uftime (event times)
  base_haz_jumps <- object$final_model_object$breslowJump[, 2]
  failure_times <- object$final_model_object$uftime

  # Lambda_0(t) = Sum of jumps where time <= prediction_time
  valid_indices <- which(failure_times <= prediction_time)

  if (length(valid_indices) == 0) {
    cum_base_haz <- 0
  } else {
    cum_base_haz <- sum(base_haz_jumps[valid_indices])
  }

  # 3) Calculate Absolute Risk (CIF)

  predicted_risk <- 1 - exp(- cum_base_haz * exp(lp))

  return(predicted_risk)
}


#' Predict cumulative incidence from a fitted ssl_psdh model
#'
#' Computes predicted absolute risk (cumulative incidence function, CIF)
#' at a specified time for new subjects.  The linear predictor is formed from
#' the fitted coefficients and the Breslow baseline cumulative hazard is
#' integrated up to \code{prediction_time} to produce subject-specific
#' CIF values via \eqn{1 - \exp(-\hat{\Lambda}_0(t) \exp(\hat{\eta}_i))}.
#'
#' @param object A fitted \code{"ssl_psdh"} object returned by
#'   \code{\link{fit_ssl_psdh}}.  Must contain
#'   \code{$final_model_object$coef} (coefficient vector),
#'   \code{$final_model_object$breslowJump} (matrix of baseline hazard
#'   increments, column 2), and \code{$final_model_object$uftime}
#'   (unique failure times).
#' @param newx Numeric matrix of dimensions \eqn{m \times p}.  Design matrix
#'   for the \eqn{m} subjects to predict; columns must correspond to the same
#'   features used during fitting.
#' @param prediction_time Numeric scalar.  Time horizon at which to evaluate
#'   the cumulative incidence; only baseline hazard jumps at or before this
#'   time are summed.
#' @param ... Unused; required for compatibility with the
#'   \code{\link[stats]{predict}} generic.
#'
#' @returns Numeric vector of length \eqn{m}.  Predicted cumulative incidence
#'   (absolute risk) in \eqn{[0, 1]} for each row of \code{newx} at
#'   \code{prediction_time}.
#'
#' @seealso \code{\link{fit_ssl_psdh}}, \code{\link{wolbers_c}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fit  <- fit_ssl_psdh(x_train, y_train, init_method = "LASSO_cv")
#' risk <- predict(fit, newx = x_test, prediction_time = 30)
#' hist(risk)
#' }
predict.ssl_psdh <- function(object, newx, prediction_time, ...) {
  if (!inherits(object, "ssl_psdh"))
    stop("'object' must be of class \"ssl_psdh\"")
  .predict_cif(object, newx, prediction_time)
}
