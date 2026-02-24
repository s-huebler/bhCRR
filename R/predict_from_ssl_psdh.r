#' predict_from_ssl_psdh
#'
#' @param object
#' @param newx
#' @param prediction_time
#'
#' @returns
#' @export
#'
#' @examples
predict_from_ssl_psdh <- function(object, newx, prediction_time) {

  # 1) Calculate linear predictor
  beta <- as.vector(object$coef)
  lp <- as.vector(newx %*% beta)

  # 2) Extract baseline cumulative hazard
  # from fastCrrp: breslowJump (hazard increments) and uftime (event times)
  base_haz_jumps <- object$breslowJump[, 2]
  failure_times <- object$uftime

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
