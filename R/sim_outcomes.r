#' Simulate Weibull proportional-hazards event times
#'
#' Wraps \code{\link[simsurv]{simsurv}} to generate individual event times
#' from a Weibull proportional-hazards model.  The linear predictor
#' \eqn{\eta_i} enters the model with a fixed coefficient of 1, so the
#' person-specific hazard is
#' \deqn{h_i(t) = \gamma \lambda t^{\gamma - 1} \exp(\eta_i).}
#'
#' @param shape Numeric scalar. Weibull shape parameter \eqn{\gamma > 0}.
#'   \eqn{\gamma = 1} gives a constant (exponential) hazard;
#'   \eqn{\gamma > 1} an increasing hazard; \eqn{\gamma < 1} a decreasing
#'   hazard.
#' @param scale Numeric scalar. Weibull rate parameter \eqn{\lambda > 0}.
#'   Higher values produce shorter event times.  Note this is the rate
#'   parameterisation: the corresponding \pkg{stats} scale is
#'   \eqn{\lambda^{-1/\gamma}}.
#' @param eta_val Numeric vector of length \eqn{n}. Pre-computed linear
#'   predictors \eqn{\eta_i = x_i^\top \beta} for each subject.
#'
#' @return A data frame with one row per subject as returned by
#'   \code{\link[simsurv]{simsurv}}: columns \code{id} (integer) and
#'   \code{eventtime} (numeric).
#'
#' @importFrom simsurv simsurv
#'
#' @seealso \code{\link{sim_generate}}, \code{\link[simsurv]{simsurv}}
#'
#' @export
sim_weibull_ph <- function(shape, scale, eta_val) {
  x_frame     <- data.frame(eta = eta_val)
  betas_frame <- c(eta = 1)

  simsurv::simsurv(dist    = "weibull",
                   lambdas = scale,
                   gammas  = shape,
                   x       = x_frame,
                   betas   = betas_frame)
}


#' Generate censoring times
#'
#' Produces a vector of \eqn{n} latent censoring times according to the
#' censoring specification in \code{cspec} (the \code{censor} sub-list of a
#' \code{\link{sim_spec}} object).  Four censoring types are supported:
#'
#' \describe{
#'   \item{\code{"none"}}{Returns \code{maxtime + 10} for every subject —
#'     effectively no censoring before the last observed event.}
#'   \item{\code{"administrative"}}{Returns \code{cspec$administrative_time}
#'     for every subject.}
#'   \item{\code{"random"}}{Draws from the distribution specified by
#'     \code{cspec$random$dist} (currently only \code{"exponential"} is
#'     implemented, with rate \code{cspec$random$rate}).  Times exceeding
#'     \code{cspec$administrative_time} are capped at that value.}
#'   \item{\code{"administrative_plus_random"}}{Same as \code{"random"}: draws
#'     from the random distribution and caps at \code{administrative_time}.}
#' }
#'
#' @param cspec Named list; the \code{censor} sub-list of a
#'   \code{\link{sim_spec}} object.  Used fields: \code{type},
#'   \code{administrative_time}, \code{random$dist}, \code{random$rate}.
#' @param n Integer. Number of subjects.
#' @param maxtime Numeric. Maximum observed latent event time across all
#'   causes; used only when \code{cspec$type = "none"} to push the censoring
#'   time safely past all events.
#'
#' @return Numeric vector of length \eqn{n} containing the latent censoring
#'   time for each subject.
#'
#' @importFrom stats rexp
#'
#' @seealso \code{\link{sim_spec}}, \code{\link{sim_generate}}
#'
#' @export
sim_generate_censor_time <- function(cspec, n, maxtime) {
  if (cspec$type == "none") {
    return(rep(maxtime + 10, n))
  }

  if (cspec$type == "administrative") {
    return(rep(cspec$administrative_time, n))
  }

  latent_censor_time <- switch(
    cspec$random$dist,
    exponential = {
      temp <- rexp(n, rate = cspec$random$rate)
      temp[temp >= cspec$administrative_time] <- cspec$administrative_time
      temp
    },
    stop("Unknown random censor distribution: ", cspec$random$dist)
  )

  latent_censor_time
}
