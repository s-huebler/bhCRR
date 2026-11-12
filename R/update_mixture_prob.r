#' Update the global mixture probability
#'
#' Updates the global spike-and-slab mixture probability \eqn{\theta} under a
#' \eqn{Beta(a, b)} prior.  One of three rules is applied:
#'
#' \describe{
#'   \item{\code{"count"}}{
#'     \eqn{\theta^{(t+1)} = (a + \|\hat{\beta}\|_0) / (a + b + p)} from
#'     Rockova and George (2018, JASA 113:431--444). A closed-form
#'     approximation to \eqn{E[\theta \mid \hat{\beta}]} based on the
#'     observed support size; not an EM maximizer.
#'   }
#'   \item{\code{"expected"}}{
#'     \eqn{\theta^{(t+1)} = (\sum_j p_j + a - 1) / (p + a + b - 2)}, the
#'     EM M-step for \eqn{\theta} under the \eqn{Beta(a, b)} prior: the
#'     posterior mode given the current posterior inclusion probabilities
#'     \eqn{p_j}.
#'   }
#'   \item{\code{"expected_clamp"}}{
#'     The \code{"expected"} update floored at \code{0.001}. Analogous to
#'     \code{BhGLM::update.ptheta.group()}, which clamps \eqn{\theta} to
#'     \eqn{[0.01, 0.99]}; here the floor is \code{0.001} and there is no
#'     upper clamp.
#'   }
#' }
#'
#' @param pips Numeric vector of length \eqn{p}. Posterior inclusion
#'   probabilities \eqn{E[\gamma_j \mid \beta_j]}, as returned by
#'   \code{\link{expected_inclusion_probs}}. Required for \code{"expected"}
#'   and \code{"expected_clamp"}; ignored by \code{"count"}.
#' @param betas Numeric vector of length \eqn{p}. Current coefficient
#'   estimates \eqn{\hat{\beta}} from the previous M-step. Required for
#'   \code{"count"}; ignored by the expected rules.
#' @param a Numeric. Shape parameter \eqn{a} of the \eqn{Beta(a, b)} prior.
#'   Default \code{1}.
#' @param b Numeric. Shape parameter \eqn{b} of the \eqn{Beta(a, b)} prior.
#'   Default \code{1}.
#' @param theta_update Character string, one of \code{"count"} (default),
#'   \code{"expected"}, or \code{"expected_clamp"}. Selects the update rule.
#'
#' @returns Numeric scalar in \eqn{(0, 1)}. Updated global mixture
#'   probability \eqn{\theta}.
#'
#' @seealso \code{\link{fit_ssl_psdh}}, \code{\link{expected_inclusion_probs}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' betas <- c(1.2, 0, -0.5, 0, 0.3)
#' update_mixture_prob(betas = betas, a = 1, b = length(betas),
#'                     theta_update = "count")
#'
#' probs <- expected_inclusion_probs(s0 = 0.04, s1 = 0.5,
#'                                   betas = rnorm(100), theta = 0.1)
#' update_mixture_prob(pips = probs, theta_update = "expected")
#' update_mixture_prob(pips = probs, theta_update = "expected_clamp")
#' }
update_mixture_prob <- function(pips = NULL, betas = NULL, a = 1, b = 1,
                                theta_update = c("count", "expected", "expected_clamp")) {
  theta_update <- match.arg(theta_update)
  switch(theta_update,
    count = {
      if (is.null(betas)) stop("theta_update = 'count' requires 'betas'.")
      (a + sum(betas != 0)) / (a + b + length(betas))
    },
    expected = {
      if (is.null(pips)) stop("theta_update = 'expected' requires 'pips'.")
      (sum(pips) + a - 1) / (length(pips) + a + b - 2)
    },
    expected_clamp = {
      if (is.null(pips)) stop("theta_update = 'expected_clamp' requires 'pips'.")
      max((sum(pips) + a - 1) / (length(pips) + a + b - 2), 0.001)
    }
  )
}
