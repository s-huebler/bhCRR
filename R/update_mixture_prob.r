#' M-step: update the global mixture probability
#'
#' Updates the global spike-and-slab mixture probability \eqn{\pi} as the
#' mean of the current posterior inclusion probabilities.  This is the
#' closed-form M-step solution for \eqn{\pi} under a flat Beta(1,1) prior.
#'
#' @param p Numeric vector of length \eqn{p}. Posterior inclusion
#'   probabilities \eqn{E[\gamma_j \mid \beta_j]} for all features, as
#'   returned by \code{\link{expected_inclusion_probs}}.
#'
#' @returns Numeric scalar in \eqn{[0, 1]}.  Updated global mixture
#'   probability \eqn{\hat{\pi} = \frac{1}{p} \sum_j E[\gamma_j]}.
#'
#' @seealso \code{\link{fit_ssl_psdh}}, \code{\link{expected_inclusion_probs}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' probs <- expected_inclusion_probs(s1 = 0.5, s0 = 0.04, pi = 0.1,
#'                                   betas = rnorm(100))
#' update_mixture_prob(probs)
#' }
update_mixture_prob <- function(p){mean(p)}
