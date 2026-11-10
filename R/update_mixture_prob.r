#' M-step: update the global mixture probability
#'
#' Updates the global spike-and-slab mixture probability \eqn{\pi} as the
#' mean of the current posterior inclusion probabilities.  This is the
#' closed-form M-step solution for \eqn{\pi} under a flat Beta(1,1) prior.
#'
#' @param p Numeric vector of length \eqn{p}. Posterior inclusion
#'   probabilities \eqn{E[\gamma_j \mid \beta_j]} for all features, as
#'   returned by \code{\link{expected_inclusion_probs}}.
#' @param a Numeric. Shape parameter \eqn{a} of the \eqn{Beta(a, b)} prior.
#'   Default \code{1}.
#' @param b Numeric. Shape parameter \eqn{b} of the \eqn{Beta(a, b)} prior.
#'   Default \code{1}.
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
update_mixture_prob <- function(p, a = 1, b = 1){
  #mean(p)
  (sum(p) + a - 1) / (length(p) + a + b - 2)
  }

#' Count-based update for the global mixture probability
#'
#' Computes \eqn{\theta^{(t+1)} = (a + \|\hat{\beta}\|_0) / (a + b + p)},
#' the coordinate-ascent spike-and-slab LASSO update of Rockova and George
#' (2018, JASA 113:431--444, equation for \eqn{\pi}).  This is a closed-form
#' approximation to \eqn{E[\theta \mid \hat{\beta}]} based on the observed
#' support size \eqn{\|\hat{\beta}\|_0 = \sum_j \mathbf{1}(\hat{\beta}_j \neq 0)}.
#'
#' @param betas Numeric vector of length \eqn{p}. Current coefficient
#'   estimates \eqn{\hat{\beta}} from the previous M-step.
#' @param a Numeric. Shape parameter \eqn{a} of the \eqn{Beta(a, b)} prior
#'   on \eqn{\theta}. Default \code{1}.
#' @param b Numeric. Shape parameter \eqn{b} of the \eqn{Beta(a, b)} prior
#'   on \eqn{\theta}. Default \code{1}.
#'
#' @returns Numeric scalar in \eqn{(0, 1)}. Updated global mixture
#'   probability \eqn{\theta}.
#'
#' @seealso \code{\link{update_mixture_prob}}, \code{\link{fit_ssl_psdh}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' betas <- c(1.2, 0, -0.5, 0, 0.3)
#' update_mixture_prob_count(betas, a = 1, b = length(betas))
#' }
update_mixture_prob_count <- function(betas, a = 1, b = 1) {
  (a + sum(betas != 0)) / (a + b + length(betas))
}
