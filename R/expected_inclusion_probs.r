#' E-step: posterior inclusion probabilities
#'
#' Computes \eqn{E[\gamma_j \mid \beta_j, \pi, s_0, s_1]} for each feature
#' \eqn{j} — the posterior probability that feature \eqn{j} belongs to the
#' slab (active) component of the spike-and-slab prior.  Uses Bayes' theorem
#' with Laplace spike (\code{s0}) and slab (\code{s1}) densities and global
#' mixture weight \code{pi}.
#'
#' @param s1 Numeric scalar. Slab scale parameter of the Laplace prior
#'   (larger value; controls the spread of the active-feature distribution).
#' @param s0 Numeric scalar. Spike scale parameter of the Laplace prior
#'   (small value; concentrates inactive features near zero).
#' @param pi Numeric scalar or vector broadcastable to length \eqn{p}.
#'   Current global mixture probability — the prior probability that any
#'   given feature is active (\eqn{\pi = \Pr(\gamma_j = 1)}).
#' @param betas Numeric vector of length \eqn{p}. Current coefficient
#'   estimates \eqn{\hat{\beta}} from the previous M-step.
#'
#' @returns Numeric vector of length \eqn{p}. Posterior inclusion
#'   probability \eqn{E[\gamma_j \mid \beta_j]} for each feature; values
#'   are in \eqn{[0, 1]}.
#'
#' @seealso \code{\link{fit_ssl_psdh}}, \code{\link{expected_penalty_weights}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' betas <- c(0.5, 0.0, -0.3, 0.0, 0.1)
#' expected_inclusion_probs(s1 = 0.5, s0 = 0.04, pi = 0.1, betas = betas)
#' }
expected_inclusion_probs <- function(s1, s0, pi, betas) {

  #p(betaj | gammaj = 1, s1)
  dens_Slab <- dlaplace(betas, mu = 1, b = s1)
  #p(betaj | gammaj = 0, s0)
  dens_Spike <- dlaplace(betas, mu = 0, b = s0)

  #p(gammaj = 1 | pi)
  #p(gammaj = 0 | pi)
  prior_Slab <- pi
  prior_Spike <- 1-pi

  dens_Slab * prior_Slab / (dens_Spike * prior_Spike + dens_Slab * prior_Slab)
}

