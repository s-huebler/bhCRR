#' E-step: expected penalty weights for the M-step LASSO
#'
#' Computes the per-feature penalty weights \eqn{1 / E[S_j]} used in the
#' LASSO M-step of \code{\link{fit_ssl_psdh}}.  Each weight is the
#' posterior expectation of the inverse prior scale:
#' \deqn{\frac{1}{E[S_j]} = \frac{1 - p_j}{s_0} + \frac{p_j}{s_1},}
#' where \eqn{p_j} is the posterior inclusion probability for feature
#' \eqn{j}.  Features with low inclusion probability receive a large penalty
#' (spike-like shrinkage); features with high inclusion probability receive a
#' small penalty (slab-like shrinkage).
#'
#' @param s1 Numeric scalar. Slab scale parameter of the Laplace prior.
#' @param s0 Numeric scalar. Spike scale parameter of the Laplace prior.
#' @param p Numeric vector of length \eqn{p}. Posterior inclusion
#'   probabilities \eqn{E[\gamma_j \mid \beta_j]}, as returned by
#'   \code{\link{expected_inclusion_probs}}.
#'
#' @returns Numeric vector of length \eqn{p}. Per-feature penalty weights
#'   passed as \code{penalty.factor} to the LASSO solver in
#'   \code{\link{update_betas}}.
#'
#' @seealso \code{\link{fit_ssl_psdh}}, \code{\link{expected_inclusion_probs}},
#'   \code{\link{update_betas}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' probs <- expected_inclusion_probs(s1 = 0.5, s0 = 0.04, pi = 0.1,
#'                                   betas = c(0.5, 0.0, -0.3))
#' expected_penalty_weights(s1 = 0.5, s0 = 0.04, p = probs)
#' }
expected_penalty_weights <- function(s1, s0, p){
  (1 - p) / s0 + p / s1
}

