#' E-step: posterior inclusion probabilities
#'
#' Computes \eqn{E[\gamma_j \mid \beta_j, \theta, s_0, s_1]} for each feature
#' \eqn{j} — the posterior probability that feature \eqn{j} belongs to the
#' slab (active) component of the spike-and-slab prior.  Uses Bayes' theorem
#' with Laplace spike (\code{s0}) and slab (\code{s1}) densities and prior
#' inclusion probability \code{theta}.
#'
#' @param s0 Numeric scalar. Spike scale parameter of the Laplace prior
#'   (small value; concentrates inactive features near zero).
#' @param s1 Numeric scalar. Slab scale parameter of the Laplace prior
#'   (larger value; controls the spread of the active-feature distribution).
#' @param betas Numeric vector of length \eqn{p}. Current coefficient
#'   estimates \eqn{\hat{\beta}} from the previous M-step.
#' @param theta Numeric. Prior inclusion probability \eqn{\Pr(\gamma_j = 1)}.
#'   Either a single value (shared across features) or a vector of length
#'   \code{length(betas)} (per-feature). All values must be strictly between
#'   0 and 1.
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
#' expected_inclusion_probs(s0 = 0.04, s1 = 0.5, betas = betas, theta = 0.1)
#' }
expected_inclusion_probs <- function(s0, s1, betas, theta) {

  if (!is.numeric(theta) ||
      !(length(theta) == 1L || length(theta) == length(betas)) ||
      any(theta <= 0) || any(theta >= 1))
    stop("'theta' must be numeric of length 1 or length(betas), ",
         "with all values strictly between 0 and 1.")

  dens_Slab  <- dlaplace(betas, mu = 0, b = s1)
  dens_Spike <- dlaplace(betas, mu = 0, b = s0)

  theta * dens_Slab / (theta * dens_Slab + (1 - theta) * dens_Spike)
}
