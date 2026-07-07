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
#' @param pips Numeric vector of length \eqn{p}.
#'   Set of current prior probability that any
#'   given feature is active (\eqn{\pi = \Pr(\gamma_j = 1)}).
#' @param betas Numeric vector of length \eqn{p}. Current coefficient
#'   estimates \eqn{\hat{\beta}} from the previous M-step.
#' @param exact If false, overall mixing probability is substituted for the calculated leave one out mixing probabilities. If true, individual leave one out mixing probabilities used.
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
expected_inclusion_probs <- function(s0, s1, pips, betas,
                                     exact = FALSE,
                                     a = 1, b = 1) {


  #if(exact == FALSE){
  #p(betaj | gammaj = 1, s1)
  dens_Slab <- dlaplace(betas, mu = 0, b = s1)
  #p(betaj | gammaj = 0, s0)
  dens_Spike <- dlaplace(betas, mu = 0, b = s0)

  #p(gammaj = 1 | pi)
  #p(gammaj = 0 | pi)

  # if(exact){
  # prior_Slab <- leave_one_out_mean(pips)
  # prior_Spike <- 1-leave_one_out_mean(pips)
  # }else{
  #   prior_Slab <- mean(pips)
  #   prior_Spike <- 1-mean(pips)
  # }

  p <- length(pips)

  if(exact){
    # leave-one-out Beta(a, b) mode: drop coordinate j from the numerator
    # sum and from the coefficient count in the denominator
    prior_Slab  <- (sum(pips) - pips + a - 1) / ((p - 1) + a + b - 2)
  }else{
    # global Beta(a, b) mode: (sum(pips) + a - 1) / (p + a + b - 2)
    prior_Slab  <- (sum(pips) + a - 1) / (p + a + b - 2)
  }
  prior_Spike <- 1 - prior_Slab

  ret <- dens_Slab * prior_Slab / (dens_Spike * prior_Spike + dens_Slab * prior_Slab)
  # }else if(exact == TRUE){
  #
  # thetas <- leave_one_out_mean(pips)
  #   ret <- NULL
  #   for(x in 1:length(betas)){
  #     denom <- 1 + (1-thetas[x])/thetas[x] *
  #       (s1/s0)*exp(-abs(betas[x])*(1/s0-1/s1))
  #     ret[x]<- 1/denom
  #   }
  # }
    ret
}



