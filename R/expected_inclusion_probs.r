#' expected_inclusion_probs
#'
#' @param s1
#' @param s0
#' @param pi
#' @param betas
#'
#' @returns
#' @export
#'
#' @examples
expected_inclusion_probs <- function(s1, s0, pi, betas) {

  #p(betaj | gammaj = 1, s1)
  dens_Slab <- dlaplace(betas, mu = 1, b = s1)
  #p(betaj | gammaj = 0, s0)
  dens_Spike <- dlaplace(betas, mu = 0, b = s0)

  #p(gammaj = 1 | pi)
  #p(gammaj = 0 | pi)
  prior_Slab <- pi
  prior_Spike <- 1 - pi

  dens_Slab * prior_Slab / (dens_Spike * prior_Spike + dens_Slab * prior_Slab)
}
