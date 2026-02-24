#' expected_penalty_weights
#'
#' @param s1
#' @param s0
#' @param p
#'
#' @returns
#' @export
#'
#' @examples
expected_penalty_weights <- function(s1, s0, p){
  (1 - p) / s0 + p / s1
}
