#' update_betas
#'
#' @param penalty_weights
#' @param timing_vector
#' @param status_vector
#' @param feature_matrix
#' @param cencode
#' @param failcode
#' @param lambda
#'
#' @returns
#' @export
#'
#' @examples
update_betas <- function(penalty_weights,
                         timing_vector,
                         status_vector,
                         feature_matrix,
                         cencode_num = 0,
                         failcode_num = 1,
                         lambda = mean(penalty_weights)/length(status_vector)
){
  fastcmprsk::fastCrrp(Crisk(timing_vector, status_vector,
                             cencode = cencode_num,
                             failcode = failcode_num) ~ feature_matrix,
                       penalty.factor = penalty_weights,
                       penalty = "LASSO",
                       lambda = lambda)
}
