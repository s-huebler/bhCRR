#' Fit a spike-and-slab lasso competing-risks model
#'
#' Estimates regression coefficients for cause-specific competing risks via
#' an EM algorithm that alternates between updating posterior inclusion
#' probabilities (E-step) and refitting a penalised subdistribution hazard
#' model (M-step).  The penalty uses a Laplace spike-and-slab prior
#' parameterised by scale pair \code{ss = c(s0, s1)}.
#'
#' @param x Numeric matrix of dimensions \eqn{n \times p}. Feature/design
#'   matrix (no intercept column).
#' @param y Two-column numeric matrix of dimensions \eqn{n \times 2}.
#'   Column 1 contains observed event/censoring times; column 2 contains
#'   event status codes (\code{0} = censored, \code{1} = event of interest,
#'   \code{2} = competing event).
#' @param ss Length-2 numeric vector \code{c(s0, s1)}. Spike (\code{s0})
#'   and slab (\code{s1}) scale parameters of the Laplace prior.
#'   Must satisfy \code{s1 > s0 > 0}. Default \code{c(0.04, 0.5)}.
#' @param initial_sparsity Numeric in \code{(0, 1)}. Starting value for the
#'   global mixture probability (prior proportion of active features).
#'   Default \code{0.05}.
#' @param maxit Integer. Maximum number of EM iterations. Default \code{50}.
#' @param epsilon Numeric. Convergence threshold: iteration stops when the
#'   relative change in log-likelihood falls below this value (after at
#'   least 5 iterations). Default \code{1e-04}.
#'
#' @returns A \code{fastCrrp} model object augmented with additional fields:
#'   \describe{
#'     \item{\code{$x}}{The feature matrix supplied as \code{x}.}
#'     \item{\code{$y}}{The outcome matrix supplied as \code{y}.}
#'     \item{\code{$coefficients}}{Data frame with columns \code{Variable}
#'       and \code{Estimate} giving the final coefficient for each feature.}
#'     \item{\code{$penalty.factor}}{Numeric vector of final per-feature
#'       penalty weights from the last EM iteration.}
#'     \item{\code{$lambda}}{Numeric. The LASSO tuning parameter used at
#'       convergence.}
#'     \item{\code{$ss}}{The \code{ss} argument as supplied.}
#'   }
#'
#' @seealso \code{\link{tune_ssl_psdh}}, \code{\link{cv_ssl_psdh}},
#'   \code{\link{predict_from_ssl_psdh}}, \code{\link{update_betas}},
#'   \code{\link{expected_inclusion_probs}},
#'   \code{\link{expected_penalty_weights}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- fit_ssl_psdh(x, y, ss = c(0.04, 0.5), initial_sparsity = 0.05)
#' fit$coefficients
#' }
fit_ssl_psdh <- function(x, y,
                         ss=c(0.04, 0.5),
                         initial_sparsity = 0.05,
                         maxit = 50,
                         epsilon=1e-04){

  ss0 <- ss[1]
  ss1 <- ss[2]

  #Initial penalty weights
  current_mixture_prob <- rep(initial_sparsity, nrow(x))
  init_mixture_scale <- (1-current_mixture_prob)*ss0+current_mixture_prob*ss1

  current_penalty_weights <- 1/init_mixture_scale


  #Initial betas
  init_lambda <- mean(current_penalty_weights)/length(y[,2])/ncol(x)
  #print(init_lambda)
  init_mod <- update_betas(penalty_weights = current_penalty_weights,
                           timing_vector = y[,1],
                           status_vector = y[,2],
                           feature_matrix = x,
                           cencode_num = 0,
                           failcode_num = 1,
                           lambda = init_lambda)
  current_betas <- init_mod$coef
  #print(current_betas)

  #Initial devold
  devold <- 0

  for(iter in 1:maxit){
    # E-Step
    # Update inclusion probabilites (gamma_j)


    current_inclusion_probs <- suppressWarnings(expected_inclusion_probs(ss1, ss0,
                                                        current_mixture_prob,
                                                        current_betas))

    # Update penalty weights (inverse S_j)
    current_penalty_weights <- expected_penalty_weights(ss0, ss1,
                                                        current_inclusion_probs)

    # M-Step
    # Update mixture prob (pi)
    current_mixture_prob <- mean(current_inclusion_probs)

    # Update betas
    Pf <- ifelse(abs(current_penalty_weights) < 1e-10, 1e-10, current_penalty_weights)
    current_lambda <- sum(Pf)/(nrow(x) * ncol(x))
    mod <- update_betas(Pf,
                        y[,1], y[,2], x,
                        cencode_num = 0,
                        failcode_num = 1,
                        lambda = current_lambda)






    # # Convergence check (using log-likelihood)
    logLik <- mod$logLik
    if(abs(logLik - devold)/(0.1 + abs(logLik)) < epsilon & iter > 5) {
      conv <- TRUE
      break
    }
    devold <- logLik
  }

  coefficients_df <- data.frame("Variable" = colnames(x),
                                "Estimate" = mod$coef)


  mod$x <- x
  mod$y <- y
  mod$coefficients <- coefficients_df
  mod$penalty.factor <- Pf
  mod$lambda <- mod$lambda.path
  mod$ss <- ss

  return(mod)

}



