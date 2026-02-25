#' fit_ssl_psdh
#'
#' @param x
#' @param y
#' @param ss
#' @param initial_sparsity
#' @param maxit
#' @param epsilon
#'
#' @returns
#' @export
#'
#' @examples
#'
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
                           cencode = 0,
                           failcode = 1,
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
                        cencode = 0,
                        failcode = 1,
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



