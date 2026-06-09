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
#' @param init_lam_path Numeric sequence vector. Path to define global shrinkage parameter for the initializing model.
#' @param fixed_global_shrinkage NULL or numeric. If not null, the global
#'   shrinkage parameter is coerced to the supplied value for each step in the
#'   model, AND the (expensive) \code{cv_fastCrrp} search for the initializing
#'   model is skipped: the supplied value is used directly as the initial
#'   lambda. The tuned initial lambda (whether searched or supplied) is
#'   returned in \code{$tuned_lambda}, so callers can capture it from one fit
#'   and reuse it via \code{fixed_global_shrinkage} on subsequent fits to avoid
#'   re-running the search.
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
#'     \item{\code{$tuned_lambda}}{Numeric. The initial global shrinkage
#'       lambda used to seed the EM algorithm: the \code{cv_fastCrrp}
#'       \code{lambda_min} when \code{fixed_global_shrinkage} is \code{NULL},
#'       otherwise the supplied \code{fixed_global_shrinkage}. Intended for
#'       caching and reuse across fits with the same data folds.}
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
                         epsilon=1e-04,
                         init_lam_path = 10^seq(log10(0.1),
                                                log10(0.001),
                                                length = 25),
                         fixed_global_shrinkage = NULL){

  .dedupe_warnings({

    ss0 <- ss[1]
    ss1 <- ss[2]

    #Initial penalty weights
     ##expected inclusion probs all equal to initial sparsity
    current_mixture_prob <- initial_sparsity
    current_inclusion_probs <- rep(initial_sparsity, ncol(x))
    current_penalty_weights <- expected_penalty_weights(ss1,
                                                        ss0,
                                                        current_inclusion_probs)



    #Initial model with non-adaptive lasso
    # When a single fixed global shrinkage value is supplied, skip the
    # expensive cv_fastCrrp search and use it directly as the initial lambda.
    # A vector fixed_global_shrinkage (e.g. a lambda path) preserves the old
    # behaviour: the search still seeds the initial model.
    if(is.null(fixed_global_shrinkage) || length(fixed_global_shrinkage) > 1){
      init_mod_search <- cv_fastCrrp(x, y[,1], y[,2], k = 5,
                                     penalty = "LASSO",
                                     lambda_path = init_lam_path,
                                     tuning = "wolbers",
                                     eval_quantile = 0.5)
      current_lambda <- init_mod_search$lambda_min
    } else {
      current_lambda <- fixed_global_shrinkage
    }

    # The tuned (or supplied) initial lambda; returned so it can be cached
    # and reused as fixed_global_shrinkage on subsequent fits.
    tuned_lambda <- current_lambda

    init_mod <- fastcmprsk::fastCrrp(
      Crisk(
        y[,1],
        y[,2],
        cencode = 0,
        failcode = 1
      ) ~ x,
      penalty = "LASSO",
      lambda = current_lambda
    )
    current_betas <- init_mod$coef


    #Initial likelihood
    devold <- 0

    outer_convergence <- FALSE
    iterations <- 0
    for(iter in 1:maxit){



      # E-Step
      # Update inclusion probabilites (gamma_j)


      current_inclusion_probs <- expected_inclusion_probs(ss1, ss0,
                                                        current_inclusion_probs,
                                                        current_betas)



      # Update penalty weights (inverse S_j)
      current_penalty_weights <- expected_penalty_weights(ss1,
                                                       ss0,
                                                       current_inclusion_probs)
      penalties <- ifelse(abs(current_penalty_weights) < 1e-10,
                          1e-10,
                          current_penalty_weights)




      # M-Step
      # Update mixture prob (pi)
      current_mixture_prob <- mean(current_inclusion_probs)


      # Update shrinkage (lambda)
      current_lambda <- 1/(ncol(x)*nrow(x))
      if(!is.null(fixed_global_shrinkage)){
        current_lambda <- fixed_global_shrinkage}



      # Update regression coefficients (betas)
      mod <- update_betas(penalties,
                          y[,1], y[,2], x,
                          cencode_num = 0,
                          failcode_num = 1,
                          lambda = current_lambda)

      current_betas <- mod$coef

      logLik <- mod$logLik

      # # Convergence check (using log-likelihood)

      if(iter > 5 && abs(logLik - devold)/(0.1 + abs(logLik)) < epsilon) {
        outer_convergence <- TRUE
        break
      }
      devold <- logLik
      iterations <- iterations+1


    }#end iter loop

    coefficients_df <- data.frame("Variable" = colnames(x),
                                  "Estimate" = mod$coef)

    ret <- list()
    ret$x <- x
    ret$y <- y
    ret$coefficients <- coefficients_df
    ret$penalty.factor <- penalties
    ret$lambda <- mod$lambda.path
    ret$tuned_lambda <- tuned_lambda
    ret$ss <- ss
    ret$conv <- outer_convergence
    ret$iterations <- iterations
    ret$pips <- current_inclusion_probs
    ret$final_model_object <-mod

    return(ret)

  })

}
