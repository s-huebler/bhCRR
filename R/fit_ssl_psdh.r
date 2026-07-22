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
#' @param theta_a Numeric. Shape parameter \eqn{a} of the \eqn{Beta(a, b)}
#'   prior on the global mixture probability \eqn{\theta}. Default \code{1}.
#' @param theta_b Numeric. Shape parameter \eqn{b} of the \eqn{Beta(a, b)}
#'   prior on \eqn{\theta}. The default \code{ncol(x)} (i.e. \eqn{p}) gives
#'   the paper's non-separable \eqn{Beta(1, p)} calibration, which drives
#'   \eqn{\theta} toward the true sparsity \eqn{q/p} and adaptively shrinks
#'   noise coefficients. Set \code{theta_b = 1} to recover the flat
#'   \eqn{Beta(1, 1)} (separable-like) behaviour.
#' @param maxit Integer. Maximum number of EM iterations. Default \code{50}.
#' @param epsilon Numeric. Convergence threshold: iteration stops when the
#'   relative change in log-likelihood falls below this value (after at
#'   least 5 iterations). Default \code{1e-04}.
#' @param init Length-p numeric vector. Provide initial estimates. Otherwise leave as null to get initial estimates from non-adaptive LASSO.
#' @param init_lam_path Numeric sequence vector. Path to define LASSO shrinkage parameter for the initializing model. Required if init is NULL.
#' @param inner_maxit_start Integer. Starting value for the inner
#'   \code{fastCrrp} iteration budget (\code{max.iter}) used in the M-step.
#'   If the inner solver fails to converge, this is escalated by 200 per
#'   refit up to a ceiling of 2000; if it still has not converged the EM
#'   stops and \code{conv} is returned as \code{FALSE}. The escalated value
#'   is carried across EM iterations. Default \code{1000}.

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
                         theta_a = 1,
                         theta_b = ncol(x),
                         maxit = 50,
                         epsilon=1e-04,
                         init = NULL,
                         init_lam_path = 10^seq(log10(0.1),
                                                log10(0.001),
                                                length = 10),
                         inner_maxit_start = 1000){

  #Requirements (to add once fastCrrp back on CRAN)
  #if (!requireNamespace("fastcmprsk")) install.packages("glmnet")
  #require(fastcmprsk)

  #Input checks
  if(maxit < 1) stop("maxit must be positive integer.")
  if(epsilon <= 0) stop("epsilon must be a positive number.")
  if(length(ss)!=2) stop("ss must be a vector of length 2")
  if(ss[1]>ss[2]){
    ss<- sort(ss)
    print("ss warning, scale values supplied out of order. Spike scale taken as ss[1] and slab scale taken as ss[2].")
  }
  if(!is.null(init)){
    if(!is.numeric(init) | length(init)!=ncol(x)){
      stop("provide initial value to each coefficient (no intercept) or leave NULL")
    }
  }
  if(nrow(x) != nrow(y)) stop("x and y must have the same number of rows")

  .dedupe_warnings({

    ss0 <- ss[1]
    ss1 <- ss[2]

    #Initial penalty weights
     ##expected inclusion probs all equal to initial sparsity
    current_mixture_prob <- initial_sparsity
    current_inclusion_probs <- rep(initial_sparsity, ncol(x))
    current_penalty_weights <- expected_penalty_weights(ss0,
                                                        ss1,
                                                        current_inclusion_probs)



    #Initial estimates
    if(!is.null(init)){
      #If provided
      current_betas <- init

    }else{
      #If init==NULL

      ## to-do: add in overall arguments to the larger function that can be supplied here (tuning type and eval_quantile)
      init_mod_search <- cv_fastCrrp_cpp(x, y[,1], y[,2], k = 5,
                                     penalty = "LASSO",
                                     lambda_path = init_lam_path,
                                     tuning = "wolbers",
                                     eval_quantile = 0.5)
      current_betas <- init_mod_search$full_model$coef[,
        init_mod_search$lambda == init_mod_search$lambda_min]

    betas_path <- data.frame("Initial" = current_betas)
    pips_path <- data.frame("Initial" = current_inclusion_probs)

    inner_maxit <- inner_maxit_start









    #Initial deviance
    previous_dev <- 0
    outer_convergence <- FALSE
    iterations <- 0
    for(iter in 1:maxit){
      previous_betas <- current_betas



      # E-Step
      # Update inclusion probabilites (gamma_j)


      current_inclusion_probs <- expected_inclusion_probs(ss0, ss1,
                                                        current_inclusion_probs,
                                                        current_betas,
                                                        exact = ncol(x)<nrow(x),
                                                        a = theta_a,
                                                        b = theta_b)



      # Update penalty weights (inverse S_j)
      current_penalty_weights <- expected_penalty_weights(ss0,
                                                       ss1,
                                                       current_inclusion_probs)
      penalties <- ifelse(abs(current_penalty_weights) < 1e-10,
                          1e-10,
                          current_penalty_weights)




      # M-Step
      # Update mixture prob (pi) under the Beta(theta_a, theta_b) prior
      current_mixture_prob <- update_mixture_prob(current_inclusion_probs,
                                                  a = theta_a,
                                                  b = theta_b)


      # Update shrinkage (lambda)
      current_lambda <- 1/nrow(x)

      # Update regression coefficients (betas)
      # Fit with the current inner iteration budget.
      mod <- update_betas(penalties,
                          y[,1], y[,2], x,
                          cencode_num = 0,
                          failcode_num = 1,
                          lambda = current_lambda,
                          max.iter = inner_maxit)
                          #lambda = current_lambda)

      # If fastCrrp did not converge, escalate the inner iteration budget by
      # 200 and refit, up to a ceiling of 2000. inner_maxit is defined
      # outside the EM loop, so once escalated it stays escalated for
      # subsequent iterations (which tend to need at least as many inner
      # iterations as the one that first stalled).
      while(mod$converged == 0 && inner_maxit < 20000){
        inner_maxit <- inner_maxit + 200
        mod <- update_betas(penalties,
                            y[,1], y[,2], x,
                            cencode_num = 0,
                            failcode_num = 1,
                            lambda = current_lambda,
                            max.iter = inner_maxit)
      }

      current_betas <- mod$coef

      betas_path[,iter+1] <- current_betas
      pips_path[,iter+1]<- current_inclusion_probs

      # Inner solver hit the iteration ceiling without converging: the
      # M-step estimate can't be trusted, so stop the EM and report
      # non-convergence of the outer loop.
      if(mod$converged == 0){
        outer_convergence <- FALSE
        break
      }



      dev <- -2*mod$logLik

      ##Convergence check
      ##using deviance as necessary and stability as sufficient

      #necessary condition
      if(iter > 5 && abs(dev - previous_dev) < epsilon) {
        #sufficient condition

        #to-do: optimize efficiency
        unstable <- sum(abs(current_betas - previous_betas) > epsilon)

        if(unstable < 1){
        outer_convergence <- TRUE
        break
        }
      }
      previous_dev <- dev
      iterations <- iterations+1


    }#end iter loop

    if(!outer_convergence){
      warning(sprintf(
        "fit_ssl_psdh did not converge within maxit = %d EM iterations; consider increasing maxit or relaxing epsilon.",
        maxit))
    }

    coefficients_df <- data.frame("Variable" = colnames(x),
                                  "Estimate" = mod$coef)

    ret <- list()
    ret$x <- x
    ret$y <- y
    ret$coefficients <- coefficients_df
    ret$penalty.factor <- penalties
    ret$ss <- ss
    ret$conv <- outer_convergence
    ret$iterations <- iterations
    ret$pips <- current_inclusion_probs
    ret$final_model_object <-mod
    ret$betas_path <- betas_path
    ret$pips_path <- pips_path

    return(ret)

  })

}


