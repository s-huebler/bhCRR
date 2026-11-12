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
#' @param initial_sparsity Numeric in \code{(0, 1)}. The initial mixing
#'   probability \eqn{\theta^{(0)}}, used in the first E-step under both
#'   \code{theta_update} rules. Default \code{0.5}, following Rockova and
#'   George (2018).
#' @param theta_a Numeric. Shape parameter \eqn{a} of the \eqn{Beta(a, b)}
#'   prior on the global mixture probability \eqn{\theta}. Default \code{1}.
#' @param theta_b Numeric. Shape parameter \eqn{b} of the \eqn{Beta(a, b)}
#'   prior on \eqn{\theta}. The default \code{ncol(x)} (i.e. \eqn{p}) gives
#'   the paper's non-separable \eqn{Beta(1, p)} calibration, which drives
#'   \eqn{\theta} toward the true sparsity \eqn{q/p} and adaptively shrinks
#'   noise coefficients. Set \code{theta_b = 1} to recover the flat
#'   \eqn{Beta(1, 1)} (separable-like) behaviour.
#' @param theta_update Character string, one of \code{"count"} (default) or
#'   \code{"expected"}. Controls how the global mixture probability
#'   \eqn{\theta} is updated after each M-step.
#'   \code{"count"}: the coordinate-ascent update
#'   \eqn{\theta^{(t+1)} = (a + \|\hat{\beta}^{(t)}\|_0) / (a + b + p)} from
#'   Rockova and George (2018), based on the observed support size.
#'   \code{"expected"}: the EM update
#'   \eqn{\theta^{(t+1)} = (\sum_j p_j^{(t)} + a - 1) / (p + a + b - 2)},
#'   based on the current posterior inclusion probabilities.
#'   The updated \eqn{\theta^{(t+1)}} is used as the prior inclusion
#'   probability in the next E-step.
#' @param maxit Integer. Maximum number of EM iterations. Default \code{50}.
#' @param epsilon Numeric. Convergence threshold. After at least 5 EM
#'   iterations the loop stops when (a)
#'   \eqn{|\text{dev}_t - \text{dev}_{t-1}| / (0.1 + |\text{dev}_t|) < \epsilon},
#'   where \eqn{\text{dev} = -2 \times \log} pseudo-likelihood (matches
#'   \code{BhGLM::bmlasso}), AND (b) no coefficient moved between zero and
#'   nonzero, AND (c) among coefficients nonzero in both iterations,
#'   \eqn{\max_j |\beta_j^t - \beta_j^{t-1}| / |\beta_j^{t-1}| < \epsilon};
#'   condition (c) is vacuously satisfied when that shared-support set is
#'   empty. Default \code{1e-04}.
#' @param init Length-\eqn{p} numeric vector of initial coefficient values.
#'   Exactly one of \code{init} or \code{init_method} must be supplied.
#' @param init_lam_path Defunct. Passing any value other than \code{NULL}
#'   raises an error. Pass \code{init_args = list(lambda_path = ...)} instead.
#' @param inner_maxit_start Integer. Starting value for the inner
#'   \code{fastCrrp} iteration budget (\code{max.iter}) used in the M-step.
#'   If the inner solver fails to converge, this is escalated by 200 per
#'   refit up to a ceiling of 20000; if it still has not converged the EM
#'   stops and \code{conv} is returned as \code{FALSE}. The escalated value
#'   is carried across EM iterations. Default \code{1000}.
#' @param init_method Character string naming a built-in initialization method
#'   (\code{"LASSO_cv"}, \code{"LASSO_bic"}, \code{"zero"}), or a function
#'   with signature \code{function(x, y, ...)}, or a string naming such a
#'   function in the calling environment. Exactly one of \code{init} or
#'   \code{init_method} must be supplied.
#' @param init_args Named list of extra arguments forwarded to the
#'   \code{init_method} function via \code{do.call}. Ignored when \code{init}
#'   is supplied. Use this to override method defaults, e.g.
#'   \code{init_args = list(lambda_path = my_path, k = 10)}.
#'
#' @returns A list with the following fields:
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
#'     \item{\code{$conv}}{Logical. \code{TRUE} if the EM converged within
#'       \code{maxit} iterations.}
#'     \item{\code{$iterations}}{Integer. Number of completed EM iterations.}
#'     \item{\code{$pips}}{Numeric vector of final posterior inclusion
#'       probabilities.}
#'     \item{\code{$init}}{Numeric vector. The initial coefficient vector
#'       actually used to start the EM.}
#'     \item{\code{$init_method}}{Character. Label of the initialization used:
#'       \code{"supplied"} when \code{init} was passed directly, the built-in
#'       name (\code{"LASSO_cv"}, \code{"LASSO_bic"}, \code{"zero"}), or
#'       \code{"custom"} for a user-supplied function.}
#'     \item{\code{$init_meta}}{List. Provenance returned by the
#'       initialization function (e.g. selected lambda, CV scores, elapsed
#'       time). \code{NULL} when \code{init} was supplied directly.}
#'     \item{\code{$theta}}{Numeric scalar. The final value of the global
#'       mixture probability \eqn{\theta} at convergence.}
#'     \item{\code{$theta_path}}{Numeric vector of length \code{iterations + 1}.
#'       The trajectory of \eqn{\theta}: first entry is \code{initial_sparsity};
#'       subsequent entries are the values after each M-step update.}
#'     \item{\code{$theta_update}}{Character. The \code{theta_update} rule used:
#'       \code{"count"} or \code{"expected"}.}
#'   }
#'
#' @seealso \code{\link{predict.ssl_psdh}}, \code{\link{update_betas}},
#'   \code{\link{expected_inclusion_probs}},
#'   \code{\link{expected_penalty_weights}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Built-in CV-LASSO initialization (count-based theta update, default)
#' fit <- fit_ssl_psdh(x, y, ss = c(0.04, 0.5), init_method = "LASSO_cv")
#' fit$coefficients
#'
#' # Supply a precomputed starting vector
#' init_vec <- rep(0, ncol(x))
#' fit2 <- fit_ssl_psdh(x, y, ss = c(0.04, 0.5), init = init_vec)
#'
#' # User-defined initialization function
#' custom_init <- function(x, y, ...) rep(0, ncol(x))
#' fit3 <- fit_ssl_psdh(x, y, ss = c(0.04, 0.5), init_method = custom_init)
#' }
fit_ssl_psdh <- function(x, y,
                         ss = c(0.04, 0.5),
                         initial_sparsity = 0.5,
                         theta_a = 1,
                         theta_b = ncol(x),
                         theta_update = c("count", "expected"),
                         maxit = 50,
                         epsilon = 1e-04,
                         init = NULL,
                         init_lam_path = NULL,
                         inner_maxit_start = 1000,
                         init_method = NULL,
                         init_args = list()) {

  # Defunct parameter
  if (!is.null(init_lam_path))
    stop("init_lam_path is defunct; pass init_args = list(lambda_path = ...) instead.")

  # Initialization contract: exactly one of init / init_method must be supplied
  both_given    <- !is.null(init) && !is.null(init_method)
  neither_given <- is.null(init)  && is.null(init_method)
  if (both_given)
    stop("Supply either 'init' or 'init_method', not both.")
  if (neither_given)
    stop(
      "fit_ssl_psdh() now requires an initialization: supply either ",
      "init = <length-p numeric> or init_method = one of 'LASSO_cv', ",
      "'LASSO_bic', 'zero' (or your own function)."
    )

  theta_update <- match.arg(theta_update)

  # Standard input checks
  if (!is.numeric(initial_sparsity) || length(initial_sparsity) != 1 ||
      initial_sparsity <= 0 || initial_sparsity >= 1)
    stop("initial_sparsity must be a single number strictly between 0 and 1.")
  if (maxit < 1) stop("maxit must be positive integer.")
  if (epsilon <= 0) stop("epsilon must be a positive number.")
  if (length(ss) != 2) stop("ss must be a vector of length 2")
  if (ss[1] > ss[2]) {
    ss <- sort(ss)
    print("ss warning, scale values supplied out of order. Spike scale taken as ss[1] and slab scale taken as ss[2].")
  }
  if (nrow(x) != nrow(y)) stop("x and y must have the same number of rows")

  # Capture caller's frame before entering .dedupe_warnings so that
  # .resolve_init_method can find user-named functions in the calling scope.
  caller_env <- parent.frame()

  .dedupe_warnings({

    ss0 <- ss[1]
    ss1 <- ss[2]

    #Initial theta (theta^(0)); used by the first E-step under both rules.
    current_theta <- initial_sparsity
    theta_path <- initial_sparsity
    current_inclusion_probs <- rep(initial_sparsity, ncol(x))

    #Initial estimates
    if (!is.null(init)) {
      res <- .validate_init(init, ncol(x), "supplied")
      init_label <- "supplied"
    } else {
      resolved  <- .resolve_init_method(init_method, envir = caller_env)
      raw       <- do.call(resolved$fn, c(list(x = x, y = y), init_args))
      res       <- .validate_init(raw, ncol(x), resolved$label)
      init_label <- resolved$label
    }
    current_betas <- res$init

    betas_path <- data.frame("Initial" = current_betas)
    pips_path <- data.frame("Initial" = current_inclusion_probs)

    inner_maxit <- inner_maxit_start


    #Initial deviance
    previous_dev <- 0
    outer_convergence <- FALSE
    iterations <- 0
    rel_dev_change     <- NA_real_
    n_support_changes  <- NA_integer_
    rel_beta_change    <- NA_real_
    for(iter in 1:maxit){
      previous_betas <- current_betas



      # E-Step
      # Update inclusion probabilites (gamma_j)
      current_inclusion_probs <- expected_inclusion_probs(ss0, ss1,
                                                          current_betas,
                                                          theta = current_theta)

      # Update penalty weights (inverse S_j)
      current_penalty_weights <- expected_penalty_weights(ss0,
                                                       ss1,
                                                       current_inclusion_probs)
      penalties <- ifelse(abs(current_penalty_weights) < 1e-10,
                          1e-10,
                          current_penalty_weights)



      # M-Step
      # Update shrinkage (lambda)
      current_lambda <- 1/nrow(x)

      # Update regression coefficients (betas)
      # Fit with the current inner iteration budget.
      mod <- update_betas(penalties,
                          y[,1], y[,2], x,
                          cencode_num = 0,
                          failcode_num = 1,
                          lambda = current_lambda,
                          getBreslowJumps = FALSE,
                          max.iter = inner_maxit)

      # If fastCrrp did not converge, escalate the inner iteration budget by
      # 200 and refit, up to a ceiling of 20000. inner_maxit is defined
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
                            getBreslowJumps = FALSE,
                            max.iter = inner_maxit)
      }

      current_betas <- mod$coef

      if (theta_update == "count") {
        current_theta <- update_mixture_prob_count(current_betas, theta_a, theta_b)
      } else {
        current_theta <- update_mixture_prob(current_inclusion_probs, theta_a, theta_b)
      }
      theta_path <- c(theta_path, current_theta)

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
      ##relative deviance change is the necessary condition; coefficient stability is the sufficient condition

      #necessary condition
      rel_dev_change <- abs(dev - previous_dev) / (0.1 + abs(dev))
      if (iter > 5 && rel_dev_change < epsilon) {
        #sufficient condition
        nz_prev <- previous_betas != 0
        nz_cur  <- current_betas  != 0
        n_support_changes <- sum(nz_prev != nz_cur)
        both_nz <- nz_prev & nz_cur
        rel_beta_change <- if (any(both_nz)) {
          max(abs(current_betas[both_nz] - previous_betas[both_nz]) /
              abs(previous_betas[both_nz]))
        } else 0

        if (n_support_changes == 0L && rel_beta_change < epsilon) {
          outer_convergence <- TRUE
          break
        }
      }
      previous_dev <- dev
      iterations <- iterations+1


    }#end iter loop

    if(!outer_convergence){
      warning(sprintf(
        "fit_ssl_psdh did not converge within maxit = %d EM iterations; consider increasing maxit or relaxing epsilon. last relative deviance change = %g; support changes = %d; max relative coefficient change on shared support = %g",
        maxit, rel_dev_change, n_support_changes, rel_beta_change))
    }

    # One final refit at the converged state with Breslow jumps enabled.
    # This is the only call that computes breslowJump / uftime; all EM
    # M-steps above ran with getBreslowJumps = FALSE to avoid the n x p
    # allocation on every iteration.
    final_mod <- update_betas(penalties,
                              y[,1], y[,2], x,
                              cencode_num = 0,
                              failcode_num = 1,
                              lambda = current_lambda,
                              getBreslowJumps = TRUE,
                              max.iter = inner_maxit)

    coefficients_df <- data.frame("Variable" = colnames(x),
                                  "Estimate" = as.numeric(mod$coef))

    ret <- list()
    ret$x <- x
    ret$y <- y
    ret$coefficients <- coefficients_df
    ret$penalty.factor <- penalties
    ret$ss <- ss
    ret$conv <- outer_convergence
    ret$iterations <- iterations
    ret$pips <- current_inclusion_probs
    ret$final_model_object <- final_mod
    ret$betas_path <- betas_path
    ret$pips_path <- pips_path
    ret$init        <- res$init
    ret$init_method <- init_label
    ret$init_meta   <- res$meta
    ret$theta       <- current_theta
    ret$theta_path  <- theta_path
    ret$theta_update <- theta_update

    class(ret) <- "ssl_psdh"
    return(ret)

  })

}


