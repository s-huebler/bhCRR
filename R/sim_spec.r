#' Build a simulation specification
#'
#' Constructs a fully-specified \code{sim_spec} object that controls every
#' aspect of a bhCRR simulation: design matrix generation, variable selection
#' truth, competing-risk baseline hazards, and censoring.  User-supplied
#' sub-lists are merged recursively into the defaults via
#' \code{\link{sim_merge_lists}}, so you only need to specify the fields you
#' want to override.
#'
#' @param n Integer. Number of observations (samples). Default \code{50}.
#' @param p Integer. Number of features (columns of the design matrix).
#'   Default \code{100}.
#' @param x Named list of design-matrix options merged into the defaults
#'   below.  Recognised keys:
#'   \describe{
#'     \item{\code{family}}{Character. Distribution family for generating
#'       covariates. One of \code{"gaussian"} (default), \code{"binary_latent"},
#'       \code{"zip"}, or \code{"mixed"}.}
#'     \item{\code{corr}}{Named list with keys \code{type} (one of
#'       \code{"indep"} (default), \code{"ar1"}, \code{"block"},
#'       \code{"factor"}) and \code{rho} (correlation strength,
#'       default \code{0.0}).}
#'     \item{\code{zero_inflation}}{Numeric in \code{[0, 1)}. Proportion of
#'       entries zeroed out after family generation. Default \code{0.0}.}
#'     \item{\code{standardize}}{Logical. Whether to column-standardize
#'       \code{X} to zero mean and unit variance. Default \code{TRUE}.}
#'     \item{\code{binary_prob}}{Numeric. Marginal probability of 1 when
#'       \code{family = "binary_latent"}. Default \code{0.5}.}
#'     \item{\code{zip_rate}}{Numeric. Base Poisson rate for
#'       \code{family = "zip"}. Default \code{1.0}.}
#'     \item{\code{factor_k}}{Integer. Number of latent factors when
#'       \code{corr$type = "factor"}. Default \code{5}.}
#'   }
#' @param truth Named list of variable-selection truth options merged into
#'   the defaults below.  Recognised keys:
#'   \describe{
#'     \item{\code{sparsity}}{Numeric. Either a proportion in \code{(0, 1)}
#'       (fraction of \code{p} active predictors per cause) or an integer
#'       \code{>= 1} (absolute count). Default \code{0.05}.}
#'     \item{\code{share_active_set}}{Integer. Number of predictors shared
#'       between cause 1 and cause 2 active sets. Default \code{0}.}
#'     \item{\code{shared_multiplier}}{Numeric scalar or vector of length
#'       \code{share_active_set}. Ratio of cause-1 to cause-2 beta for each
#'       shared predictor. Default \code{1}.}
#'     \item{\code{beta}}{Named list describing the beta prior:
#'       \code{dist} (\code{"normal"}), \code{mean} (\code{0}),
#'       \code{sd} (length-2 numeric giving per-cause standard deviation,
#'       default \code{c(0.3, 0.3)}).}
#'     \item{\code{target_eta_sd}}{Length-2 numeric. Desired standard
#'       deviation of the linear predictor \eqn{\eta = X\beta} for each
#'       cause. Betas are rescaled post-draw to hit this target exactly.
#'       Default \code{c(0.5, 0.5)}.}
#'   }
#' @param risks Named list of baseline hazard options.  The default contains
#'   two causes (\code{cause1}, \code{cause2}); additional causes can be
#'   added.  Each cause is itself a named list with:
#'   \describe{
#'     \item{\code{baseline}}{Named list with \code{dist} (\code{"weibull"}),
#'       \code{shape} (Weibull shape \eqn{\gamma}, default \code{1}), and
#'       \code{scale} (Weibull rate \eqn{\lambda}, default \code{0.011}).}
#'     \item{\code{link}}{Character. Link function. Currently only
#'       \code{"PH"} (proportional hazards) is implemented.}
#'   }
#' @param censor Named list of censoring options merged into the defaults
#'   below.  Recognised keys:
#'   \describe{
#'     \item{\code{type}}{Character. One of \code{"none"},
#'       \code{"administrative"} (default), \code{"random"}, or
#'       \code{"administrative_plus_random"}.}
#'     \item{\code{administrative_time}}{Numeric. Study end time used for
#'       administrative censoring. Default \code{65}.}
#'     \item{\code{random}}{Named list with \code{dist}
#'       (\code{"exponential"}) and \code{rate} (default \code{0.15}).
#'       Only used when \code{type} includes random censoring.}
#'   }
#' @param seed Integer or \code{NULL}. Random seed stored in the spec and
#'   forwarded to \code{\link{sim_generate}}. Default \code{NULL}.
#'
#' @return A named list of class \code{"sim_spec"} with components:
#'   \describe{
#'     \item{\code{n}}{Integer. Number of observations.}
#'     \item{\code{p}}{Integer. Number of features.}
#'     \item{\code{x}}{Named list. Fully-resolved design-matrix options.}
#'     \item{\code{truth}}{Named list. Fully-resolved variable-selection
#'       truth options.}
#'     \item{\code{risks}}{Named list. Fully-resolved per-cause baseline
#'       hazard options.}
#'     \item{\code{censor}}{Named list. Fully-resolved censoring options.}
#'     \item{\code{seed}}{Integer or \code{NULL}.}
#'   }
#'
#' @seealso \code{\link{sim_generate}}, \code{\link{sim_merge_lists}}
#'
#' @examples
#' # Default specification
#' spec <- sim_spec()
#'
#' # Override sample size and use AR(1) correlation
#' spec2 <- sim_spec(
#'   n = 200,
#'   p = 50,
#'   x = list(corr = list(type = "ar1", rho = 0.6)),
#'   truth = list(sparsity = 0.1, target_eta_sd = c(0.3, 0.3)),
#'   seed = 42
#' )
#'
#' @export
sim_spec <- function(
  n = 50,
  p = 100,
  x = list(),
  truth = list(),
  risks = list(),
  censor = list(),
  seed = NULL
) {
  defaults <- list(
    n = n,
    p = p,
    x = list(
      family = "gaussian",
      corr = list(type = "indep", rho = 0.0),
      zero_inflation = 0.0,
      standardize = TRUE,
      binary_prob = 0.5,
      zip_rate = 1.0,
      factor_k = 5
    ),
    truth = list(
      share_active_set = 0,
      shared_multiplier = 1,
      sparsity = 0.05,
      beta = list(dist = "normal", mean = 0, sd = c(0.3, 0.3)),
      target_eta_sd = c(0.5, 0.5)
    ),
    risks = list(
      cause1 = list(
        baseline = list(dist = "weibull", shape = 1, scale = 0.011),
        link = "PH"
      ),
      cause2 = list(
        baseline = list(dist = "weibull", shape = 1, scale = 0.011),
        link = "PH"
      )
    ),
    censor = list(
      type = "administrative",
      administrative_time = 65,
      random = list(dist = "exponential", rate = 0.15)
    ),
    seed = seed
  )

  spec <- defaults
  spec$x      <- sim_merge_lists(defaults$x,      x)
  spec$truth  <- sim_merge_lists(defaults$truth,  truth)
  spec$risks  <- sim_merge_lists(defaults$risks,  risks)
  spec$censor <- sim_merge_lists(defaults$censor, censor)

  class(spec) <- "sim_spec"

  spec
}
