#' Bayesian Bootstrap Spike-and-Slab LASSO for competing risks (BB-SSL)
#'
#' Approximate posterior sampling for the penalised subdistribution hazard
#' (Fine-Gray) spike-and-slab model via the Bayesian Bootstrap Spike-and-Slab
#' LASSO of Nie & Rockova (2023, JASA).  This is a competing-risks adaptation
#' of \pkg{BBSSL}'s \code{BB_SSL}: the inner optimiser (their \code{SSLASSO_2},
#' a Gaussian-likelihood SSL solver) is replaced by \code{\link{fit_ssl_psdh}},
#' which maximises the Fine-Gray partial likelihood under a Laplace
#' spike-and-slab prior.
#'
#' @details
#' Implements Algorithm 2 (BB-SSL Sampling) of Nie & Rockova (2023).  For each
#' of \code{NSample} draws:
#' \enumerate{
#'   \item \strong{(a) sample weights} \eqn{w_t \sim \pi(w)}: drawn as
#'     \eqn{w \sim n \times \mathrm{Dirichlet}(\alpha, \dots, \alpha)} via
#'     normalised \code{Gamma(alpha, 1)} variates.
#'   \item \strong{(b) sample jitter} \eqn{\mu_t}: each \eqn{\mu_j} is drawn
#'     i.i.d. from the spike \eqn{\psi_0(\mu)}, a mean-zero Laplace with scale
#'     \eqn{s_0 = 1/\lambda_0}.
#'   \item \strong{(c) optimise} \eqn{\tilde\beta_t} from the reweighted,
#'     recentred objective (their Eq. 13).
#' }
#'
#' \strong{Competing-risks adaptations.}  The original BB-SSL implements the
#' reweighting (\eqn{X^* = \sqrt{w}\,X}) and prior recentring
#' (\eqn{Y^* = \sqrt{w}(Y - X\mu)}, then \eqn{\tilde\beta = \hat\beta^* + \mu})
#' as a Gaussian least-squares trick.  That trick does not transfer to the
#' Fine-Gray partial likelihood, and \code{fastcmprsk::fastCrrp} (called inside
#' \code{\link{fit_ssl_psdh}}) exposes neither case weights nor an offset.  Two
#' substitutions are therefore used:
#' \itemize{
#'   \item \strong{Weights via the Bayesian bootstrap.}  Rather than row-scaling
#'     the design matrix, the weights \eqn{w} define a resampling distribution
#'     \eqn{w / n} over the \eqn{n} observations; a size-\eqn{n} resample of
#'     rows is drawn and \code{\link{fit_ssl_psdh}} is refit on it.  This is the
#'     classical Bayesian-bootstrap reading of \eqn{\pi(w)} and leaves
#'     \code{fit_ssl_psdh} untouched as a black-box optimiser.
#'   \item \strong{Jitter as a post-hoc location shift.}  \eqn{\mu_t} is sampled
#'     from the spike as in step (b) and added back to the fitted coefficients
#'     (\eqn{\tilde\beta_t = \hat\beta_t + \mu_t}), approximating the paper's
#'     recentre-then-add-back step since the partial likelihood cannot be
#'     recentred directly.  Set \code{jitter = FALSE} to drop step (b) and run a
#'     pure weighted Bayesian bootstrap.
#' }
#'
#' @param x Numeric matrix \eqn{n \times p}.  Feature/design matrix, no
#'   intercept column.
#' @param y Two-column numeric matrix \eqn{n \times 2}: column 1 event/censoring
#'   times, column 2 status codes (\code{0} censored, \code{1} event of
#'   interest, \code{2} competing event).  Same convention as
#'   \code{\link{fit_ssl_psdh}}.
#' @param NSample Integer.  Number of posterior samples to generate.
#' @param ss Length-2 numeric \code{c(s0, s1)}: spike and slab \emph{scale}
#'   parameters of the Laplace prior, passed straight to
#'   \code{\link{fit_ssl_psdh}}.  Must satisfy \code{s1 > s0 > 0}.
#'   Default \code{c(0.04, 0.5)}.
#' @param lambda Length-2 numeric \code{c(lambda0, lambda1)}: spike and slab
#'   \emph{rate} (penalty) parameters, \code{lambda0 >> lambda1}.  Used only for
#'   (i) sampling the jitter \eqn{\mu_j \sim} Laplace(scale \code{1/lambda0})
#'   and (ii) the \eqn{\gamma} inclusion threshold.  Defaults to \code{1 / ss}
#'   (i.e. \code{lambda0 = 1/s0}, \code{lambda1 = 1/s1}), the rate parameterised
#'   form of \code{ss}; supply it explicitly to decouple the jitter/threshold
#'   scale from the fitting scale.
#' @param alpha Numeric > 0.  Dirichlet concentration for the bootstrap weights;
#'   \eqn{w \sim n \times \mathrm{Dir}(\alpha, \dots, \alpha)}.  Larger
#'   \code{alpha} => weights closer to uniform.  Default \code{3}.
#' @param theta Numeric in \code{(0, 1)}.  Prior mixing proportion used in the
#'   \eqn{\gamma} threshold.  Default \code{0.5}.
#' @param jitter Logical.  If \code{TRUE} (default) apply the step-(b) location
#'   shift; if \code{FALSE} run a pure weighted Bayesian bootstrap.
#' @param discard Logical.  If \code{TRUE}, redraw any sample whose inner
#'   \code{fit_ssl_psdh} EM did not converge (uses the \code{$conv} flag).
#'   Default \code{FALSE}.
#' @param progress Logical.  Print a progress bar.  Default \code{TRUE}.
#' @param ... Further arguments forwarded to \code{\link{fit_ssl_psdh}}
#'   (e.g. \code{initial_sparsity}, \code{maxit}, \code{epsilon},
#'   \code{init_lam_path}, \code{fixed_global_shrinkage}).
#'
#' @returns A list with:
#'   \describe{
#'     \item{\code{beta}}{\code{NSample x p} matrix of sampled coefficients.}
#'     \item{\code{gamma}}{\code{NSample x p} 0/1 matrix of inclusion
#'       indicators from thresholding each \code{beta} draw.}
#'     \item{\code{inclusion_prob}}{Length-\eqn{p} vector of marginal inclusion
#'       probabilities (column means of \code{gamma}).}
#'     \item{\code{post_mean}}{Length-\eqn{p} posterior mean of \code{beta}.}
#'     \item{\code{ss}, \code{lambda}, \code{alpha}}{Settings as used.}
#'   }
#'
#' @seealso \code{\link{fit_ssl_psdh}}, \code{\link{tune_ssl_psdh}},
#'   \code{\link{expected_inclusion_probs}}
#'
#' @references
#' Nie, L. & Rockova, V. (2023). Bayesian Bootstrap Spike-and-Slab LASSO.
#' \emph{Journal of the American Statistical Association}, 118(543), 2013-2028.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- bb_ssl(x, y, NSample = 100, ss = c(0.04, 0.5))
#' fit$inclusion_prob
#' }
bb_ssl <- function(x, y,
                   NSample,
                   ss = c(0.04, 0.5),
                   lambda = 1 / ss,
                   alpha = 3,
                   theta = 0.5,
                   jitter = TRUE,
                   discard = FALSE,
                   progress = TRUE,
                   ...) {

  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(x)
  p <- ncol(x)

  lambda0 <- lambda[1]   # spike rate (large)
  lambda1 <- lambda[2]   # slab  rate (small)
  if (lambda0 <= lambda1) {
    warning("Expected lambda0 (spike) >> lambda1 (slab); got lambda = c(",
            lambda0, ", ", lambda1, ").")
  }

  # Mean-zero Laplace sampler with scale b  (difference of two Exp(rate = 1/b))
  rlaplace0 <- function(m, b) {
    stats::rexp(m, rate = 1 / b) - stats::rexp(m, rate = 1 / b)
  }

  beta  <- matrix(NA_real_, nrow = NSample, ncol = p)
  gamma <- matrix(0L,       nrow = NSample, ncol = p)

  pb <- if (progress) utils::txtProgressBar(min = 0, max = NSample, style = 3)

  i <- 1L
  while (i <= NSample) {

    # (a) sample weights  w ~ n * Dirichlet(alpha, ..., alpha)
    w <- stats::rgamma(n, shape = alpha, rate = 1)
    w <- w * (n / sum(w))

    # (b) sample jitter   mu_j ~ Laplace(0, scale = 1/lambda0)  [the spike]
    mu <- if (jitter) rlaplace0(p, b = 1 / lambda0) else rep(0, p)

    # (c) optimise tilde-beta_t.
    #     Bayesian-bootstrap resample of rows with probabilities w / n,
    #     then refit fit_ssl_psdh (replacing SSLASSO_2) on the resample.
    idx     <- sample.int(n, size = n, replace = TRUE, prob = w / n)
    fit     <- fit_ssl_psdh(x[idx, , drop = FALSE],
                            y[idx, , drop = FALSE],
                            ss = ss, ...)

    if (isTRUE(discard) && isFALSE(fit$conv)) next  # redraw if EM didn't converge

    beta_hat <- fit$coefficients$Estimate
    beta[i, ] <- beta_hat + mu                      # post-hoc location shift

    # threshold beta -> gamma  (posterior inclusion under the SSL prior)
    #   p* = 1 / (1 + (1-theta)/theta * (lambda0/lambda1) * exp(-(lambda0-lambda1)|beta|))
    odds <- (1 - theta) / theta *
            (lambda0 / lambda1) * exp(-(lambda0 - lambda1) * abs(beta[i, ]))
    pstar <- 1 / (1 + odds)
    pstar[is.na(pstar)] <- 1
    gamma[i, ] <- as.integer(pstar > 0.5)

    if (progress) {
      utils::setTxtProgressBar(pb, i)
      if (i == NSample) close(pb)
    }
    i <- i + 1L
  }

  colnames(beta) <- colnames(gamma) <-
    if (is.null(colnames(x))) paste0("V", seq_len(p)) else colnames(x)

  list(beta           = beta,
       gamma          = gamma,
       inclusion_prob = colMeans(gamma),
       post_mean      = colMeans(beta),
       ss             = ss,
       lambda         = c(lambda0 = lambda0, lambda1 = lambda1),
       alpha          = alpha)
}
