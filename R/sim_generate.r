#' Generate a simulated competing-risks dataset
#'
#' The top-level simulation driver.  Given a \code{\link{sim_spec}} object,
#' generates the design matrix \eqn{X}, the true active predictor sets and
#' coefficients, latent Weibull event times for each cause, censoring times,
#' and the observed \code{(time, status)} outcome.  Any of the three
#' intermediate objects (\eqn{X}, predictor sets, betas) can be replaced by
#' user-supplied values.
#'
#' @param spec A \code{sim_spec} object created by \code{\link{sim_spec}}, or
#'   a plain list with the same structure (coerced silently).
#' @param custom_design_matrix Optional numeric matrix or data frame
#'   (\eqn{n \times p}) to use in place of the generated design matrix.
#'   \code{NULL} (default) generates \eqn{X} from \code{spec$x}.
#' @param custom_predictor_set Optional length-2 list of integer index
#'   vectors giving the active predictor columns for cause 1 and cause 2
#'   respectively.  Both vectors must contain integers in \code{[1, p]}.
#'   If the validation check fails a warning is issued and the argument is
#'   ignored.  \code{NULL} (default) samples predictors from
#'   \code{spec$truth}.
#' @param custom_betas Optional length-2 list of numeric vectors of
#'   coefficients, parallel to \code{custom_predictor_set}.  Lengths must
#'   match the corresponding predictor sets; a length mismatch triggers a
#'   warning and falls back to generated betas.  \code{NULL} (default)
#'   generates betas via \code{\link{sim_compute_effects}}.
#' @param seed Integer or \code{NULL}.  When non-\code{NULL} this seed takes
#'   precedence over \code{spec$seed}.  \code{NULL} (default) falls back to
#'   \code{spec$seed}; if that is also \code{NULL} the RNG is not reset.
#'
#' @return A named list of class \code{"sim_data"} with components:
#'   \describe{
#'     \item{\code{id}}{Integer vector \code{1:n}.}
#'     \item{\code{time}}{Numeric vector of observed event/censoring times.}
#'     \item{\code{status}}{Integer vector: \code{0} censored, \code{1} cause
#'       1 event, \code{2} cause 2 event.}
#'     \item{\code{latent_time_event1}}{Numeric. Latent cause-1 event times.}
#'     \item{\code{latent_time_event2}}{Numeric. Latent cause-2 event times.}
#'     \item{\code{latent_time_censor}}{Numeric. Latent censoring times.}
#'     \item{\code{X}}{Numeric matrix \eqn{n \times p}. Design matrix used.}
#'     \item{\code{spec_call}}{The \code{sim_spec} object used.}
#'     \item{\code{truth}}{Named list of resolved truth: \code{cause1},
#'       \code{cause2}, \code{shared_predictors}, \code{share_active_set},
#'       \code{shared_multiplier}, \code{sparsity}, and \code{effects}
#'       (output of \code{\link{sim_compute_effects}}).}
#'   }
#'
#' @seealso \code{\link{sim_spec}}, \code{\link{sim_generate_x}},
#'   \code{\link{sim_generate_truth}}, \code{\link{sim_compute_effects}},
#'   \code{\link{sim_weibull_ph}}, \code{\link{sim_generate_censor_time}}
#'
#' @export
sim_generate <- function(spec,
                         custom_design_matrix = NULL,
                         custom_predictor_set = NULL,
                         custom_betas = NULL,
                         seed = NULL) {
  if (!inherits(spec, "sim_spec")) {
    class(spec) <- unique(c("sim_spec", class(spec)))
  }

  if (!is.null(seed)) {
    set.seed(seed)
  } else if (!is.null(spec$seed)) {
    set.seed(spec$seed)
  }

  n <- as.integer(spec$n)
  p <- as.integer(spec$p)

  # Design matrix
  if (!is.null(custom_design_matrix)) {
    X <- data.matrix(custom_design_matrix)
  } else {
    X <- sim_generate_x(spec$x, n, p)
  }

  print("Success in generating X")

  truth <- spec$truth

  # Variable selection truth
  if (!is.null(custom_predictor_set)) {
    if (!(typeof(custom_predictor_set) == "list" &
          length(custom_predictor_set) == 2 &
          sum(custom_predictor_set[[1]] %% 1) == 0 &
          sum(custom_predictor_set[[2]] %% 1) == 0 &
          max(custom_predictor_set[[1]]) <= p &
          max(custom_predictor_set[[2]]) <= p)) {
      warning("custom_predictor_set should be a list of 2 sets of index integers (with max integer <= p). Ignoring custom predictor sets and custom_betas. Generating from truth defaults. ")
      custom_predictor_set <- NULL
      custom_betas <- NULL
    }
  }

  if (!is.null(custom_predictor_set)) {
    truth$shared_multiplier  <- spec$truth$shared_multiplier
    truth$cause1             <- custom_predictor_set[[1]]
    truth$cause2             <- custom_predictor_set[[2]]
    truth$shared_predictors  <- intersect(custom_predictor_set[[1]],
                                          custom_predictor_set[[2]])
    truth$share_active_set   <- length(truth$shared_predictors)
    truth$sparsity           <- min(c(length(truth$cause1), length(truth$cause1))) / p
    truth$beta               <- spec$truth$beta
    truth$target_eta_sd      <- spec$truth$target_eta_sd
  } else {
    truth <- sim_generate_truth(spec$truth, p, seed = seed)
  }

  print("Success in generating truth")

  # Covariate effects
  if (!is.null(custom_betas)) {
    if (!(length(custom_betas[[1]]) == length(truth$cause1) &
          length(custom_betas[[2]]) == length(truth$cause2))) {
      warning("Lengh of custom betas does not match length of predictor sets. Generating from truth defaults.")
      custom_betas <- NULL
    }
  }

  if (!is.null(custom_betas)) {
    effects1 <- list()
    effects2 <- list()

    effects1$betas1    <- custom_betas[[1]]
    effects1$betas1_sd <- sd(effects1$betas1)
    effects1$eta1      <- X[, truth$cause1] %*% custom_betas[[1]]
    effects1$eta1_mean <- mean(effects1$eta1)
    effects1$eta1_sd   <- sd(effects1$eta1)

    effects2$betas2    <- custom_betas[[2]]
    effects2$betas2_sd <- sd(effects2$betas2)
    effects2$eta2      <- X[, truth$cause2] %*% custom_betas[[2]]
    effects2$eta2_mean <- mean(effects2$eta2)
    effects2$eta2_sd   <- sd(effects2$eta2)

    effects <- list(effects1 = effects1, effects2 = effects2)
  } else {
    effects <- sim_compute_effects(X, truth)
  }

  if (truth$share_active_set > 0) {
    i1 <- match(truth$shared_predictors, truth$cause1)
    i2 <- match(truth$shared_predictors, truth$cause2)
    truth$shared_multiplier <- effects$effects1$betas1[i1] / effects$effects2$betas2[i2]
  } else {
    truth$shared_multiplier <- NA
  }

  truth$beta         <- NULL
  truth$target_eta_sd <- NULL
  truth$effects      <- effects

  print("Success in generating truth with outcomes")

  # Latent event times
  b1 <- spec$risks$cause1$baseline
  b2 <- spec$risks$cause2$baseline

  event1 <- sim_weibull_ph(shape   = b1$shape %||% 1,
                           scale   = b1$scale %||% 0.011,
                           eta_val = truth$effects$effects1$eta1)
  T1 <- event1$eventtime

  event2 <- sim_weibull_ph(shape   = b2$shape %||% 1,
                           scale   = b2$scale %||% 0.011,
                           eta_val = truth$effects$effects2$eta2)
  T2 <- event2$eventtime

  print("Success in generating latent event times")

  maxtime <- max(c(T1, T2))

  # Censoring
  C    <- sim_generate_censor_time(spec$censor, n, maxtime)
  time <- pmin(T1, T2, C)

  event1 <- (T1 <= T2) & (T1 <= C)
  event2 <- (T2 <  T1) & (T2 <= C)
  status <- integer(n)
  status[event1] <- 1L
  status[event2] <- 2L

  print("Success in generating censoring")

  out <- list(
    id                  = seq_len(n),
    time                = as.numeric(time),
    status              = as.integer(status),
    latent_time_event1  = T1,
    latent_time_event2  = T2,
    latent_time_censor  = C,
    X                   = X,
    spec_call           = spec,
    truth               = truth
  )

  class(out) <- "sim_data"
  out
}


#' Generate a design matrix
#'
#' Generates an \eqn{n \times p} numeric design matrix according to the
#' options in \code{xspec} (the \code{x} sub-list of a \code{\link{sim_spec}}
#' object).  Generation proceeds in four steps: (1) draw a latent Gaussian
#' matrix with the requested correlation structure; (2) transform to the
#' requested marginal family; (3) optionally apply zero-inflation; (4)
#' optionally column-standardize.
#'
#' @param xspec Named list of design-matrix options; see the \code{x}
#'   parameter of \code{\link{sim_spec}} for all recognised keys.
#' @param n Integer. Number of rows (observations).
#' @param p Integer. Number of columns (features).
#'
#' @return A numeric matrix of dimensions \eqn{n \times p}.
#'
#' @importFrom stats rnorm qnorm rpois runif sd
#'
#' @seealso \code{\link{sim_spec}}, \code{\link{sim_generate}}
#'
#' @export
sim_generate_x <- function(xspec, n, p) {
  corr_type <- xspec$corr$type %||% "indep"
  rho       <- xspec$corr$rho  %||% 0

  # 1) latent Gaussian core with requested correlation structure
  Z <- switch(
    corr_type,
    indep = matrix(stats::rnorm(n * p), nrow = n),
    ar1 = {
      if (!is.finite(rho) || abs(rho) >= 1)
        stop("For ar1, corr$rho must satisfy abs(rho) < 1.")
      eps <- matrix(stats::rnorm(n * p), nrow = n)
      out <- eps
      if (p >= 2) {
        s <- sqrt(1 - rho^2)
        for (j in 2:p) out[, j] <- rho * out[, j - 1] + s * eps[, j]
      }
      out
    },
    block = {
      block_size <- as.integer(xspec$corr$block_size %||% 50L)
      if (block_size < 1)
        stop("For block correlation, corr$block_size must be >= 1.")
      out <- matrix(0, nrow = n, ncol = p)
      nb  <- ceiling(p / block_size)
      for (b in seq_len(nb)) {
        j0 <- (b - 1) * block_size + 1
        j1 <- min(b * block_size, p)
        m  <- j1 - j0 + 1
        z_shared <- stats::rnorm(n)
        eps      <- matrix(stats::rnorm(n * m), nrow = n)
        out[, j0:j1] <- sqrt(max(rho, 0)) * z_shared +
                        sqrt(max(1 - max(rho, 0), 0)) * eps
      }
      out
    },
    factor = {
      k      <- max(1L, as.integer(xspec$factor_k %||% 5L))
      Zf     <- matrix(stats::rnorm(n * k), nrow = n)
      Lambda <- matrix(stats::rnorm(k * p), nrow = k) / sqrt(k)
      E      <- matrix(stats::rnorm(n * p), nrow = n)
      Zf %*% Lambda + E
    },
    stop("Unknown corr$type: ", corr_type)
  )

  # 2) transform latent Z to requested marginal family
  family <- xspec$family %||% "gaussian"
  X <- switch(
    family,
    gaussian = Z,
    binary_latent = {
      prob <- min(max(xspec$binary_prob %||% 0.5, 1e-6), 1 - 1e-6)
      thr  <- stats::qnorm(1 - prob)
      (Z > thr) * 1.0
    },
    zip = {
      rate0 <- xspec$zip_rate %||% 1.0
      mu    <- log(pmax(rate0, 1e-8)) + 0.3 * Z
      stats::rpois(n * p, lambda = exp(mu)) |> matrix(nrow = n)
    },
    mixed = {
      # 70% Gaussian, 30% binary_latent (first pass)
      p_bin   <- as.integer(round(0.3 * p))
      idx_bin <- if (p_bin > 0) sample.int(p, p_bin) else integer(0)
      Xmix    <- Z
      if (length(idx_bin) > 0) Xmix[, idx_bin] <- (Z[, idx_bin] > 0) * 1.0
      Xmix
    },
    stop("Unknown x$family: ", family)
  )

  # 3) zero-inflation overlay
  zi <- xspec$zero_inflation %||% 0
  if (zi > 0) {
    mask    <- matrix(stats::runif(n * p) < zi, nrow = n)
    X[mask] <- 0
  }

  # 4) optional column standardization
  if (isTRUE(xspec$standardize)) {
    cm      <- colMeans(X)
    cs      <- apply(X, 2, stats::sd)
    cs[cs == 0] <- 1
    X <- sweep(X, 2, cm, "-")
    X <- sweep(X, 2, cs, "/")
  }

  X
}


#' Sample the true active predictor sets
#'
#' Determines which columns of \eqn{X} are truly associated with each
#' competing cause, respecting the requested sparsity level and shared active
#' set size.  The active sets are drawn once and stored so that downstream
#' calls to \code{\link{sim_compute_effects}} produce matching coefficients.
#'
#' @param trspec Named list; the \code{truth} sub-list of a
#'   \code{\link{sim_spec}} object.  Used fields: \code{sparsity},
#'   \code{share_active_set}, \code{shared_multiplier}, \code{beta},
#'   \code{target_eta_sd}.
#' @param p Integer. Total number of features (columns of \eqn{X}).
#' @param seed Integer or \code{NULL}.  When \code{NULL} the fixed fallback
#'   seed \code{9134} is used to ensure reproducible predictor sampling even
#'   without an explicit seed.
#'
#' @return A named list with components:
#'   \describe{
#'     \item{\code{cause1}}{Integer vector of active column indices for
#'       cause 1.}
#'     \item{\code{cause2}}{Integer vector of active column indices for
#'       cause 2.}
#'     \item{\code{shared_predictors}}{Integer vector of column indices
#'       present in both \code{cause1} and \code{cause2}.}
#'     \item{\code{share_active_set}}{Integer. Number of shared predictors
#'       (copied from \code{trspec}).}
#'     \item{\code{shared_multiplier}}{Numeric. Copied from
#'       \code{trspec$shared_multiplier}; interpreted as the ratio of
#'       cause-1 to cause-2 betas for shared predictors.}
#'     \item{\code{sparsity}}{Copied from \code{trspec$sparsity}.}
#'     \item{\code{target_eta_sd}}{Copied from \code{trspec$target_eta_sd}.}
#'     \item{\code{beta}}{Copied from \code{trspec$beta}.}
#'   }
#'
#' @seealso \code{\link{sim_spec}}, \code{\link{sim_generate}},
#'   \code{\link{sim_compute_effects}}
#'
#' @export
sim_generate_truth <- function(trspec, p, seed) {

  # Number of active predictors per cause
  if (trspec$sparsity < 1) {
    s <- ceiling(p * (trspec$sparsity %||% 0.05))
    s <- max(0L, min(p, as.integer(s)))
  } else if (trspec$sparsity >= 1 & trspec$sparsity <= p) {
    s <- trspec$sparsity
  }

  if (!is.null(seed)) {
    set.seed(seed)
  } else {
    set.seed(9134)
  }

  active_predictors       <- sample(1:p, 2 * s - trspec$share_active_set)
  active_predictors_cause1 <- active_predictors[1:s]
  active_predictors_cause2 <- active_predictors[((s - trspec$share_active_set) + 1):length(active_predictors)]
  shared_predictors        <- intersect(active_predictors_cause1,
                                        active_predictors_cause2)

  list(shared_multiplier  = trspec$shared_multiplier,
       cause1             = active_predictors_cause1,
       cause2             = active_predictors_cause2,
       shared_predictors  = shared_predictors,
       share_active_set   = trspec$share_active_set,
       sparsity           = trspec$sparsity,
       target_eta_sd      = trspec$target_eta_sd,
       beta               = trspec$beta)
}


#' Compute true coefficients and linear predictors
#'
#' Given the design matrix \eqn{X} and the resolved truth list from
#' \code{\link{sim_generate_truth}}, draws regression coefficients for each
#' cause and rescales them so that the resulting linear predictor
#' \eqn{\eta = X\beta} achieves the target standard deviation specified in
#' \code{truth$target_eta_sd}.  When causes share active predictors the
#' non-shared betas for cause 2 are further scaled via \code{optimize} to
#' hit the target \eqn{\eta} SD while preserving the specified
#' \code{shared_multiplier} ratio.
#'
#' @param X Numeric matrix \eqn{n \times p}. Design matrix (will be
#'   mean-centered internally; the original \code{X} is not modified outside
#'   this function).
#' @param truth Named list as returned by \code{\link{sim_generate_truth}}.
#'   Used fields: \code{cause1}, \code{cause2}, \code{share_active_set},
#'   \code{shared_predictors}, \code{shared_multiplier}, \code{beta},
#'   \code{target_eta_sd}.
#'
#' @return A named list with two sub-lists, \code{effects1} and
#'   \code{effects2}, each containing:
#'   \describe{
#'     \item{\code{betas1} / \code{betas2}}{Numeric vector of scaled
#'       coefficients for the active predictors of that cause.}
#'     \item{\code{betas1_sd} / \code{betas2_sd}}{Standard deviation of
#'       those coefficients.}
#'     \item{\code{eta1} / \code{eta2}}{Numeric vector (\eqn{n \times 1})
#'       of linear predictor values.}
#'     \item{\code{eta1_mean} / \code{eta2_mean}}{Mean of \eqn{\eta}
#'       (approximately zero by construction).}
#'     \item{\code{eta1_sd} / \code{eta2_sd}}{Standard deviation of
#'       \eqn{\eta} (equals \code{target_eta_sd} for cause 1, and
#'       approximately so for cause 2).}
#'   }
#'
#' @importFrom stats optimize
#'
#' @seealso \code{\link{sim_generate_truth}}, \code{\link{sim_generate}}
#'
#' @export
sim_compute_effects <- function(X, truth) {

  # Mean-center upfront so all eta have mean 0 by construction
  X <- scale(X, center = TRUE, scale = FALSE)

  # --- Cause 1 ---
  X1             <- X[, truth$cause1, drop = FALSE]
  target_beta_sd <- truth$beta$sd[1]
  target_eta_sd  <- truth$target_eta_sd[1]
  p1             <- ncol(X1)

  betas1      <- rnorm(p1, 0, target_beta_sd)
  eta1_0      <- as.numeric(X1 %*% betas1)
  print(mean(eta1_0))
  scale_factor <- target_eta_sd / sd(eta1_0)
  betas1       <- betas1 * scale_factor
  eta1         <- as.numeric(X1 %*% betas1)
  print(mean(eta1))

  # --- Cause 2 ---
  X2 <- X[, truth$cause2, drop = FALSE]

  if (length(truth$beta$sd) >= 2)       target_beta_sd <- truth$beta$sd[2]
  if (length(truth$target_eta_sd) >= 2) target_eta_sd  <- truth$target_eta_sd[2]

  p2     <- ncol(X2)
  betas2 <- numeric(p2)

  if (truth$share_active_set < 1) {
    betas2       <- rnorm(p2, 0, target_beta_sd)
    eta2_0       <- as.numeric(X2 %*% betas2)
    scale_factor <- target_eta_sd / sd(eta2_0)
    betas2       <- betas2 * scale_factor
  } else {
    i1 <- match(truth$shared_predictors, truth$cause1)
    i2 <- match(truth$shared_predictors, truth$cause2)

    if (length(truth$shared_multiplier) == 1)
      truth$shared_multiplier <- rep(truth$shared_multiplier, truth$share_active_set)

    betas2[i2]  <- betas1[i1] * (1 / truth$shared_multiplier)
    betas2[-i2] <- rnorm(p2 - length(i2), 0, target_beta_sd)

    eta_overlap        <- as.numeric(X2[,  i2, drop = FALSE] %*% betas2[i2])
    eta_nonoverlap_raw <- as.numeric(X2[, -i2, drop = FALSE] %*% betas2[-i2])

    # Solve for the scalar c that makes sd(eta_overlap + c * eta_nonoverlap_raw)
    # equal to target_eta_sd while preserving the shared-predictor betas exactly.
    objective <- function(c) {
      (sd(eta_overlap + c * eta_nonoverlap_raw) - target_eta_sd)^2
    }
    opt_res     <- optimize(objective, interval = c(-100, 100))
    betas2[-i2] <- opt_res$minimum * betas2[-i2]
  }

  eta2 <- as.numeric(X2 %*% betas2)

  list(
    effects1 = list(
      betas1     = betas1,
      betas1_sd  = sd(betas1),
      eta1       = eta1,
      eta1_mean  = mean(eta1),
      eta1_sd    = sd(eta1)
    ),
    effects2 = list(
      betas2     = betas2,
      betas2_sd  = sd(betas2),
      eta2       = eta2,
      eta2_mean  = mean(eta2),
      eta2_sd    = sd(eta2)
    )
  )
}
