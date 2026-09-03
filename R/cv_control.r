#' Build a CV control object for bhCRR cross-validation
#'
#' Bundles every setting that governs cross-validation and per-fold fitting
#' into a single validated object.  Pass this to \code{cv_ssl_psdh},
#' \code{tune_ssl_psdh}, and \code{bhcrr_autotune} instead of threading
#' settings through \code{...}.
#'
#' @param nfolds Integer \eqn{\ge 2}.  Number of CV folds per repetition.
#'   Default \code{10L}.
#' @param ncv Integer \eqn{\ge 1}.  Number of independent CV repetitions.
#'   Each repetition uses a fresh random fold assignment.  Default \code{2L}.
#' @param foldid Integer matrix of dimensions \eqn{n \times \code{ncv}} with
#'   pre-specified fold assignments (values in \eqn{1, \ldots, \code{nfolds}}).
#'   When supplied, \code{nfolds}, \code{ncv} and \code{strata} are ignored for
#'   fold generation; the column count of \code{foldid} becomes the effective
#'   \code{ncv}.  The row count (n) is verified at fit time when n is known.
#'   Default \code{NULL} (generate folds randomly).
#' @param strata Character.  How to stratify fold assignment.
#'   \code{"cause1"} (default) stratifies on the binary indicator
#'   \code{status == 1}; \code{"status"} stratifies on the full three-level
#'   status code; \code{"none"} performs a plain unstratified shuffle.
#'   Partially matched via \code{match.arg()}.
#' @param pool Numeric in \eqn{(0, 0.5)}.  Passed to \pkg{rsample}'s
#'   stratification as the minimum stratum fraction before pooling.  Default
#'   \code{0.01}.  Note: \pkg{rsample} defaults to \code{0.1}; if you set
#'   \code{pool >= 0.1} with a stratified \code{strata} value and your cause-1
#'   rate is below 10\%, folds will be silently unstratified — a warning is
#'   issued to flag this.
#' @param eval_quantile Numeric in \eqn{(0, 1)}.  Fallback evaluation horizon:
#'   when \code{eval_time} is \code{NULL}, \eqn{\tau} is derived as this
#'   quantile of the cause-1 event times in the \emph{full} data passed to the
#'   CV function.  This is a marginal, model-free quantity — it touches no
#'   fitted model — but it does consume all of \eqn{y}, so supplying
#'   \code{eval_time} explicitly is preferred for any result going into the
#'   manuscript.  Default \code{0.5}.
#' @param eval_time Numeric scalar or \code{NULL}.  A fixed, a-priori evaluation
#'   horizon \eqn{\tau} (in the same time units as \code{y[, 1]}) used for both
#'   the per-fold and the pooled Wolbers C-index.  When supplied, \code{eval_time}
#'   takes precedence over \code{eval_quantile}.  Validate: must be a single
#'   finite positive number or \code{NULL}.  Default \code{NULL}.
#' @param init_method Character string naming a built-in initialization
#'   (\code{"LASSO_cv"}, \code{"LASSO_bic"}, \code{"zero"}), a function with
#'   signature \code{function(x, y, ...)}, or a string naming such a function
#'   in the calling environment.  Validated eagerly at control construction so
#'   a typo is caught here rather than inside hundreds of fold fits.
#'   Default \code{"LASSO_cv"}.
#' @param init_args Named list.  Extra arguments forwarded to the
#'   \code{init_method} function via \code{do.call}.  Default \code{list()}.
#' @param warm_start Logical.  If \code{TRUE} (default), per-fold coefficients
#'   from one \code{(s0, s1)} pair are used as the starting point for the next
#'   pair in the grid traversal order, avoiding re-initialisation from scratch.
#' @param fit_args Named list.  Additional arguments forwarded to
#'   \code{\link{fit_ssl_psdh}} for every fold fit (e.g.\ \code{maxit},
#'   \code{epsilon}, \code{theta_a}, \code{theta_b}, \code{initial_sparsity},
#'   \code{inner_maxit_start}).  This is the only escape hatch: CV functions
#'   do not accept \code{...}.  Supplying \code{"init"}, \code{"init_method"},
#'   \code{"init_args"}, \code{"x"} or \code{"y"} here is an error; use the
#'   dedicated fields instead.
#' @param parallel Logical.  When \code{TRUE}, fold tasks are dispatched via
#'   \code{parallel::mclapply} (forking).  Width is \code{nfolds * ncv} under
#'   \code{warm_start = TRUE} (one task per fold chain) and
#'   \code{nfolds * ncv * npairs} under \code{warm_start = FALSE} (one task
#'   per pair, no chain state to preserve).  Set \code{control$workers} to cap
#'   the core count; on CHPC, set this from the SLURM allocation as
#'   \code{detectCores()} reports the whole node.  Forking is unreliable inside
#'   RStudio and Positron — a message is emitted when \code{parallel = TRUE}
#'   and \code{interactive()} is \code{TRUE}.  Default \code{FALSE}.
#' @param workers Integer or \code{NULL}.  Number of parallel workers.
#'   Ignored when \code{parallel = FALSE}.  Default \code{NULL} (auto).
#' @param seed Integer or \code{NULL}.  RNG seed set before fold generation.
#'   Default \code{NULL} (no seed).
#' @param keep_coefs Logical.  Whether \code{bhcrr_cv()} returns fitted
#'   coefficients for every (pair, fold) combination.  Default \code{FALSE}.
#'   At \eqn{p = 24618} a 24-pair grid at \code{nfolds = 10}, \code{ncv = 2}
#'   is roughly 95 MB of doubles, so the default suppresses that output.
#'
#' @returns A list of class \code{"bhcrr_cv_control"} containing all validated
#'   settings.  The field \code{$init_label} stores the resolved method label
#'   (as returned by \code{.resolve_init_method()}) alongside the original
#'   \code{$init_method}.
#'
#' @seealso \code{\link{cv_ssl_psdh}}, \code{\link{tune_ssl_psdh}},
#'   \code{\link{bhcrr_autotune}}, \code{\link{fit_ssl_psdh}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ctrl <- bhcrr_cv_control(nfolds = 5, ncv = 3, init_method = "LASSO_bic")
#' ctrl
#'
#' # Pass extra fit_ssl_psdh arguments via fit_args
#' ctrl2 <- bhcrr_cv_control(
#'   fit_args = list(maxit = 100, theta_b = 1)
#' )
#' }
bhcrr_cv_control <- function(
    nfolds        = 10L,
    ncv           = 2L,
    foldid        = NULL,
    strata        = c("cause1", "none", "status"),
    pool          = 0.01,
    eval_quantile = 0.5,
    init_method   = "LASSO_cv",
    init_args     = list(),
    warm_start    = TRUE,
    fit_args      = list(),
    parallel      = FALSE,
    workers       = NULL,
    seed          = NULL,
    eval_time     = NULL,
    keep_coefs    = FALSE
) {
  # ---- strata ----
  strata <- match.arg(strata)

  # ---- nfolds / ncv ----
  nfolds <- as.integer(nfolds)
  ncv    <- as.integer(ncv)
  if (is.na(nfolds) || nfolds < 2L)
    stop("'nfolds' must be an integer >= 2.")
  if (is.na(ncv) || ncv < 1L)
    stop("'ncv' must be an integer >= 1.")

  # ---- foldid ----
  if (!is.null(foldid)) {
    foldid <- as.matrix(foldid)
    if (!is.numeric(foldid) && !is.integer(foldid))
      stop("'foldid' must be a numeric or integer matrix.")
    if (any(!is.finite(foldid)) || any(foldid < 1L))
      stop("'foldid' must contain positive integers (fold indices >= 1).")
    # n check deferred to fold-generation time when n is known
  }

  # ---- pool ----
  if (!is.numeric(pool) || length(pool) != 1L || is.na(pool) ||
      pool <= 0 || pool >= 0.5)
    stop("'pool' must be a single number in (0, 0.5).")
  if (pool >= 0.1 && strata != "none")
    warning(
      "'pool' >= 0.1 with strata = '", strata, "': rsample pools any stratum ",
      "whose fraction is below 'pool' into the remaining observations. ",
      "If the cause-1 rate in your data is below 10%, fold assignment will be ",
      "silently unstratified. Consider reducing 'pool' (e.g. pool = 0.01) or ",
      "setting strata = 'none'.",
      call. = FALSE
    )

  # ---- eval_quantile ----
  if (!is.numeric(eval_quantile) || length(eval_quantile) != 1L ||
      is.na(eval_quantile) || eval_quantile <= 0 || eval_quantile >= 1)
    stop("'eval_quantile' must be a single number in (0, 1).")

  # ---- init_method: eager validation ----
  # Resolve now so a typo surfaces at control construction, not during fitting.
  caller_env  <- parent.frame()
  resolved    <- .resolve_init_method(init_method, envir = caller_env)
  init_label  <- resolved$label

  if (!is.list(init_args))
    stop("'init_args' must be a named list.")

  # ---- warm_start ----
  if (!is.logical(warm_start) || length(warm_start) != 1L || is.na(warm_start))
    stop("'warm_start' must be TRUE or FALSE.")

  # ---- fit_args ----
  if (!is.list(fit_args))
    stop("'fit_args' must be a named list.")

  reserved <- c("init", "init_method", "init_args", "x", "y")
  bad_reserved <- intersect(names(fit_args), reserved)
  if (length(bad_reserved) > 0L)
    stop(
      "fit_args must not contain: ", paste(bad_reserved, collapse = ", "), ". ",
      "Supply 'init' / 'init_method' / 'init_args' as dedicated fields of ",
      "bhcrr_cv_control(), and 'x' / 'y' directly to the CV function."
    )

  allowed_fit_formals <- setdiff(
    names(formals(fit_ssl_psdh)),
    c("x", "y", "init", "init_method", "init_args", "init_lam_path")
  )
  bad_names <- setdiff(names(fit_args), allowed_fit_formals)
  if (length(bad_names) > 0L)
    stop(
      "fit_args contains name(s) that are not formals of fit_ssl_psdh(): ",
      paste(bad_names, collapse = ", "), ". ",
      "Allowed names: ", paste(sort(allowed_fit_formals), collapse = ", "), "."
    )

  # ---- parallel / workers ----
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel))
    stop("'parallel' must be TRUE or FALSE.")
  if (!is.null(workers)) {
    workers <- as.integer(workers)
    if (is.na(workers) || workers < 1L)
      stop("'workers' must be a positive integer or NULL.")
  }

  # ---- seed ----
  if (!is.null(seed)) {
    seed <- as.integer(seed)
    if (is.na(seed))
      stop("'seed' must be an integer or NULL.")
  }

  # ---- eval_time ----
  if (!is.null(eval_time)) {
    if (!is.numeric(eval_time) || length(eval_time) != 1L ||
        !is.finite(eval_time) || eval_time <= 0)
      stop("'eval_time' must be a single finite positive number or NULL.")
  }

  # ---- keep_coefs ----
  if (!is.logical(keep_coefs) || length(keep_coefs) != 1L || is.na(keep_coefs))
    stop("'keep_coefs' must be TRUE or FALSE.")

  structure(
    list(
      nfolds        = nfolds,
      ncv           = ncv,
      foldid        = foldid,
      strata        = strata,
      pool          = pool,
      eval_quantile = eval_quantile,
      init_method   = init_method,
      init_label    = init_label,
      init_args     = init_args,
      warm_start    = warm_start,
      fit_args      = fit_args,
      parallel      = parallel,
      workers       = workers,
      seed          = seed,
      eval_time     = eval_time,
      keep_coefs    = keep_coefs
    ),
    class = "bhcrr_cv_control"
  )
}

#' Print a bhcrr_cv_control object
#'
#' @param x A \code{bhcrr_cv_control} object.
#' @param ... Ignored.
#'
#' @returns \code{x}, invisibly.
#' @export
print.bhcrr_cv_control <- function(x, ...) {
  fold_desc <- if (!is.null(x$foldid)) {
    sprintf("user-supplied (%d x %d)", nrow(x$foldid), ncol(x$foldid))
  } else {
    sprintf("%d folds x %d reps, strata = '%s', pool = %g",
            x$nfolds, x$ncv, x$strata, x$pool)
  }

  init_desc <- if (is.function(x$init_method)) {
    sprintf("custom function [label: %s]", x$init_label)
  } else {
    sprintf("'%s' [label: %s]", x$init_method, x$init_label)
  }

  fit_desc <- if (length(x$fit_args) == 0L) {
    "none"
  } else {
    paste(names(x$fit_args), collapse = ", ")
  }

  exec_parts <- sprintf("parallel = %s", x$parallel)
  if (x$parallel && !is.null(x$workers))
    exec_parts <- sprintf("%s, workers = %d", exec_parts, x$workers)
  if (!is.null(x$seed))
    exec_parts <- sprintf("%s, seed = %d", exec_parts, x$seed)

  horizon_desc <- if (!is.null(x$eval_time)) {
    sprintf("eval_time = %g (fixed)", x$eval_time)
  } else {
    sprintf("quantile %g (derived)", x$eval_quantile)
  }

  cat("<bhcrr_cv_control>\n")
  cat("  resampling : ", fold_desc, "\n", sep = "")
  cat("  horizon    : ", horizon_desc, "\n", sep = "")
  cat("  init       : ", init_desc,
      if (length(x$init_args)) sprintf(" + %d extra arg(s)", length(x$init_args)) else "",
      "\n", sep = "")
  cat("  warm_start : ", x$warm_start, "\n", sep = "")
  cat("  fit_args   : ", fit_desc, "\n", sep = "")
  cat("  execution  : ", exec_parts, "\n", sep = "")
  invisible(x)
}
