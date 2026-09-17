# ---- CV control --------------------------------------------------------------

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
#' @seealso \code{\link{fit_ssl_psdh}}
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


# ---- fold generation ---------------------------------------------------------

#' Generate stratified cross-validation folds for competing-risks data
#'
#' Creates an \eqn{n \times \code{ncv}} integer matrix of fold assignments,
#' optionally stratifying on the cause-1 indicator or full status code via
#' \pkg{rsample}.  Supersedes \code{generate_foldid}.
#'
#' @param y Two-column numeric matrix of dimensions \eqn{n \times 2}.
#'   Column 1 is observed time; column 2 is status (0 censored, 1 cause 1,
#'   2 competing event).
#' @param control A \code{\link{bhcrr_cv_control}} object.  All folding
#'   parameters (\code{nfolds}, \code{ncv}, \code{foldid}, \code{strata},
#'   \code{pool}) are read from here.
#'
#' @returns A list with:
#'   \describe{
#'     \item{\code{foldid}}{Integer matrix \eqn{n \times \code{ncv}}. Entry
#'       \code{[i, k]} is the fold index of observation \eqn{i} in repetition
#'       \eqn{k}, in \eqn{1, \ldots, \code{nfolds}}.}
#'     \item{\code{nfolds}}{Integer. Effective number of folds.}
#'     \item{\code{ncv}}{Integer. Effective number of repetitions.}
#'     \item{\code{strata}}{Character. The stratification method used.}
#'     \item{\code{fold_event_counts}}{Integer matrix \eqn{\code{ncv} \times
#'       \code{nfolds}} giving the number of cause-1 events in each fold of
#'       each repetition.  Use this to spot degenerate splits before fitting.}
#'   }
#'
#' @seealso \code{\link{bhcrr_cv_control}}
#'
#' @importFrom rsample vfold_cv assessment
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ctrl  <- bhcrr_cv_control(nfolds = 5, ncv = 2, strata = "cause1")
#' folds <- bhcrr_make_folds(y, ctrl)
#' dim(folds$foldid)   # n x 2
#' folds$fold_event_counts
#' }
bhcrr_make_folds <- function(y, control) {
  if (!inherits(control, "bhcrr_cv_control"))
    stop("'control' must be a bhcrr_cv_control object.")
  if (!is.matrix(y) || ncol(y) != 2L)
    stop("'y' must be an n x 2 matrix (time, status).")

  n <- nrow(y)

  # ---- passthrough: user supplied foldid ----
  if (!is.null(control$foldid)) {
    fid <- control$foldid
    if (nrow(fid) != n)
      stop("nrow(control$foldid) = ", nrow(fid), " but nrow(y) = ", n, ".")
    effective_nfolds <- max(fid, na.rm = TRUE)
    effective_ncv    <- ncol(fid)
    foldid_mat   <- fid
    was_supplied <- TRUE
  } else {
    was_supplied <- FALSE

    # ---- cap and LOO adjustment ----
    effective_nfolds <- min(as.integer(control$nfolds), n)
    effective_ncv    <- if (effective_nfolds == n) 1L else as.integer(control$ncv)
    is_loo           <- (effective_nfolds == n)

    # ---- build minimal data frame for rsample (no x, no large columns) ----
    strat_col <- switch(
      control$strata,
      cause1 = factor(as.integer(y[, 2L] == 1L)),
      status = factor(y[, 2L]),
      none   = NULL
    )

    df <- if (is.null(strat_col)) {
      data.frame(.row = seq_len(n))
    } else {
      data.frame(.row = seq_len(n), .strat = strat_col)
    }

    # ---- LOO: rsample does not support vfold_cv with v == n ----
    if (is_loo) {
      foldid_mat <- matrix(seq_len(n), nrow = n, ncol = 1L)
    } else {
      # ---- run rsample, suppressing repeated low-stratum-size warnings ----
      rset <- .dedupe_warnings(
        if (is.null(strat_col)) {
          rsample::vfold_cv(df,
                            v       = effective_nfolds,
                            repeats = effective_ncv)
        } else {
          rsample::vfold_cv(df,
                            v       = effective_nfolds,
                            repeats = effective_ncv,
                            strata  = ".strat",
                            pool    = control$pool)
        }
      )

      # ---- convert rset to n x ncv integer matrix ----
      # With repeats > 1: id = Repeat*, id2 = Fold*
      # With repeats == 1: id = Fold*, no id2 column
      has_id2 <- "id2" %in% names(rset)

      if (has_id2) {
        rep_labels  <- unique(rset$id)    # "Repeat1", "Repeat2", ...
        fold_labels <- unique(rset$id2)   # "Fold1",   "Fold2", ...
      } else {
        rep_labels  <- "Repeat1"
        fold_labels <- unique(rset$id)    # "Fold1", ..., "FoldV"
      }

      foldid_mat <- matrix(NA_integer_, nrow = n, ncol = length(rep_labels))

      for (k in seq_along(rep_labels)) {
        sub <- if (has_id2) rset[rset$id == rep_labels[k], ] else rset

        for (j in seq_along(fold_labels)) {
          fold_lbl <- fold_labels[j]
          row_sub  <- if (has_id2) sub[sub$id2 == fold_lbl, ] else sub[sub$id == fold_lbl, ]
          assess_rows <- rsample::assessment(row_sub$splits[[1L]])$.row
          foldid_mat[assess_rows, k] <- j
        }
      }
    }
  }

  # ---- compute per-fold cause-1 event counts ----
  cause1        <- as.integer(y[, 2L] == 1L)
  ec_nfolds     <- max(foldid_mat, na.rm = TRUE)
  ec_ncv        <- ncol(foldid_mat)
  fold_event_counts <- matrix(NA_integer_, nrow = ec_ncv, ncol = ec_nfolds)

  for (k in seq_len(ec_ncv)) {
    for (j in seq_len(ec_nfolds)) {
      fold_event_counts[k, j] <- sum(cause1[foldid_mat[, k] == j])
    }
  }

  # ---- diagnostics ----
  if (any(fold_event_counts == 0L)) {
    zero_idx <- which(fold_event_counts == 0L, arr.ind = TRUE)
    warning(
      "One or more folds have zero cause-1 events: ",
      paste(apply(zero_idx, 1L, function(r)
        sprintf("rep %d fold %d", r[1L], r[2L])), collapse = "; "),
      ". wolbers_c() will return NA for those folds.",
      call. = FALSE
    )
  }

  if (!was_supplied && control$strata != "none") {
    spread <- max(fold_event_counts) - min(fold_event_counts)
    if (spread > 1L) {
      warning(
        "Stratification by '", control$strata, "' with pool = ", control$pool,
        " produced uneven cause-1 fold counts (spread = ", spread, "): ",
        paste(as.vector(fold_event_counts), collapse = ", "),
        ". True stratification cannot differ by more than 1. ",
        "The pool value may be too high, causing rsample to fall back to ",
        "an unstratified shuffle.",
        call. = FALSE
      )
    }
  }

  list(
    foldid            = foldid_mat,
    nfolds            = max(foldid_mat, na.rm = TRUE),
    ncv               = ncol(foldid_mat),
    strata            = control$strata,
    fold_event_counts = fold_event_counts
  )
}


# ---- hyperparameter grid -----------------------------------------------------

#' Build and order the (s0, s1) hyperparameter grid
#'
#' Internal helper used by \code{bhcrr_cv()}.  Computes the Cartesian product
#' of \code{s0_seq} and \code{s1_seq}, drops pairs where \code{s1 <= s0}, and
#' orders the survivors in warm-start traversal order: unique \code{s1} values
#' in the order they appear in \code{s1_seq}, and within each group \code{s0}
#' ascending.  This matches the traversal in \code{tune_ssl_psdh}.
#'
#' @param s0_seq Numeric vector of candidate spike scale values.
#' @param s1_seq Numeric vector of candidate slab scale values.
#'
#' @return A \code{data.frame} with columns \code{s0}, \code{s1}, and
#'   \code{pair} (integer traversal index, 1..npairs).
#'
#' @keywords internal
.cv_grid <- function(s0_seq, s1_seq) {
  grid  <- expand.grid(s0 = s0_seq, s1 = s1_seq, stringsAsFactors = FALSE)
  valid <- grid[grid$s1 > grid$s0, , drop = FALSE]

  if (nrow(valid) == 0L)
    stop(
      "No valid (s0, s1) pairs after dropping s1 <= s0. ",
      "s0 range: [", min(s0_seq), ", ", max(s0_seq), "]; ",
      "s1 range: [", min(s1_seq), ", ", max(s1_seq), "]."
    )

  # Traversal order: unique s1 in s1_seq order, s0 ascending within each group.
  # expand.grid() places s1 in s1_seq order (s0 varies fastest), so
  # unique(valid$s1) preserves that order — matching tune_ssl_psdh exactly.
  s1_levels <- unique(valid$s1)

  rows <- lapply(s1_levels, function(s1_val) {
    s0_vals <- sort(valid$s0[valid$s1 == s1_val])
    data.frame(s0 = s0_vals, s1 = s1_val, stringsAsFactors = FALSE)
  })

  result       <- do.call(rbind, rows)
  result$pair  <- seq_len(nrow(result))
  rownames(result) <- NULL
  result
}


# ---- per-fold workers --------------------------------------------------------

#' Compute the initialization vector for one (repetition, fold)
#'
#' Extracts the training subset, resolves the init method from \code{control},
#' and returns the validated initial coefficient vector.  Factored out of
#' \code{.cv_fold_path()} so it can be dispatched independently (e.g. in
#' Phase 1 of the wide parallel path where warm-start state cannot cross
#' worker boundaries).
#'
#' @param x Full design matrix, \eqn{n \times p}.
#' @param y Full outcome matrix, \eqn{n \times 2} (time, status).
#' @param train_idx Integer vector.  Row indices of the training fold.
#' @param control A \code{\link{bhcrr_cv_control}} object.
#'
#' @return Numeric vector of length \eqn{p}: the validated initial coefficients.
#'
#' @keywords internal
.cv_fold_init <- function(x, y, train_idx, control) {
  x_train    <- x[train_idx, , drop = FALSE]
  y_train    <- y[train_idx, , drop = FALSE]
  caller_env <- parent.frame()
  resolved   <- .resolve_init_method(control$init_method, envir = caller_env)
  raw        <- do.call(resolved$fn,
                        c(list(x = x_train, y = y_train), control$init_args))
  .validate_init(raw, ncol(x), resolved$label)$init
}


#' Per-fold grid worker for bhcrr_cv()
#'
#' Runs one (repetition, fold) across the entire hyperparameter grid in
#' traversal order, optionally carrying a warm-start coefficient chain from
#' pair to pair.  This is the parallelisable unit of \code{bhcrr_cv()}: it is
#' self-contained, has no side effects, and never calls \code{message()} or
#' \code{warning()}.
#'
#' @param x Full design matrix, \eqn{n \times p}.
#' @param y Full outcome matrix, \eqn{n \times 2} (time, status).
#' @param train_idx Integer vector.  Row indices of the training fold.
#' @param test_idx Integer vector.  Row indices of the test fold.
#' @param grid Data frame from \code{\link{.cv_grid}} with columns
#'   \code{s0}, \code{s1}, \code{pair}.
#' @param control A \code{\link{bhcrr_cv_control}} object.
#' @param eval_time Numeric scalar.  Already-resolved evaluation horizon
#'   (the worker never derives it from \code{control}).
#' @param init Numeric vector of length \eqn{p}, or \code{NULL}.  When
#'   supplied, used directly as the starting coefficients and the
#'   \code{init_method} in \code{control} is never called.  When \code{NULL}
#'   the init method is invoked exactly once on the training fold.
#'
#' @return A list:
#'   \describe{
#'     \item{\code{lp}}{Numeric matrix, \code{length(test_idx)} rows by
#'       \code{nrow(grid)} columns.  Entry \code{[i, j]} is the predicted
#'       absolute risk for test observation \code{i} at pair \code{j},
#'       computed by \code{\link{predict.ssl_psdh}}.  \code{NA} for
#'       pairs where the fit failed.}
#'     \item{\code{init}}{Numeric vector of length \eqn{p}.  The initial
#'       coefficient vector actually used.}
#'     \item{\code{iterations}}{Integer vector, one entry per grid pair.
#'       \code{NA} where the fit errored.  When \code{conv[j]} is
#'       \code{FALSE} and \code{iterations[j] == maxit} (from
#'       \code{control\$fit_args}), the outer EM exhausted its iteration
#'       budget; when \code{conv[j]} is \code{FALSE} and
#'       \code{iterations[j] < maxit}, the inner \pkg{fastcmprsk} solver
#'       hit its escalation ceiling and broke early — the only way to
#'       distinguish the two non-convergence modes.}
#'     \item{\code{conv}}{Logical vector, one entry per grid pair.
#'       \code{TRUE} if the EM converged within \code{maxit} iterations,
#'       \code{FALSE} if it did not, \code{NA} where the fit errored.}
#'     \item{\code{errors}}{Data frame with columns \code{pair}, \code{s0},
#'       \code{s1}, \code{message}; zero rows when all fits succeeded.}
#'     \item{\code{n_failed}}{Integer.  Number of failed fits.}
#'     \item{\code{coefs}}{List of length-\eqn{p} numeric vectors (one per
#'       pair) when \code{control$keep_coefs} is \code{TRUE}; otherwise
#'       \code{NULL}.}
#'   }
#'
#' @keywords internal
.cv_fold_path <- function(x, y, train_idx, test_idx, grid, control, eval_time,
                          init = NULL) {
  x_train <- x[train_idx, , drop = FALSE]
  y_train <- y[train_idx, , drop = FALSE]
  x_test  <- x[test_idx,  , drop = FALSE]

  n_test  <- length(test_idx)
  n_pairs <- nrow(grid)
  p       <- ncol(x)

  # ---- initialisation: computed at most once, never per pair ----
  if (!is.null(init)) {
    b_init <- .validate_init(init, p, "cached")$init
  } else {
    b_init <- .cv_fold_init(x, y, train_idx, control)
  }

  # ---- output containers ----
  lp_mat      <- matrix(NA_real_, nrow = n_test, ncol = n_pairs)
  iters       <- rep(NA_integer_,  n_pairs)
  conv        <- rep(NA,           n_pairs)   # logical; NA where the fit errored
  errors_list <- vector("list",    n_pairs)
  coef_list   <- if (isTRUE(control$keep_coefs)) vector("list", n_pairs) else NULL

  # b   — coefficients passed to the next fit
  # b_c — last SUCCESSFUL coefficients (chain repair fallback)
  b   <- b_init
  b_c <- b_init

  for (j in seq_len(n_pairs)) {
    ss_j <- c(grid$s0[j], grid$s1[j])

    fit_result <- tryCatch(
      do.call(fit_ssl_psdh,
              c(list(x = x_train, y = y_train, ss = ss_j, init = b),
                control$fit_args)),
      error = function(e) e
    )

    if (inherits(fit_result, "error")) {
      errors_list[[j]] <- data.frame(
        pair    = grid$pair[j],
        s0      = grid$s0[j],
        s1      = grid$s1[j],
        message = conditionMessage(fit_result),
        stringsAsFactors = FALSE
      )
      # Chain repair: revert to the last good state without re-invoking init.
      b <- b_c
    } else {
      iters[j]   <- as.integer(fit_result$iterations)
      conv[j]    <- isTRUE(fit_result$conv)
      lp_mat[, j] <- predict(fit_result,
                             newx            = x_test,
                             prediction_time = eval_time)
      if (!is.null(coef_list))
        coef_list[[j]] <- as.numeric(fit_result$final_model_object$coef)

      new_coef <- as.numeric(fit_result$final_model_object$coef)
      if (isTRUE(control$warm_start)) {
        b   <- new_coef
        b_c <- new_coef
      } else {
        b   <- b_init
        b_c <- b_init
      }
    }
  }

  non_null <- Filter(Negate(is.null), errors_list)
  errors_df <- if (length(non_null) > 0L) {
    do.call(rbind, non_null)
  } else {
    data.frame(pair    = integer(0L),
               s0      = numeric(0L),
               s1      = numeric(0L),
               message = character(0L),
               stringsAsFactors = FALSE)
  }

  list(
    lp         = lp_mat,
    init       = b_init,
    iterations = iters,
    conv       = conv,
    errors     = errors_df,
    n_failed   = nrow(errors_df),
    coefs      = coef_list
  )
}
