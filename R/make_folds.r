#' Generate stratified cross-validation folds for competing-risks data
#'
#' Creates an \eqn{n \times \code{ncv}} integer matrix of fold assignments,
#' optionally stratifying on the cause-1 indicator or full status code via
#' \pkg{rsample}.  Supersedes \code{\link{generate_foldid}}.
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
#' @seealso \code{\link{generate_foldid}} (older unstratified generator,
#'   retained for back-compatibility), \code{\link{bhcrr_cv_control}},
#'   \code{\link{cv_ssl_psdh}}
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
