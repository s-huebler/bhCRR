#' Cross-validate bhCRR over a hyperparameter grid
#'
#' Runs K-fold cross-validation (with optional repetitions) over the
#' \code{(s0, s1)} grid produced by \code{\link{.cv_grid}}.  Each
#' (repetition, fold) is handled by \code{\link{.cv_fold_path}}, which fits
#' the full grid in warm-start traversal order.  This function is sequential;
#' parallel execution is reserved for a future commit.
#'
#' @param x Numeric matrix (\eqn{n \times p}).
#' @param y Numeric matrix (\eqn{n \times 2}).  Column 1 is event time;
#'   column 2 is status in \eqn{\{0, 1, 2\}}.
#' @param s0_seq Numeric vector of candidate spike scale values.
#' @param s1_seq Numeric vector of candidate slab scale values.
#' @param control A \code{\link{bhcrr_cv_control}} object.
#' @param folds Output of \code{\link{bhcrr_make_folds}}, or \code{NULL} to
#'   generate folds from \code{y} and \code{control}.
#' @param fold_inits Nested list \code{[[k]][[i]]} of numeric vectors (one per
#'   (rep, fold) combination), or \code{NULL}.  When supplied the
#'   \code{init_method} is never called; pass \code{\$fold_inits} from a
#'   previous \code{bhcrr_cv()} call for exact coefficient-chain
#'   reproducibility without re-invoking initialisation.
#'
#' @return An object of class \code{"bhcrr_cv"} containing:
#'   \describe{
#'     \item{\code{tuning}}{Data frame with columns \code{s0}, \code{s1},
#'       \code{score_mean}, \code{score_sd}, \code{score_fold_sd},
#'       \code{n_failed_folds}; one row per grid pair.}
#'     \item{\code{best}}{The row of \code{tuning} with the highest
#'       \code{score_mean}, or \code{NULL} if every fit failed.}
#'     \item{\code{pooled}}{Numeric matrix (\code{n_pairs} x \code{ncv}).
#'       Entry \code{[j, k]} is the Wolbers C-index for pair \code{j} in
#'       repetition \code{k}, pooling all out-of-sample predictions.}
#'     \item{\code{fold_scores}}{Numeric array (\code{n_pairs} x
#'       \code{nfolds} x \code{ncv}).  Per-fold Wolbers C-index.}
#'     \item{\code{fold_inits}}{List \code{[[k]][[i]]} of the init vector
#'       actually used for each (rep, fold).  \code{NULL} entries mark
#'       folds that failed before initialisation.}
#'     \item{\code{eval_time}}{List with \code{\$value} (the resolved
#'       \eqn{\tau}) and \code{\$source} (\code{"supplied"} or
#'       \code{"quantile"}).}
#'     \item{\code{errors}}{Data frame with columns \code{rep}, \code{fold},
#'       \code{pair}, \code{s0}, \code{s1}, \code{message}.  Zero rows when
#'       clean; \code{pair}/\code{s0}/\code{s1} are \code{NA} for fold-level
#'       failures.}
#'     \item{\code{folds}}{The folds object used (output of
#'       \code{\link{bhcrr_make_folds}}).}
#'     \item{\code{control}}{The \code{bhcrr_cv_control} as used.}
#'     \item{\code{grid}}{The traversal-ordered grid from
#'       \code{\link{.cv_grid}}.}
#'     \item{\code{timing}}{List with \code{\$elapsed}, \code{\$workers},
#'       \code{\$parallel}.}
#'     \item{\code{coefs}}{Present only when \code{control\$keep_coefs}:
#'       list \code{[[k]][[i]]} of coef lists from each fold's
#'       \code{.cv_fold_path()} call.}
#'   }
#'
#' @section Failure policy:
#' If any fits fail a single \code{warning()} is emitted at the end
#' summarising the total count, distinct pairs affected, and the first error
#' message.  If \emph{every} fit fails \code{stop()} is called; an all-NA
#' tuning frame is never returned.
#'
#' @seealso \code{\link{bhcrr_cv_control}}, \code{\link{bhcrr_make_folds}},
#'   \code{\link{.cv_fold_path}}, \code{\link{.cv_grid}},
#'   \code{\link{wolbers_c}}
#'
#' @export
bhcrr_cv <- function(x, y, s0_seq, s1_seq, control = bhcrr_cv_control(),
                     folds = NULL, fold_inits = NULL) {

  t_start <- proc.time()

  # ---- 1. Validate inputs ----
  if (!is.matrix(x) || !is.numeric(x))
    stop("'x' must be a numeric matrix.")
  n <- nrow(x)
  p <- ncol(x)

  if (!is.matrix(y) || ncol(y) != 2L)
    stop("'y' must be a matrix with 2 columns (time, status).")
  if (nrow(y) != n)
    stop("nrow(x) and nrow(y) must be equal.")
  if (!all(y[, 2] %in% c(0L, 1L, 2L)))
    stop("y[, 2] (status) must contain only 0, 1, or 2.")

  if (!inherits(control, "bhcrr_cv_control"))
    stop("'control' must be a bhcrr_cv_control object.")

  grid    <- .cv_grid(s0_seq, s1_seq)
  n_pairs <- nrow(grid)

  # ---- 2. Folds ----
  if (!is.null(folds)) {
    if (!is.list(folds) || !is.matrix(folds$foldid))
      stop("'folds' must be the output of bhcrr_make_folds().")
  } else {
    folds <- bhcrr_make_folds(y, control)
  }
  foldid <- folds$foldid
  nfolds <- folds$nfolds
  ncv    <- ncol(foldid)

  # ---- 3. Resolve horizon once ----
  if (!is.null(control$eval_time)) {
    tau        <- control$eval_time
    tau_source <- "supplied"
  } else {
    cause1_t <- y[y[, 2] == 1, 1]
    if (length(cause1_t) == 0L)
      stop("No cause-1 events in y; cannot derive eval_time from quantile.")
    tau        <- quantile(cause1_t, control$eval_quantile, names = FALSE)
    tau_source <- "quantile"
  }

  # ---- 4. Sequential loop over (rep, fold) ----
  fold_results   <- vector("list", ncv)
  test_idx_store <- vector("list", ncv)
  fold_init_out  <- vector("list", ncv)

  for (k in seq_len(ncv)) {
    fold_results[[k]]   <- vector("list", nfolds)
    test_idx_store[[k]] <- vector("list", nfolds)
    fold_init_out[[k]]  <- vector("list", nfolds)

    for (i in seq_len(nfolds)) {
      train_idx <- which(foldid[, k] != i)
      test_idx  <- which(foldid[, k] == i)
      test_idx_store[[k]][[i]] <- test_idx

      supplied_init <- if (!is.null(fold_inits)) fold_inits[[k]][[i]] else NULL

      result <- tryCatch(
        .cv_fold_path(x, y, train_idx, test_idx, grid, control, tau,
                      init = supplied_init),
        error = function(e) e
      )
      fold_results[[k]][[i]]  <- result
      fold_init_out[[k]][[i]] <- if (inherits(result, "error")) NULL
                                  else result$init
    }
  }

  # ---- 5. Assemble scores and errors ----
  pooled_mat     <- matrix(NA_real_, nrow = n_pairs, ncol = ncv)
  fold_score_arr <- array(NA_real_,  dim  = c(n_pairs, nfolds, ncv))
  n_failed_vec   <- integer(n_pairs)
  all_errors     <- list()

  coefs_out <- if (isTRUE(control$keep_coefs)) {
    lapply(seq_len(ncv), function(k) vector("list", nfolds))
  } else NULL

  for (k in seq_len(ncv)) {
    lp_rep <- matrix(NA_real_, nrow = n, ncol = n_pairs)

    for (i in seq_len(nfolds)) {
      result   <- fold_results[[k]][[i]]
      test_idx <- test_idx_store[[k]][[i]]

      if (inherits(result, "error")) {
        # Fold-level failure (e.g. init method threw): record one row; all
        # pairs in this fold count as failed.
        all_errors <- c(all_errors, list(data.frame(
          rep     = k,
          fold    = i,
          pair    = NA_integer_,
          s0      = NA_real_,
          s1      = NA_real_,
          message = conditionMessage(result),
          stringsAsFactors = FALSE
        )))
        n_failed_vec <- n_failed_vec + 1L

      } else {
        # Scatter lp into the rep-level matrix for pooled scoring.
        lp_rep[test_idx, ] <- result$lp

        # Pair-level errors from within the fold worker.
        if (nrow(result$errors) > 0L) {
          err <- result$errors
          all_errors <- c(all_errors, list(data.frame(
            rep     = k,
            fold    = i,
            pair    = err$pair,
            s0      = err$s0,
            s1      = err$s1,
            message = err$message,
            stringsAsFactors = FALSE
          )))
          n_failed_vec[err$pair] <- n_failed_vec[err$pair] + 1L
        }

        # Per-fold Wolbers C for each pair.
        for (j in seq_len(n_pairs)) {
          lp_j <- result$lp[, j]
          if (any(is.finite(lp_j))) {
            fold_score_arr[j, i, k] <- tryCatch(
              wolbers_c(y[test_idx, , drop = FALSE], lp_j, tau),
              error = function(e) NA_real_
            )
          }
        }

        if (!is.null(coefs_out))
          coefs_out[[k]][[i]] <- result$coefs
      }
    }

    # Pooled Wolbers C for rep k: scatter all folds' lp into lp_rep first,
    # then score.  NAs at failed-fold positions cause wolbers_c to return NA.
    for (j in seq_len(n_pairs)) {
      lp_j <- lp_rep[, j]
      if (any(is.finite(lp_j))) {
        pooled_mat[j, k] <- tryCatch(
          wolbers_c(y, lp_j, tau),
          error = function(e) NA_real_
        )
      }
    }
  }

  # ---- Build tuning data frame ----
  score_mean <- rowMeans(pooled_mat, na.rm = TRUE)
  score_mean[!is.finite(score_mean)] <- NA_real_

  score_sd <- apply(pooled_mat, 1L, function(v) sd(v, na.rm = TRUE))
  score_sd[!is.finite(score_sd)] <- NA_real_

  # fold_sds[j, k] = within-rep sd of per-fold scores for pair j, rep k
  fold_sds      <- apply(fold_score_arr, c(1L, 3L), function(v) sd(v, na.rm = TRUE))
  score_fold_sd <- rowMeans(fold_sds, na.rm = TRUE)
  score_fold_sd[!is.finite(score_fold_sd)] <- NA_real_

  tuning <- data.frame(
    s0             = grid$s0,
    s1             = grid$s1,
    score_mean     = score_mean,
    score_sd       = score_sd,
    score_fold_sd  = score_fold_sd,
    n_failed_folds = n_failed_vec,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  finite_idx <- which(is.finite(tuning$score_mean))
  best <- if (length(finite_idx) > 0L) {
    row <- tuning[which.max(tuning$score_mean), , drop = FALSE]
    rownames(row) <- NULL
    row
  } else NULL

  errors_df <- if (length(all_errors) > 0L) {
    do.call(rbind, all_errors)
  } else {
    data.frame(
      rep = integer(0), fold = integer(0), pair = integer(0),
      s0  = numeric(0), s1  = numeric(0), message = character(0),
      stringsAsFactors = FALSE
    )
  }

  # ---- 7. Failure reporting ----
  total_failed <- sum(n_failed_vec)
  if (total_failed > 0L) {
    first_msg <- errors_df$message[1L]
    if (is.null(best)) {
      stop(sprintf(
        paste0("Every fit failed; no tuning scores could be computed. ",
               "First error: %s"),
        first_msg
      ), call. = FALSE)
    }
    # Count pairs with at least one failure; fold-level failures affect all pairs.
    n_affected <- if (any(is.na(errors_df$pair))) {
      n_pairs
    } else {
      length(unique(errors_df$pair))
    }
    warning(
      sprintf("%d fit failure(s) across %d pair(s). First message: %s",
              total_failed, n_affected, first_msg),
      call. = FALSE
    )
  }

  # ---- 6. Return ----
  elapsed <- proc.time() - t_start

  out <- list(
    tuning      = tuning,
    best        = best,
    pooled      = pooled_mat,
    fold_scores = fold_score_arr,
    fold_inits  = fold_init_out,
    eval_time   = list(value = tau, source = tau_source),
    errors      = errors_df,
    folds       = folds,
    control     = control,
    grid        = grid,
    timing      = list(elapsed = elapsed, workers = 1L, parallel = FALSE)
  )
  if (!is.null(coefs_out)) out$coefs <- coefs_out

  structure(out, class = "bhcrr_cv")
}


#' Print a bhcrr_cv object
#'
#' @param x A \code{bhcrr_cv} object.
#' @param ... Ignored.
#'
#' @return \code{x}, invisibly.
#' @export
print.bhcrr_cv <- function(x, ...) {
  g  <- x$grid
  fl <- x$folds

  cat("<bhcrr_cv>\n")
  cat(sprintf("  grid     : %d pair(s) (s0: [%.4g, %.4g], s1: [%.4g, %.4g])\n",
              nrow(g), min(g$s0), max(g$s0), min(g$s1), max(g$s1)))

  fold_desc <- sprintf("%d folds \u00d7 %d reps, strata = '%s'",
                       fl$nfolds, fl$ncv, fl$strata)
  cat("  folds    :", fold_desc, "\n")

  tau_desc <- if (x$eval_time$source == "supplied") {
    sprintf("%.4g (supplied)", x$eval_time$value)
  } else {
    sprintf("%.4g (quantile %.2f)", x$eval_time$value, x$control$eval_quantile)
  }
  cat("  horizon  : \u03c4 =", tau_desc, "\n")

  if (!is.null(x$best)) {
    b <- x$best
    cat(sprintf("  best     : s0 = %.4g, s1 = %.4g, score = %.4f\n",
                b$s0, b$s1, b$score_mean))
  } else {
    cat("  best     : (none \u2014 all fits failed)\n")
  }

  cat(sprintf("  failures : %d fit(s)\n", sum(x$tuning$n_failed_folds)))

  elapsed_s <- x$timing$elapsed["elapsed"]
  cat(sprintf("  elapsed  : %.2f s\n", elapsed_s))

  invisible(x)
}
