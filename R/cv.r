#' Cross-validate bhCRR over a hyperparameter grid
#'
#' Runs K-fold cross-validation (with optional repetitions) over the
#' \code{(s0, s1)} grid produced by \code{\link{.cv_grid}}.  Each
#' (repetition, fold) is handled by \code{\link{.cv_fold_path}}, which fits
#' the full grid in warm-start traversal order.
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
#' @section Parallelism:
#' When \code{control$parallel} is \code{TRUE}, fold tasks are dispatched via
#' \code{parallel::mclapply} (forking).  Under \code{warm_start = TRUE} the
#' width is \code{nfolds x ncv} (one chain per fold).  Under
#' \code{warm_start = FALSE} chains carry no state so the width expands to
#' \code{nfolds x ncv x npairs} — but cold starts may converge less often
#' within \code{maxit}, so compare \code{$tuning$n_not_converged} before
#' choosing.
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

  # ---- 4. Dispatch (rep, fold) tasks ----
  tasks <- expand.grid(i = seq_len(nfolds), k = seq_len(ncv),
                       KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  n_tasks <- nrow(tasks)

  if (isTRUE(control$parallel)) {
    n_workers <- if (!is.null(control$workers)) {
      as.integer(control$workers)
    } else {
      min(n_tasks, max(1L, parallel::detectCores() - 1L))
    }
    if (interactive()) {
      message(sprintf(
        "bhcrr_cv: forking %d worker(s) for %d (rep, fold) task(s). ",
        n_workers, n_tasks),
        "Forking can hang in RStudio/Positron; use control$parallel = FALSE if the session freezes.",
        appendLF = TRUE
      )
    }
  } else {
    n_workers <- 1L
  }

  # Reproducibility: switch to L'Ecuyer-CMRG before fork so per-worker
  # streams are independent and, when seed is set, reproducible.
  old_rngkind <- RNGkind()[1L]
  on.exit(RNGkind(old_rngkind), add = TRUE)
  RNGkind("L'Ecuyer-CMRG")
  if (!is.null(control$seed)) set.seed(control$seed)

  run_fn <- if (isTRUE(control$parallel)) {
    function(X, FUN) parallel::mclapply(X, FUN,
                                        mc.cores       = n_workers,
                                        mc.preschedule = FALSE)
  } else {
    lapply
  }

  if (isTRUE(control$warm_start)) {
    # --- narrow path: one task per (rep, fold); entire grid per task ---
    dispatch_label <- "fold"

    res_flat <- run_fn(seq_len(n_tasks), function(t) {
      k_t <- tasks$k[t]; i_t <- tasks$i[t]
      train_idx <- which(foldid[, k_t] != i_t)
      test_idx  <- which(foldid[, k_t] == i_t)
      supplied  <- if (!is.null(fold_inits)) fold_inits[[k_t]][[i_t]] else NULL
      list(
        k        = k_t,
        i        = i_t,
        test_idx = test_idx,
        result   = tryCatch(
          .cv_fold_path(x, y, train_idx, test_idx, grid, control, tau,
                        init = supplied),
          error = function(e) e
        )
      )
    })

    # --- worker-death check ---
    if (length(res_flat) != n_tasks)
      stop(sprintf("mclapply returned %d results for %d (rep,fold) tasks; a worker process died.",
                   length(res_flat), n_tasks))
    for (t in seq_len(n_tasks)) {
      el <- res_flat[[t]]
      if (inherits(el, "try-error"))
        stop(sprintf("Worker died for rep %d fold %d: %s",
                     tasks$k[t], tasks$i[t], as.character(el)))
    }

    # --- unflatten into fold_results / test_idx_store / fold_init_out ---
    fold_results   <- vector("list", ncv)
    test_idx_store <- vector("list", ncv)
    fold_init_out  <- vector("list", ncv)
    for (k in seq_len(ncv)) {
      fold_results[[k]]   <- vector("list", nfolds)
      test_idx_store[[k]] <- vector("list", nfolds)
      fold_init_out[[k]]  <- vector("list", nfolds)
    }
    for (t in seq_len(n_tasks)) {
      el <- res_flat[[t]]; k_t <- el$k; i_t <- el$i
      test_idx_store[[k_t]][[i_t]] <- el$test_idx
      fold_results[[k_t]][[i_t]]   <- el$result
      fold_init_out[[k_t]][[i_t]]  <- if (inherits(el$result, "error")) NULL
                                       else el$result$init
    }

  } else {
    # --- wide path: warm_start = FALSE, parallel across (rep, fold, pair) ---
    dispatch_label <- "fold-pair"

    # Phase 1: compute one init per (k, i), in parallel.
    # Skip entirely if fold_inits was supplied.
    if (!is.null(fold_inits)) {
      init_mat <- fold_inits  # already [[k]][[i]]
    } else {
      init_flat <- run_fn(seq_len(n_tasks), function(t) {
        k_t <- tasks$k[t]; i_t <- tasks$i[t]
        train_idx <- which(foldid[, k_t] != i_t)
        tryCatch(
          .cv_fold_init(x, y, train_idx, control),
          error = function(e) e
        )
      })
      if (length(init_flat) != n_tasks)
        stop(sprintf("mclapply returned %d results for %d init tasks; a worker process died.",
                     length(init_flat), n_tasks))
      for (t in seq_len(n_tasks)) {
        el <- init_flat[[t]]
        if (inherits(el, "try-error"))
          stop(sprintf("Worker died computing init for rep %d fold %d: %s",
                       tasks$k[t], tasks$i[t], as.character(el)))
        if (inherits(el, "error"))
          stop(sprintf("Init failed for rep %d fold %d: %s",
                       tasks$k[t], tasks$i[t], conditionMessage(el)))
      }
      init_mat <- vector("list", ncv)
      for (k in seq_len(ncv)) init_mat[[k]] <- vector("list", nfolds)
      for (t in seq_len(n_tasks)) {
        init_mat[[tasks$k[t]]][[tasks$i[t]]] <- init_flat[[t]]
      }
    }

    # Phase 2: fan out nfolds * ncv * n_pairs tasks.
    pair_tasks <- expand.grid(j = seq_len(n_pairs), i = seq_len(nfolds),
                              k = seq_len(ncv),
                              KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    n_pair_tasks <- nrow(pair_tasks)

    pair_res_flat <- run_fn(seq_len(n_pair_tasks), function(t) {
      k_t <- pair_tasks$k[t]; i_t <- pair_tasks$i[t]; j_t <- pair_tasks$j[t]
      train_idx <- which(foldid[, k_t] != i_t)
      test_idx  <- which(foldid[, k_t] == i_t)
      init_vec  <- init_mat[[k_t]][[i_t]]
      list(
        k = k_t, i = i_t, j = j_t,
        test_idx = test_idx,
        result   = tryCatch(
          .cv_fold_path(x, y, train_idx, test_idx, grid[j_t, , drop = FALSE],
                        control, tau, init = init_vec),
          error = function(e) e
        )
      )
    })

    if (length(pair_res_flat) != n_pair_tasks)
      stop(sprintf("mclapply returned %d results for %d (rep,fold,pair) tasks; a worker process died.",
                   length(pair_res_flat), n_pair_tasks))
    for (t in seq_len(n_pair_tasks)) {
      el <- pair_res_flat[[t]]
      if (inherits(el, "try-error"))
        stop(sprintf("Worker died for rep %d fold %d pair %d: %s",
                     pair_tasks$k[t], pair_tasks$i[t], pair_tasks$j[t], as.character(el)))
    }

    # Stitch per-pair results back into per-fold shape.
    fold_results   <- vector("list", ncv)
    test_idx_store <- vector("list", ncv)
    fold_init_out  <- vector("list", ncv)
    for (k in seq_len(ncv)) {
      fold_results[[k]]   <- vector("list", nfolds)
      test_idx_store[[k]] <- vector("list", nfolds)
      fold_init_out[[k]]  <- vector("list", nfolds)
    }

    # Pre-fill test_idx_store and fold_init_out.
    for (t in seq_len(n_tasks)) {
      k_t <- tasks$k[t]; i_t <- tasks$i[t]
      test_idx_store[[k_t]][[i_t]] <- which(foldid[, k_t] == i_t)
      fold_init_out[[k_t]][[i_t]]  <- init_mat[[k_t]][[i_t]]
    }

    # Build per-fold synthetic result objects from the per-pair pieces.
    for (k in seq_len(ncv)) {
      for (i in seq_len(nfolds)) {
        n_test_ki <- length(test_idx_store[[k]][[i]])

        lp_out    <- matrix(NA_real_,    nrow = n_test_ki, ncol = n_pairs)
        iters_out <- rep(NA_integer_,    n_pairs)
        conv_out  <- rep(NA,             n_pairs)
        errs_list <- vector("list",      n_pairs)

        for (t in seq_len(n_pair_tasks)) {
          el <- pair_res_flat[[t]]
          if (el$k != k || el$i != i) next
          j_t <- el$j
          r_t <- el$result
          if (inherits(r_t, "error")) {
            errs_list[[j_t]] <- data.frame(
              pair    = grid$pair[j_t],
              s0      = grid$s0[j_t],
              s1      = grid$s1[j_t],
              message = conditionMessage(r_t),
              stringsAsFactors = FALSE
            )
          } else {
            # r_t is from .cv_fold_path with a 1-row grid; all vectors length 1.
            lp_out[, j_t]    <- r_t$lp[, 1L]
            iters_out[j_t]   <- r_t$iterations[1L]
            conv_out[j_t]    <- r_t$conv[1L]
            if (!is.null(r_t$coefs) && isTRUE(control$keep_coefs)) {
              # coef_list from a 1-row grid has 1 element; insert at position j_t
              # — handled later via fold_results coef_list stitching if needed.
            }
          }
        }

        non_null <- Filter(Negate(is.null), errs_list)
        errs_df  <- if (length(non_null) > 0L) {
          do.call(rbind, non_null)
        } else {
          data.frame(pair = integer(0), s0 = numeric(0), s1 = numeric(0),
                     message = character(0), stringsAsFactors = FALSE)
        }

        fold_results[[k]][[i]] <- list(
          lp         = lp_out,
          init       = init_mat[[k]][[i]],
          iterations = iters_out,
          conv       = conv_out,
          errors     = errs_df,
          n_failed   = nrow(errs_df),
          coefs      = NULL   # keep_coefs not supported on wide path for now
        )
      }
    }
  }

  # ---- 5. Assemble scores and errors ----
  pooled_mat     <- matrix(NA_real_,   nrow = n_pairs, ncol = ncv)
  fold_score_arr <- array(NA_real_,    dim  = c(n_pairs, nfolds, ncv))
  iters_arr      <- array(NA_integer_, dim  = c(n_pairs, nfolds, ncv))
  conv_arr       <- array(NA,          dim  = c(n_pairs, nfolds, ncv))
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
        # Scatter lp, iterations, and convergence into rep-level stores.
        lp_rep[test_idx, ] <- result$lp
        iters_arr[, i, k]  <- result$iterations
        conv_arr[, i, k]   <- result$conv

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

  # n_not_converged[j]: count of (rep, fold) combinations where the fit
  # returned successfully but conv == FALSE.  Failed fits (NA conv) are excluded.
  n_not_converged_vec <- apply(conv_arr, 1L,
                               function(v) sum(v == FALSE, na.rm = TRUE))

  tuning <- data.frame(
    s0              = grid$s0,
    s1              = grid$s1,
    score_mean      = score_mean,
    score_sd        = score_sd,
    score_fold_sd   = score_fold_sd,
    n_failed_folds  = n_failed_vec,
    n_not_converged = n_not_converged_vec,
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

  # ---- 7. Failure and non-convergence reporting ----
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

  # Non-convergence is a separate, less-severe condition: the fit returned a
  # usable result but the EM did not satisfy the convergence criterion.
  total_nc <- sum(n_not_converged_vec)
  if (total_nc > 0L) {
    n_pairs_nc  <- sum(n_not_converged_vec > 0L)
    total_fits  <- sum(!is.na(conv_arr))
    warning(
      sprintf(
        paste0("%d of %d fit(s) did not converge across %d pair(s). ",
               "Consider increasing maxit or relaxing epsilon."),
        total_nc, total_fits, n_pairs_nc
      ),
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
    iterations  = iters_arr,
    convergence = conv_arr,
    fold_inits  = fold_init_out,
    eval_time   = list(value = tau, source = tau_source),
    errors      = errors_df,
    folds       = folds,
    control     = control,
    grid        = grid,
    timing      = list(elapsed = elapsed, workers = n_workers,
                       parallel = isTRUE(control$parallel),
                       dispatch = dispatch_label)
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

  cat(sprintf("  fits     : %d failed, %d not converged\n",
              sum(x$tuning$n_failed_folds),
              sum(x$tuning$n_not_converged)))

  elapsed_s <- x$timing$elapsed["elapsed"]
  cat(sprintf("  elapsed  : %.2f s\n", elapsed_s))

  invisible(x)
}
