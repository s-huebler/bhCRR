#' Compute the initialisation vector for one (repetition, fold)
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
#'       computed by \code{\link{predict_from_ssl_psdh}}.  \code{NA} for
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
      lp_mat[, j] <- predict_from_ssl_psdh(fit_result,
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
