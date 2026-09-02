#' Default initialization lambda path
#'
#' @keywords internal
.default_init_lam_path <- function() {
  10^seq(log10(0.1), log10(0.001), length = 10)
}

#' CV-LASSO initialization
#'
#' Selects starting coefficients by cross-validated C-index over a LASSO path
#' using \code{cv_fastCrrp_cpp}. Column selection is by index (the column
#' achieving the maximum CV C-index) rather than floating-point lambda equality,
#' so a path with duplicate lambda values always returns a vector, not a matrix.
#'
#' @param x Numeric matrix \eqn{n \times p}.
#' @param y Two-column matrix: column 1 times, column 2 status (0/1/2).
#' @param k Integer. Number of CV folds. Default \code{5}.
#' @param lambda_path Numeric vector of candidate lambdas.
#' @param tuning Character. CV metric passed to \code{cv_fastCrrp_cpp}.
#' @param eval_quantile Numeric in (0,1). Evaluation quantile for the metric.
#'
#' @return \code{list(init = <numeric p>, meta = <list>)}.
#' @keywords internal
.LASSO_cv <- function(x, y,
                      k             = 5,
                      lambda_path   = .default_init_lam_path(),
                      tuning        = "wolbers",
                      eval_quantile = 0.5) {
  t0  <- proc.time()[["elapsed"]]
  res <- cv_fastCrrp_cpp(x, y[, 1], y[, 2],
                         k             = k,
                         penalty       = "LASSO",
                         lambda_path   = lambda_path,
                         tuning        = tuning,
                         eval_quantile = eval_quantile)
  idx  <- which.max(res$cv_c_index)
  init <- as.numeric(res$full_model$coef[, idx])
  list(
    init = init,
    meta = list(
      method        = "LASSO_cv",
      lambda_min    = res$lambda_min,
      lambda_path   = lambda_path,
      cv_c_index    = res$cv_c_index,
      k             = k,
      tuning        = tuning,
      eval_quantile = eval_quantile,
      elapsed_sec   = proc.time()[["elapsed"]] - t0
    )
  )
}

#' BIC-LASSO initialization
#'
#' Fits a single LASSO path with \code{fastcmprsk::fastCrrp} and selects the
#' coefficient vector at the lowest-BIC lambda. Unlike \code{.LASSO_cv} this
#' requires no fold splitting, so it degrades more gracefully when cause-1
#' events are sparse.
#'
#' @param x Numeric matrix \eqn{n \times p}.
#' @param y Two-column matrix: column 1 times, column 2 status (0/1/2).
#' @param lambda_path Numeric vector of candidate lambdas.
#'
#' @return \code{list(init = <numeric p>, meta = <list>)}.
#' @keywords internal
.LASSO_bic <- function(x, y, lambda_path = .default_init_lam_path()) {
  n <- nrow(x)
  p <- ncol(x)
  t0 <- proc.time()[["elapsed"]]

  fit <- fastcmprsk::fastCrrp(
    fastcmprsk::Crisk(y[, 1L], y[, 2L], cencode = 0, failcode = 1) ~ x,
    penalty = "LASSO",
    lambda  = lambda_path
  )

  beta_path <- as.matrix(fit$coef)
  if (nrow(beta_path) != p && ncol(beta_path) == p) beta_path <- t(beta_path)
  if (nrow(beta_path) != p)
    stop("fastCrrp() returned a coefficient matrix of unexpected dimension: ",
         paste(dim(beta_path), collapse = " x "), " (expected ", p, " rows).")

  df_path <- colSums(beta_path != 0)
  loglik  <- as.numeric(fit$logLik)
  if (length(loglik) != ncol(beta_path))
    stop("fastCrrp() logLik has length ", length(loglik),
         " but the coefficient path has ", ncol(beta_path), " columns.")
  bic_path <- -2 * loglik + log(n) * df_path

  valid <- is.finite(bic_path) & (df_path > 0)
  conv  <- fit$converged
  if (!is.null(conv) && length(conv) == length(bic_path))
    valid <- valid & as.logical(conv)

  if (!any(valid)) {
    lam_str <- paste(format(lambda_path, digits = 3), collapse = ", ")
    stop("No finite, converged, non-null LASSO solution on the initialization ",
         "path (lambda = c(", lam_str, ")); cannot select an initial estimate ",
         "by BIC. Try a wider or denser lambda path.")
  }

  bic_sel <- bic_path
  bic_sel[!valid] <- Inf
  best_idx <- which.min(bic_sel)

  init <- as.numeric(beta_path[, best_idx])
  if (length(init) != p || any(!is.finite(init)))
    stop("BIC-selected initial estimate is not a finite length-", p, " vector.")

  lam_used <- fit$lambda %||% fit$lambda.path
  list(
    init = init,
    meta = list(
      method      = "LASSO_bic",
      index       = as.integer(best_idx),
      lambda      = if (is.null(lam_used)) NA_real_ else as.numeric(lam_used[best_idx]),
      bic         = as.numeric(bic_path[best_idx]),
      df          = as.integer(df_path[best_idx]),
      bic_path    = as.numeric(bic_path),
      df_path     = as.integer(df_path),
      elapsed_sec = proc.time()[["elapsed"]] - t0
    )
  )
}

#' Zero initialization
#'
#' Returns a zero vector of length \code{ncol(x)} as the starting point.
#'
#' @param x Numeric matrix \eqn{n \times p}.
#' @param y Two-column matrix (unused; present for a uniform call signature).
#'
#' @return \code{list(init = <numeric p>, meta = <list>)}.
#' @keywords internal
.zero_init <- function(x, y) {
  list(
    init = rep(0, ncol(x)),
    meta = list(method = "zero")
  )
}

#' Validate and normalise an initialization value
#'
#' Accepts either a bare numeric vector or a \code{list(init, meta)} and
#' returns the normalised form.  Errors if the result is not a finite numeric
#' vector of length \code{p}.
#'
#' @param value Output from an init function or a bare numeric vector.
#' @param p Expected length.
#' @param method_label Character label used in error messages.
#'
#' @return \code{list(init = <numeric p>, meta = <list or NULL>)}.
#' @keywords internal
.validate_init <- function(value, p, method_label) {
  if (is.numeric(value)) {
    out <- list(init = value, meta = NULL)
  } else if (is.list(value) && !is.null(value$init)) {
    out <- list(init = value$init, meta = value$meta)
  } else {
    stop("Init method '", method_label, "' must return a numeric vector or ",
         "list(init = <numeric>, meta = <list>); got class: ",
         paste(class(value), collapse = "/"))
  }

  v <- out$init
  if (!is.numeric(v))
    stop("Init method '", method_label, "' returned a non-numeric init vector ",
         "(class: ", paste(class(v), collapse = "/"), ").")
  if (length(v) != p)
    stop("Init method '", method_label, "' returned a vector of length ",
         length(v), " but x has ", p, " columns.")
  if (anyNA(v) || any(!is.finite(v)))
    stop("Init method '", method_label, "' returned an init vector containing ",
         "NA, NaN, or Inf.")
  out
}

#' Resolve an init_method argument to a function and a label
#'
#' @param init_method A function, a built-in name (\code{"LASSO_cv"},
#'   \code{"LASSO_bic"}, \code{"zero"}), or a character string naming a
#'   function in \code{envir}.
#' @param envir Environment to search when \code{init_method} is an
#'   unrecognised string. Defaults to the caller's environment.
#'
#' @return \code{list(fn = <function>, label = <character>)}.
#' @keywords internal
.resolve_init_method <- function(init_method, envir = parent.frame()) {
  builtin_names <- c("LASSO_cv", "LASSO_bic", "zero")
  builtins <- list(
    LASSO_cv  = .LASSO_cv,
    LASSO_bic = .LASSO_bic,
    zero      = .zero_init
  )

  if (is.function(init_method)) {
    return(list(fn = init_method, label = "custom"))
  }

  if (is.character(init_method) && length(init_method) == 1L) {
    nm <- init_method
    matched <- builtin_names[tolower(builtin_names) == tolower(nm)]
    if (length(matched) == 1L) {
      return(list(fn = builtins[[matched]], label = matched))
    }
    fn <- tryCatch(
      get(nm, mode = "function", envir = envir),
      error = function(e) NULL
    )
    if (!is.null(fn))
      return(list(fn = fn, label = "custom"))
    stop("'", nm, "' is not a built-in init method (",
         paste(builtin_names, collapse = ", "), ") and was not found as a ",
         "function in the calling environment.")
  }

  stop("'init_method' must be a function or a character string; got class: ",
       paste(class(init_method), collapse = "/"))
}
