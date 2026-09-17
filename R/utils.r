# ---- operators ---------------------------------------------------------------

#' Null coalescing operator
#'
#' Returns \code{a} if it is not \code{NULL}, otherwise returns \code{b}.
#' Useful for supplying default values when a list element may be absent.
#'
#' @param a An object to test.
#' @param b The fallback value returned when \code{a} is \code{NULL}.
#'
#' @return \code{a} if \code{!is.null(a)}, else \code{b}.
#'
#' @name grapes-or-or-grapes
#' @examples
#' NULL %||% 42      # 42
#' "x"  %||% 42     # "x"
#'
#' @export
`%||%` <- function(a, b) if (!is.null(a)) a else b


# ---- list utilities ----------------------------------------------------------

#' Recursively merge two named lists
#'
#' Merges \code{user} into \code{default} so that values in \code{user}
#' override corresponding values in \code{default}.  Keys present in
#' \code{user} but absent from \code{default} are added.  When both sides
#' have a list at the same key, the merge recurses into that sub-list rather
#' than replacing it wholesale.
#'
#' @param default A named list supplying baseline values.
#' @param user A named list of overrides, or \code{NULL} (in which case
#'   \code{default} is returned unchanged).
#'
#' @return A named list combining \code{default} and \code{user}.
#'
#' @examples
#' d <- list(a = 1, b = list(x = 10, y = 20))
#' u <- list(b = list(y = 99), c = 3)
#' sim_merge_lists(d, u)
#' # list(a = 1, b = list(x = 10, y = 99), c = 3)
#'
#' @export
sim_merge_lists <- function(default, user) {
  if (is.null(user)) return(default)
  if (!is.list(user)) stop("Expected a list, got: ", class(user)[1])
  out <- default
  for (nm in names(user)) {
    if (!nm %in% names(out)) {
      out[[nm]] <- user[[nm]]
    } else if (is.list(out[[nm]]) && is.list(user[[nm]])) {
      out[[nm]] <- sim_merge_lists(out[[nm]], user[[nm]])
    } else {
      out[[nm]] <- user[[nm]]
    }
  }
  out
}


# ---- distribution functions --------------------------------------------------

# Laplace (double exponential) distribution functions in base R

# Density: f(x) = 1/(2*b) * exp(-|x - mu| / b)
dlaplace <- function(x, mu = 0, b = 1, log = FALSE) {
  if (any(b <= 0)) stop("scale 'b' must be > 0")
  z <- abs(x - mu) / b
  logd <- -log(2 * b) - z
  if (log) logd else exp(logd)
}

leave_one_out_mean <- function(x) {
  (sum(x) - x) / (length(x) - 1)
}

# RETIRED 2026-09-17: plaplace(), qlaplace(), rlaplace() have no callers in R/,
# tests/, vignettes/, analyses/, Sim/, chpc/, or manuscript/. Commented out to
# reduce namespace clutter; kept for reference in case a future simulation study
# needs them.

# # CDF: for x < mu: 0.5 * exp((x - mu) / b)
# #      for x >= mu: 1 - 0.5 * exp(-(x - mu) / b)
# plaplace <- function(q, mu = 0, b = 1, lower.tail = TRUE, log.p = FALSE) {
#   if (any(b <= 0)) stop("scale 'b' must be > 0")
#   # ensure vector recycling like base R
#   q <- as.numeric(q)
#   mu <- as.numeric(mu)
#   b <- as.numeric(b)
#   # recycle
#   n <- max(length(q), length(mu), length(b))
#   q <- rep(q, length.out = n)
#   mu <- rep(mu, length.out = n)
#   b <- rep(b, length.out = n)
#
#   p <- numeric(n)
#   left <- q < mu
#   # left side
#   p[left] <- 0.5 * exp((q[left] - mu[left]) / b[left])
#   # right side (including equality)
#   p[!left] <- 1 - 0.5 * exp(-(q[!left] - mu[!left]) / b[!left])
#
#   if (!lower.tail) p <- 1 - p
#   if (log.p) log(p) else p
# }
#
# # Quantile function (inverse CDF)
# qlaplace <- function(p, mu = 0, b = 1, lower.tail = TRUE, log.p = FALSE) {
#   if (any(b <= 0)) stop("scale 'b' must be > 0")
#   p <- as.numeric(p)
#   mu <- as.numeric(mu)
#   b <- as.numeric(b)
#   n <- max(length(p), length(mu), length(b))
#   p <- rep(p, length.out = n)
#   mu <- rep(mu, length.out = n)
#   b <- rep(b, length.out = n)
#
#   if (log.p) p <- exp(p)
#   if (!lower.tail) p <- 1 - p
#
#   # validate p
#   if (any(p < 0 | p > 1, na.rm = TRUE)) stop("p must be in [0,1]")
#
#   q <- numeric(n)
#   # handle extremes
#   q[p == 0] <- -Inf
#   q[p == 1] <- Inf
#
#   mid <- (p > 0) & (p < 1)
#   if (any(mid)) {
#     pm <- p[mid]
#     mub <- mu[mid]
#     bb <- b[mid]
#     left <- pm < 0.5
#     # for p < 0.5: mu + b * log(2p)
#     q[mid][left] <- mub[left] + bb[left] * log(2 * pm[left])
#     # for p >= 0.5: mu - b * log(2*(1-p))
#     q[mid][!left] <- mub[!left] - bb[!left] * log(2 * (1 - pm[!left]))
#   }
#
#   q
# }
#
# # Random generation via inverse transform
# rlaplace <- function(n, mu = 0, b = 1) {
#   if (length(n) != 1 || n < 0) stop("'n' must be a non-negative integer scalar")
#   if (any(b <= 0)) stop("scale 'b' must be > 0")
#   n <- as.integer(n)
#   u <- stats::runif(n)
#   qlaplace(u, mu = mu, b = b)
# }


# ---- condition handling ------------------------------------------------------

#' Run an expression, emitting each unique warning at most once
#'
#' Internal utility. Wraps \code{expr} in a calling handler that records the
#' message of each warning it sees; identical warning messages raised later in
#' the same call are muffled. Used to keep repeated calls to
#' \code{fastcmprsk::Crisk()} (inside the EM loop and CV refits) from flooding
#' the console with the same "cencode is not a valid value from fstatus" note.
#'
#' @param expr An expression to evaluate.
#'
#' @return The value of \code{expr}.
#' @keywords internal
.dedupe_warnings <- function(expr) {
  seen <- character()
  withCallingHandlers(
    expr,
    warning = function(w) {
      msg <- conditionMessage(w)
      if (msg %in% seen) {
        invokeRestart("muffleWarning")
      } else {
        seen <<- c(seen, msg)
      }
    }
  )
}
