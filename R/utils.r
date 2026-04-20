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
#' @examples
#' NULL %||% 42      # 42
#' "x"  %||% 42     # "x"
#'
#' @export
`%||%` <- function(a, b) if (!is.null(a)) a else b


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


#' Permute the columns of each row independently
#'
#' For each row of \code{df}, independently shuffles the column values.
#' Used to break between-column correlations while preserving the marginal
#' distribution of each row.
#'
#' @param df A numeric data frame (or object coercible to a numeric matrix).
#' @param seed Optional integer seed passed to \code{\link{set.seed}} for
#'   reproducibility.  \code{NULL} (default) leaves the RNG state unchanged.
#'
#' @return A data frame with the same dimensions and column/row names as
#'   \code{df}, with each row's values independently permuted across columns.
#'
#' @examples
#' set.seed(1)
#' df <- as.data.frame(matrix(1:12, nrow = 3))
#' permute_rows(df, seed = 42)
#'
#' @export
permute_rows <- function(df, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  m <- as.matrix(df)
  permuted_mat <- t(apply(m, 1, function(row) row[sample.int(length(row))]))
  colnames(permuted_mat) <- colnames(m)
  rownames(permuted_mat) <- rownames(m)
  as.data.frame(permuted_mat)
}
