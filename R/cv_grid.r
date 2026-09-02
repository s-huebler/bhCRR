#' Build and order the (s0, s1) hyperparameter grid
#'
#' Internal helper used by \code{bhcrr_cv()}.  Computes the Cartesian product
#' of \code{s0_seq} and \code{s1_seq}, drops pairs where \code{s1 <= s0}, and
#' orders the survivors in warm-start traversal order: unique \code{s1} values
#' in the order they appear in \code{s1_seq}, and within each group \code{s0}
#' ascending.  This matches the traversal in \code{\link{tune_ssl_psdh}}.
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
