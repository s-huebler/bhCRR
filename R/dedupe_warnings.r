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
