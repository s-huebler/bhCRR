# ---- small utilities --------------------------------------------------------

# A check row: name, pass (TRUE/FALSE), severity, human-readable message.
.chk <- function(name, pass, severity, message) {
  data.frame(check = name, pass = pass, severity = severity,
             message = message, stringsAsFactors = FALSE)
}
s0_upper_range <- function(n, s1){
  lambda1 <- 1/s1

  lambda0 <- 2*sqrt(n)+lambda1
  1/lambda0-0.0001
}

s1_range <- function(n, p){
  max_s1 = p/sqrt(n)
  min_s1 = 1/(4*sqrt(n*log(p)))
  ret <- c(min_s1, max_s1)
  ret
}

c_plus <- function(n, s0, s1){
  lambda0 <- 1/s0
  lambda1 <- 1/s1

  rad <- 1-4*n/((lambda0-lambda1)^2)
  ret <- .5*(1+sqrt(rad))
  ret
}

zero_gap <- function(n, s0, s1, theta){
  lambda0 <- 1/s0
  lambda1 <- 1/s1

  1/(lambda0-lambda1)*log((1-theta)/theta*lambda0/lambda1*c_plus(n, s0, s1)/(1-c_plus(n, s0, s1)))
}



solve_for_s0 <- function(beta_min, n, s1, theta,
                         lower_bound = 0.001,
                         upper_bound = 0.05) {

  # The objective function: we want the result of this to be 0
  objective_fun <- function(s0) {
    zero_gap(n, s0, s1, theta) - beta_min
  }

  # Try to find the root
  # extendInt = "yes" tells R to automatically widen the search interval
  # if the root doesn't fall exactly between lower_bound and upper_bound.
  result <- uniroot(
    objective_fun,
    interval = c(lower_bound, upper_bound),
    extendInt = "yes"
  )

  # uniroot returns a list; we just extract the root (the s0 value)
  return(result$root)
}


# solve_for_s0() wraps uniroot(); it can throw or return nonsense for
# infeasible parameter combinations. Never let that kill an unmanned run.
.safe_solve_s0 <- function(beta_min, n, s1, theta, upper) {
  out <- tryCatch(
    solve_for_s0(beta_min = beta_min, n = n, s1 = s1, theta = theta,
                 lower_bound = 1e-4, upper_bound = upper),
    error = function(e) NA_real_, warning = function(w) NA_real_)
  if (length(out) != 1L || !is.finite(out) || out <= 0) return(NA_real_)
  out
}

# zero_gap() calls c_plus(), which is NaN when 4n > (lambda0 - lambda1)^2
# (i.e. s0 too large relative to the cap). Return NA rather than NaN so callers
# can distinguish "infeasible pair" cleanly.
safe_zero_gap <- function(n, s0, s1, theta) {
  g <- suppressWarnings(zero_gap(n, s0, s1, theta))
  if (length(g) != 1L || !is.finite(g)) return(NA_real_)
  g
}

# Validate a user-supplied fix_s1. NULL means "search s1 as usual"; anything
# else must be a single positive finite scalar, because the whole point of the
# argument is to collapse the s1 dimension to one point.
.check_fix_s1 <- function(fix_s1) {
  if (is.null(fix_s1)) return(NULL)
  if (!is.numeric(fix_s1) || length(fix_s1) != 1L ||
      !is.finite(fix_s1) || fix_s1 <= 0) {
    stop("fix_s1 must be a single positive finite number, or NULL to search s1.")
  }
  as.numeric(fix_s1)
}


# ---- layer 1: validate a specific grid --------------------------------------

#' Check an (s0, s1) tuning grid against the zero-gap heuristics
#'
#' Given the sample size, an estimate of the non-zero proportion \code{theta},
#' and a candidate grid, verify that the grid (a) contains only feasible pairs,
#' (b) brackets the clinically relevant minimum effect \code{beta_min} in
#' zero-gap space, and (c) does not sit entirely below the meaningful floor
#' \code{beta_floor}.
#'
#' @param n Integer. Number of observations.
#' @param theta Numeric in (0, 1). Estimated proportion of non-zero coefficients.
#' @param s0_seq,s1_seq Numeric vectors defining the candidate grid.
#' @param beta_min Numeric. Clinically relevant minimum treatment effect (the
#'   target zero gap).
#' @param beta_floor Numeric. Absolute smallest meaningful zero gap.
#'
#' @returns A list with elements \code{checks} (data frame of check rows),
#'   \code{all_pass} (logical), \code{gap_range}, \code{target}, \code{floor}.
#'
#' @seealso \code{bhcrr_tune_preflight}, \code{zero_gap}
#' @export
bhcrr_check_grid <- function(n, theta, s0_seq, s1_seq,
                             beta_min = 0.1, beta_floor = 0.01) {

  # All feasible (s1 > s0) pairs, mirroring tune_ssl_psdh()'s own filter.
  grid <- expand.grid(s0 = s0_seq, s1 = s1_seq)
  grid <- grid[grid$s1 > grid$s0, , drop = FALSE]

  checks <- list()

  if (nrow(grid) == 0L) {
    checks[[length(checks) + 1L]] <- .chk(
      "nonempty_grid", FALSE, "critical",
      "No pairs satisfy s1 > s0; tune_ssl_psdh() would error.")
    return(list(checks = do.call(rbind, checks), all_pass = FALSE,
                gap_range = c(NA, NA), target = beta_min, floor = beta_floor))
  }

  # Feasibility: s0 must sit below the cap s0_upper_range(n, s1) or c_plus is NaN.
  grid$cap <- vapply(grid$s1, function(s1) s0_upper_range(n, s1), numeric(1))
  grid$gap <- mapply(.safe_zero_gap, MoreArgs = list(n = n, theta = theta),
                     s0 = grid$s0, s1 = grid$s1)
  n_infeasible <- sum(grid$s0 >= grid$cap | is.na(grid$gap))

  checks[[length(checks) + 1L]] <- .chk(
    "feasible_pairs",
    n_infeasible == 0L,
    if (n_infeasible == nrow(grid)) "critical" else "warn",
    sprintf("%d of %d pairs are infeasible (s0 >= cap / NaN c_plus).",
            n_infeasible, nrow(grid)))

  gap_range <- if (all(is.na(grid$gap))) c(NA_real_, NA_real_) else
    range(grid$gap, na.rm = TRUE)

  # Does the grid bracket the clinical target in zero-gap space?
  brackets <- is.finite(gap_range[1]) &&
    gap_range[1] <= beta_min && beta_min <= gap_range[2]
  checks[[length(checks) + 1L]] <- .chk(
    "brackets_target", brackets, "warn",
    sprintf("Zero gap over grid spans [%.4g, %.4g]; clinical target %.4g is %s.",
            gap_range[1], gap_range[2], beta_min,
            if (isTRUE(brackets)) "inside" else "OUTSIDE"))

  # Is any part of the grid above the meaningful floor?
  above_floor <- is.finite(gap_range[2]) && gap_range[2] >= beta_floor
  checks[[length(checks) + 1L]] <- .chk(
    "above_floor", above_floor, "warn",
    sprintf("Max zero gap %.4g vs floor %.4g (%s).",
            gap_range[2], beta_floor,
            if (isTRUE(above_floor)) "ok" else "entire grid below floor"))

  checks_df <- do.call(rbind, checks)
  all_pass <- all(checks_df$pass[checks_df$severity == "critical"]) &&
    all(checks_df$pass[checks_df$severity == "warn"])

  list(checks = checks_df, all_pass = all_pass,
       gap_range = gap_range, target = beta_min, floor = beta_floor)
}








# target_zero_gap_grid <- function(s0_range,
#                                  s1_vals = 0.5,
#                                  theta_vals = 0.05,
#                                  n_val = 200,
#                                  target = 0.01) {

#   # 1. Main grid calculation (same as before)
#   grid <- expand.grid(s0 = s0_range,
#                       s1 = s1_vals,
#                       theta = theta_vals,
#                       n = n_val)
#
#   # Apply s0 upper range
#   grid$s0_max <- mapply(s0_upper_range, grid$n, grid$s1)
#
#   grid$Zero_Gap <- mapply(function(n, s0, s1, theta, s0_max) {
#     if (s0 > s0_max) return(NA_real_)
#     return(safe_zero_gap(n = n, s0 = s0, s1 = s1, theta = theta))
#   }, grid$n, grid$s0, grid$s1, grid$theta, grid$s0_max)
#
#   grid$Legend <- as.factor(paste0("s1: ", grid$s1))
#
#   # 2. Target grid for vertical drop lines
#   target_grid <- expand.grid(s1 = s1_vals, theta = theta_vals, n = n_val)
#   target_grid$Legend <- as.factor(paste0("s1: ", target_grid$s1))
#   target_grid$upper_bound <- mapply(s0_upper_range, target_grid$n, target_grid$s1)
#
#   # Find the maximum zero gap for each parameter curve to prevent solving for unreachable targets
#   max_gaps <- aggregate(Zero_Gap ~ s1 + theta, data = grid, FUN = function(x) max(x, na.rm = TRUE))
#   target_grid <- merge(target_grid, max_gaps, by = c("s1", "theta"), all.x = TRUE)
#
#   # 3. Solve for s0 conditionally using tryCatch to handle potential convergence errors
#   target_grid$target_s0 <- mapply(function(n, s1, theta, upper, max_gap) {
#     # Leave blank if the curve never reaches the target beta_min
#     if (is.na(max_gap) || max_gap < target) {
#       return(NA_real_)
#     }
#
#     tryCatch({
#       solve_for_s0(beta_min = target,
#                    n = n,
#                    s1 = s1,
#                    theta = theta,
#                    lower_bound = 0.0001,
#                    upper_bound = upper)
#     }, error = function(e) NA_real_)
#
#   }, target_grid$n, target_grid$s1, target_grid$theta, target_grid$upper_bound, target_grid$Zero_Gap)
#
#
#   ret <- target_grid |>
#     dplyr::group_by(s1, theta, n) |>
#     dplyr::summarise(max_s0 = unique(upper_bound),
#                      target_s0 = target_s0,
#                      .groups = 'drop') |>
#     dplyr::rowwise() |>
#     dplyr::mutate(max_s0_gap = safe_zero_gap(n, max_s0, s1, theta)) |>
#     dplyr::ungroup() |>
#     dplyr::mutate(target_gap = target) |>
#     dplyr::select(n, theta, s1, target_gap, target_s0, max_s0, max_s0_gap)
#
#   print(ret)
#
#   grid <- grid |>
#     dplyr::left_join(ret |>
#                        dplyr::select(n, theta, s1, target_s0),
#                      by = c("n", "theta", "s1"))
#
#   invisible(grid)
# }
#
#
# visualize_zero_gap_faceted <- function(s0_range,
#                                        s1_vals = 0.5,
#                                        theta_vals = 0.05,
#                                        n_val = 200,
#                                        target = 0.20) {
#
#
#   grid <- target_zero_gap_grid(s0_range,
#                                s1_vals,
#                                theta_vals,
#                                n_val,
#                                target)
#
#   # 4. Plotting
#   ggplot(grid, aes(x = s0, y = Zero_Gap, color = Legend, group = Legend)) +
#
#     # Dotted horizontal line using annotate instead of geom_segment
#     annotate("segment", x = 0, xend = max(grid$s0, na.rm = TRUE),
#              y = target, yend = target, color = "black", linetype = "dotted") +
#
#     geom_line(na.rm = TRUE) +
#     geom_point(na.rm = TRUE) +
#
#     # Vertical drop lines mapped to the intersection points from target_grid
#     geom_segment(aes(x = target_s0, xend = target_s0, y = target, yend = 0, color = Legend),
#                  linetype = "solid", na.rm = TRUE) +
#
#     facet_wrap(~ theta, labeller = label_both) +
#     labs(x = bquote(s[0]),
#          y = "Zero Gap",
#          color = "Parameters") +
#     theme_minimal() +
#     theme(strip.background = element_rect(fill = "grey90", color = NA))
# }
#
# visualize_zero_gap_faceted(s0_range = seq(0.001, .05, length = 20),
#                            s1_vals = c(0.2, 0.5, 0.8),
#                            theta_vals = c(0.05, 0.1, 0.25),
#                            n_val = 200,
#                            target = .05)
#
# visualize_zero_gap_faceted(s0_range = seq(0.001, .05, length = 20),
#                            s1_vals = c(.5, 1),
#                            theta_vals = c(0.0001, 0.001, 0.01, 0.1),
#                            n_val = 1309,
#                            target = .05)














