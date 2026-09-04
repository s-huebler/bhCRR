# =============================================================================
# tuning_diagnostics.r
# -----------------------------------------------------------------------------
# Automates the manual pre-flight / post-flight checks used before choosing an
# (s0, s1) tuning grid, so an unmanned simulation can size and validate its
# grid without a human eyeballing the heatmap.
#
# Builds on the internal helpers in R/threshold.R (s1_range, s0_upper_range,
# solve_for_s0, zero_gap, c_plus) and on tune_ssl_psdh().
#
# Layers:
#   bhcrr_check_grid()     - validate a specific (s0_seq, s1_seq) grid
#   bhcrr_tune_preflight() - derive theta, bounds, caps + RECOMMEND a grid
#   bhcrr_tune_validate()  - validate a completed tuning data frame
#   bhcrr_autotune()       - orchestrator: preflight -> tune -> validate,
#                            auto-widen the grid and re-tune on a boundary max
# =============================================================================


# ---- small utilities --------------------------------------------------------

# A check row: name, pass (TRUE/FALSE), severity, human-readable message.
.chk <- function(name, pass, severity, message) {
  data.frame(check = name, pass = pass, severity = severity,
             message = message, stringsAsFactors = FALSE)
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
.safe_zero_gap <- function(n, s0, s1, theta) {
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
#' @seealso \code{\link{bhcrr_tune_preflight}}, \code{\link{zero_gap}}
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


# ---- layer 2: derive parameters + recommend a grid --------------------------

#' Pre-flight tuning diagnostics and grid recommendation
#'
#' Reproduces the manual workflow: estimate \code{theta} from the untuned
#' model, find the feasible \code{s1} band and clamp the "reasonable" range into
#' it, solve for the \code{s0} values that put the zero gap at the clinical
#' target and at the floor, and assemble a recommended \code{(s0_seq, s1_seq)}
#' grid. If \code{s0_seq}/\code{s1_seq} are supplied, that grid is also checked.
#'
#' @param mod A fitted \code{\link{fit_ssl_psdh}} object (uses
#'   \code{$coefficients} and \code{$x}).
#' @param beta_min Numeric. Clinically relevant minimum treatment effect (the
#'   target zero gap).
#' @param beta_floor Numeric. Absolute smallest meaningful zero gap.
#' @param reasonable_s1 Length-2 numeric. The \code{s1} range you'd like
#'   covered, pre-clamp. Default \code{c(0.5, 1)}. Ignored when \code{fix_s1}
#'   is supplied.
#' @param n_s0,n_s1 Integer. Number of points in each recommended sequence.
#'   \code{n_s1} is ignored when \code{fix_s1} is supplied.
#' @param fix_s1 Optional single numeric. Pin the slab scale to this value
#'   instead of searching it: the recommended grid becomes a pure \code{s0}
#'   sweep at \code{s1 = fix_s1}, and \code{reasonable_s1} / \code{n_s1} are
#'   ignored. The value is used as given even if it falls outside the feasible
#'   band \code{s1_range(n, p)} - that is reported as a check row, not silently
#'   clamped, since coercing \code{s1} is a deliberate act.
#' @param s0_seq,s1_seq Optional numeric vectors. A grid to validate instead of
#'   / alongside the recommendation. When \code{fix_s1} is supplied, any
#'   \code{s1_seq} given here is replaced by \code{fix_s1}.
#' @param n,p Integer. Observations and predictors; default to the dimensions of
#'   \code{mod$x}.
#'
#' @returns An object of class \code{"bhcrr_preflight"}: a list with the derived
#'   \code{theta}, \code{s1_bounds}, per-guess caps and solve points, a
#'   \code{recommended} grid, and the grid \code{checks}.
#'
#' @seealso \code{\link{bhcrr_check_grid}}, \code{\link{bhcrr_autotune}}
#' @export
bhcrr_tune_preflight <- function(mod,
                                 beta_min = 0.1, beta_floor = 0.01,
                                 reasonable_s1 = c(0.2, 1),
                                 n_s0 = 8L, n_s1 = 3L,
                                 fix_s1 = NULL,
                                 s0_seq = NULL, s1_seq = NULL,
                                 n = NROW(mod$x), p = NCOL(mod$x)) {

  fix_s1 <- .check_fix_s1(fix_s1)

  # 1. theta: proportion of non-zero coefficients in the untuned fit.
  est   <- mod$coefficients$Estimate
  theta <- sum(est != 0) / length(est)
  theta <- min(max(theta, 1e-4), 1 - 1e-4)   # keep strictly in (0, 1)

  # 2. s1 band, and clamp the reasonable range into it.
  # With fix_s1 the band is still computed (it is reported, and tells you
  # whether the coerced value is defensible) but it does not choose s1.
  s1_bounds <- s1_range(n, p)
  if (!is.null(fix_s1)) {
    fits     <- fix_s1 >= s1_bounds[1] && fix_s1 <= s1_bounds[2]
    s1_lo    <- fix_s1
    s1_hi    <- fix_s1
    s1_guess <- fix_s1
  } else {
    fits <- !(s1_bounds[1] > reasonable_s1[1] || reasonable_s1[2] > s1_bounds[2])
    if (fits) {
      s1_lo <- reasonable_s1[1]; s1_hi <- reasonable_s1[2]
    } else {
      s1_lo <- s1_bounds[2];     s1_hi <- s1_bounds[2]   # pin to feasible max
    }
    s1_guess <- c(s1_lo, s1_hi)
  }

  # 3. Per-s1 caps and the s0 solve points (target + floor) for reference.
  caps         <- vapply(s1_guess, function(s1) s0_upper_range(n, s1), numeric(1))
  s0_at_target <- mapply(function(s1, cap) .safe_solve_s0(beta_min,   n, s1, theta, cap),
                         s1_guess, caps)
  s0_at_floor  <- mapply(function(s1, cap) .safe_solve_s0(beta_floor, n, s1, theta, cap),
                         s1_guess, caps)

  # 4. Recommended s1 sequence. A coerced s1 is exactly one point -- this is
  #    the branch that turns the recommendation into a pure s0 sweep.
  if (!is.null(fix_s1)) {
    rec_s1 <- fix_s1
  } else if (isTRUE(s1_lo == s1_hi)) {
    # Degenerate reasonable range -> small band just under the feasible max.
    lo <- max(s1_bounds[1], s1_hi * 0.8)
    rec_s1 <- sort(unique(seq(lo, s1_hi, length.out = n_s1)))
  } else {
    rec_s1 <- seq(s1_lo, s1_hi, length.out = n_s1)
  }

  # 5. Recommended s0 sequence: span the solve points (floor..target), pad a
  #    little, and clamp strictly below the tightest feasibility cap.
  anchors <- c(s0_at_target, s0_at_floor)
  anchors <- anchors[is.finite(anchors)]
  cap_min <- min(caps)
  if (length(anchors) == 0L) {
    hi <- cap_min * 0.98; lo <- cap_min * 0.10   # fall back to a sub-cap sweep
  } else {
    lo <- max(min(anchors) * 0.8, 1e-4)
    hi <- min(max(anchors) * 1.2, cap_min * 0.98)
  }
  if (hi <= lo) { lo <- cap_min * 0.10; hi <- cap_min * 0.98 }
  # Every pair must satisfy s1 > s0, so a coerced s1 is also a hard ceiling on
  # the s0 sweep. The feasibility cap is normally far below it, but clamp
  # anyway so an unusual (n, s1) combination cannot produce an empty grid.
  if (!is.null(fix_s1)) {
    hi <- min(hi, fix_s1 * 0.98)
    lo <- min(lo, hi * 0.5)
  }
  rec_s0 <- seq(lo, hi, length.out = n_s0)

  # 6. Validate the recommended grid, and any user-supplied grid. A coerced s1
  #    replaces whatever s1_seq was supplied -- it is the whole point of the
  #    argument, so honouring both would be incoherent.
  if (!is.null(fix_s1) && !is.null(s1_seq) &&
      !isTRUE(all.equal(sort(unique(s1_seq)), fix_s1))) {
    warning("fix_s1 overrides the supplied s1_seq; s1 is pinned to ",
            format(fix_s1), ".")
  }
  if (!is.null(fix_s1)) s1_seq <- fix_s1

  rec_check <- bhcrr_check_grid(n, theta, rec_s0, rec_s1, beta_min, beta_floor)
  sup_check <- if (!is.null(s0_seq) && !is.null(s1_seq))
    bhcrr_check_grid(n, theta, s0_seq, s1_seq, beta_min, beta_floor) else NULL

  # A coerced s1 outside the feasible band is reported, never silently moved.
  if (!is.null(fix_s1)) {
    fix_chk <- .chk("fix_s1_feasible", fits, "warn",
                    sprintf("Coerced s1 = %.4g is %s the feasible band [%.4g, %.4g].",
                            fix_s1, if (fits) "inside" else "OUTSIDE",
                            s1_bounds[1], s1_bounds[2]))
    rec_check$checks   <- rbind(rec_check$checks, fix_chk)
    rec_check$all_pass <- rec_check$all_pass && fits
    if (!is.null(sup_check)) {
      sup_check$checks   <- rbind(sup_check$checks, fix_chk)
      sup_check$all_pass <- sup_check$all_pass && fits
    }
  }

  structure(list(
    n = n, p = p, theta = theta,
    beta_min = beta_min, beta_floor = beta_floor,
    fix_s1 = fix_s1,
    s1_bounds = s1_bounds, s1_reasonable_fits = fits, s1_guess = s1_guess,
    s0_caps = caps, s0_at_target = s0_at_target, s0_at_floor = s0_at_floor,
    recommended = list(s0_seq = rec_s0, s1_seq = rec_s1),
    recommended_check = rec_check,
    supplied_check = sup_check,
    all_pass = rec_check$all_pass
  ), class = "bhcrr_preflight")
}

#' @export
print.bhcrr_preflight <- function(x, ...) {
  cat(sprintf("bhCRR tuning pre-flight  (n=%d, p=%d, theta=%.3f)\n",
              x$n, x$p, x$theta))
  cat(sprintf("  s1 feasible band : [%.4g, %.4g]  %s: %s\n",
              x$s1_bounds[1], x$s1_bounds[2],
              if (is.null(x$fix_s1)) "reasonable range fits" else "coerced s1 inside band",
              x$s1_reasonable_fits))
  if (!is.null(x$fix_s1))
    cat(sprintf("  s1 COERCED to    : %.4g (s0 sweep only; no s1 widening)\n",
                x$fix_s1))
  cat(sprintf("  s1 guesses       : %s\n",
              paste(sprintf("%.3f", x$s1_guess), collapse = ", ")))
  cat(sprintf("  s0 caps          : %s\n",
              paste(sprintf("%.4f", x$s0_caps), collapse = ", ")))
  cat(sprintf("  s0 @ target(%.3g): %s\n", x$beta_min,
              paste(sprintf("%.4f", x$s0_at_target), collapse = ", ")))
  cat(sprintf("  s0 @ floor(%.3g) : %s\n", x$beta_floor,
              paste(sprintf("%.4f", x$s0_at_floor), collapse = ", ")))
  cat(sprintf("  recommended s0   : seq(%.4f, %.4f, len=%d)\n",
              min(x$recommended$s0_seq), max(x$recommended$s0_seq),
              length(x$recommended$s0_seq)))
  cat(sprintf("  recommended s1   : %s\n",
              paste(sprintf("%.3f", x$recommended$s1_seq), collapse = ", ")))
  cat(sprintf("  grid checks      : %s\n",
              if (x$all_pass) "ALL PASS" else "SEE $recommended_check$checks"))
  invisible(x)
}


# ---- layer 3: validate a completed tuning run -------------------------------

#' Validate a tuning data frame from tune_ssl_psdh()
#'
#' Checks that (1) exactly one grid point attains the max score, (2) that
#' maximiser is interior to the grid (not on an edge - which would mean the
#' optimum likely lies outside the searched region), and (3) the zero gap at
#' the optimum is within a sane band around the clinical target.
#'
#' @param tuning Data frame returned by \code{\link{tune_ssl_psdh}}.
#' @param n Integer. Number of observations.
#' @param theta Numeric in (0, 1). Estimated non-zero proportion.
#' @param s0_seq,s1_seq Numeric vectors defining the grid that was searched.
#' @param beta_min Numeric. Clinically relevant minimum treatment effect.
#' @param gap_band Length-2 numeric. Multiplicative band around \code{beta_min}
#'   the optimum's zero gap should fall in. Default \code{c(0.5, 2)}.
#' @param s1_fixed Logical. Was \code{s1} coerced to a single value rather than
#'   searched? Defaults to \code{TRUE} when \code{s1_seq} holds one distinct
#'   value. A fixed \code{s1} is not a searched dimension, so the optimum can
#'   never be "on its boundary"; the edge test is skipped for it and only
#'   \code{s0} can trigger a widen.
#'
#' @returns A list: \code{checks}, \code{all_pass}, \code{best} (one-row data
#'   frame), \code{on_boundary}, \code{boundary_dims}, \code{gap_at_opt}.
#'
#' @seealso \code{\link{bhcrr_autotune}}
#' @export
bhcrr_tune_validate <- function(tuning, n, theta, s0_seq, s1_seq,
                                beta_min = 0.1, gap_band = c(0.5, 2),
                                s1_fixed = length(unique(s1_seq)) < 2L) {

  checks <- list()
  ok <- tuning[is.finite(tuning$score_mean), , drop = FALSE]

  if (nrow(ok) == 0L) {
    checks[[1L]] <- .chk("has_scores", FALSE, "critical",
                         "No finite scores; every fit failed.")
    return(list(checks = do.call(rbind, checks), all_pass = FALSE,
                best = NULL, on_boundary = NA, boundary_dims = character(0),
                gap_at_opt = NA_real_, s1_fixed = isTRUE(s1_fixed)))
  }

  # (1) single maximiser
  best_score <- max(ok$score_mean)
  max_rows   <- ok[ok$score_mean == best_score, , drop = FALSE]
  single_max <- nrow(max_rows) == 1L
  checks[[length(checks) + 1L]] <- .chk(
    "single_max", single_max, "warn",
    sprintf("%d grid point(s) tie for the max score %.4f.",
            nrow(max_rows), best_score))
  best <- max_rows[1, , drop = FALSE]

  # (2) interior optimum (not on a grid edge). A coerced s1 is a single point,
  # so it is trivially at both ends of its own "range" -- testing it would
  # report a boundary hit on every run and trigger widening that must not
  # happen. Skip the s1 edge test entirely in that case.
  s0_edge <- best$s0 <= min(s0_seq) || best$s0 >= max(s0_seq)
  s1_edge <- !isTRUE(s1_fixed) &&
    (best$s1 <= min(s1_seq) || best$s1 >= max(s1_seq))
  boundary_dims <- c(if (s0_edge) "s0", if (s1_edge) "s1")
  on_boundary <- length(boundary_dims) > 0L
  s1_note <- if (isTRUE(s1_fixed)) " [s1 fixed, not searched]" else ""
  checks[[length(checks) + 1L]] <- .chk(
    "interior_max", !on_boundary, "warn",
    if (on_boundary)
      sprintf("Max at edge in %s (s0=%.4f, s1=%.3f); optimum may lie outside grid.%s",
              paste(boundary_dims, collapse = "+"), best$s0, best$s1, s1_note)
    else sprintf("Max interior at s0=%.4f, s1=%.3f.%s", best$s0, best$s1, s1_note))

  # (3) zero gap at the optimum within band
  gap_at_opt <- .safe_zero_gap(n, best$s0, best$s1, theta)
  band <- sort(gap_band * beta_min)
  gap_ok <- is.finite(gap_at_opt) && gap_at_opt >= band[1] && gap_at_opt <= band[2]
  checks[[length(checks) + 1L]] <- .chk(
    "gap_reasonable", gap_ok, "warn",
    sprintf("Zero gap at optimum = %s; band [%.4g, %.4g] around target %.4g.",
            ifelse(is.finite(gap_at_opt), sprintf("%.4g", gap_at_opt), "NA"),
            band[1], band[2], beta_min))

  checks_df <- do.call(rbind, checks)
  all_pass <- all(checks_df$pass[checks_df$severity == "critical"]) &&
    all(checks_df$pass[checks_df$severity == "warn"])

  list(checks = checks_df, all_pass = all_pass, best = best,
       on_boundary = on_boundary, boundary_dims = boundary_dims,
       gap_at_opt = gap_at_opt, s1_fixed = isTRUE(s1_fixed))
}


# ---- orchestrator: preflight -> tune -> validate -> auto-widen --------------

# Extend a sequence outward on the low and/or high side by `add` points at the
# sequence's own spacing, clamped into (floor, ceil).
.widen_seq <- function(s, low, high, add = 2L, floor = 1e-4, ceil = Inf) {
  s <- sort(unique(s))
  step <- if (length(s) > 1L) diff(range(s)) / (length(s) - 1L) else max(s) * 0.1
  if (low)  s <- c(rev(seq(min(s) - step, by = -step, length.out = add)), s)
  if (high) s <- c(s, seq(max(s) + step, by = step, length.out = add))
  s <- s[s > floor & s < ceil]
  sort(unique(s))
}

# Did a widen actually change the grid? .widen_seq() returns a sorted, unique
# vector and can hand back exactly what it was given -- e.g. widening downward
# from a sequence already near the floor puts every new point at or below it,
# and the floor filter removes them all. Compare on the same footing (sorted,
# unique) so the mere act of sorting does not read as a change.
.seq_unchanged <- function(before, after) {
  before <- sort(unique(before))
  length(before) == length(after) && isTRUE(all.equal(before, after))
}

#' Autonomous tuning with pre-flight sizing and auto-widen retry
#'
#' Runs the whole loop for an unmanned simulation: derive/validate the grid,
#' tune, validate the result, and if the optimum lands on a grid edge, widen the
#' offending dimension(s) and re-tune (up to \code{max_widen} times). Never
#' halts; failures surface as warnings and in the returned \code{$flags}.
#'
#' A widen that adds no new grid points ends the loop rather than spending an
#' attempt re-fitting the identical grid: the attempt just completed stands as
#' the final one, and \code{$widen_stalled} is \code{TRUE}. This is the usual
#' outcome when the optimum sits on the low edge of an \code{s0} sequence that
#' already starts near the floor, since every candidate point below it is
#' clamped away.
#'
#' @param mod A fitted \code{\link{fit_ssl_psdh}} object.
#' @param beta_min,beta_floor Numeric. Clinical target and floor for the zero
#'   gap.
#' @param s0_seq,s1_seq Optional numeric vectors. Starting grid; if \code{NULL},
#'   uses the pre-flight recommendation.
#' @param fix_s1 Optional single numeric. Coerce the slab scale to this value
#'   and search \code{s0} only. It overrides \code{s1_seq}, and it pins the
#'   \code{s1} dimension out of the auto-widen loop entirely: \code{s0} still
#'   widens on a boundary maximum exactly as before, while \code{s1} is never
#'   widened, so the run's \code{s1} is the value you asked for and nothing
#'   else. Default \code{NULL} (search both dimensions).
#' @param nfolds,ncv,foldid Passed to \code{\link{tune_ssl_psdh}}.
#' @param max_widen Integer. Maximum number of auto-widen retries (0 disables).
#'   Applies to \code{s0} only when \code{fix_s1} is supplied.
#' @param gap_band Length-2 numeric. Band for \code{\link{bhcrr_tune_validate}}.
#' @param ... Additional arguments passed to \code{\link{tune_ssl_psdh}}.
#'
#' @returns An object of class \code{"bhcrr_autotune"}: a list with the
#'   \code{preflight} object, the final \code{tuning} data frame, the
#'   \code{validation} result, the \code{best} pair, the final \code{s0_seq} /
#'   \code{s1_seq}, the full \code{history}, \code{n_attempts},
#'   \code{widen_stalled} (did the loop stop because a widen added no new
#'   points?), \code{flags}, and \code{all_pass}.
#'
#' @seealso \code{\link{bhcrr_tune_preflight}}, \code{\link{tune_ssl_psdh}}
#' @export
bhcrr_autotune <- function(mod, beta_min = 0.1, beta_floor = 0.01,
                           s0_seq = NULL, s1_seq = NULL, fix_s1 = NULL,
                           nfolds = 10L, ncv = 2L, foldid = NULL,
                           max_widen = 1L, gap_band = c(0.5, 2), ...) {

  fix_s1   <- .check_fix_s1(fix_s1)
  s1_fixed <- !is.null(fix_s1)

  pf <- bhcrr_tune_preflight(mod, beta_min = beta_min, beta_floor = beta_floor,
                             fix_s1 = fix_s1,
                             s0_seq = s0_seq, s1_seq = s1_seq)
  n <- pf$n; theta <- pf$theta

  if (is.null(s0_seq)) s0_seq <- pf$recommended$s0_seq
  if (is.null(s1_seq)) s1_seq <- pf$recommended$s1_seq

  # A coerced s1 wins over anything passed in s1_seq (the preflight has already
  # warned if the two disagreed).
  if (s1_fixed) s1_seq <- fix_s1

  # Hard feasibility caps for widening: s0 must stay below the tightest cap,
  # s1 within the feasible band.
  s0_ceiling <- min(pf$s0_caps) * 0.98
  # s1 > s0 is required for every pair, so a coerced s1 also bounds how far s0
  # may widen upward.
  if (s1_fixed) s0_ceiling <- min(s0_ceiling, fix_s1 * 0.98)
  s1_floor   <- pf$s1_bounds[1]
  s1_ceiling <- pf$s1_bounds[2]

  if (!pf$recommended_check$all_pass && is.null(pf$supplied_check))
    warning("Pre-flight grid checks did not all pass; see $preflight.")

  history      <- list()
  flags        <- character(0)
  attempt      <- 0L
  widen_stalled <- FALSE
  repeat {
    attempt <- attempt + 1L
    tuning  <- tune_ssl_psdh(mod, s0_seq, s1_seq, nfolds = nfolds,
                             ncv = ncv, foldid = foldid, ...)
    val <- bhcrr_tune_validate(tuning, n, theta, s0_seq, s1_seq,
                               beta_min = beta_min, gap_band = gap_band,
                               s1_fixed = s1_fixed)
    history[[attempt]] <- list(s0_seq = s0_seq, s1_seq = s1_seq,
                               tuning = tuning, validation = val)

    if (!val$on_boundary || attempt > max_widen) break

    # Widen only the dimension(s) that hit an edge, in the direction of the edge.
    prev_s0_seq <- s0_seq
    prev_s1_seq <- s1_seq

    if ("s0" %in% val$boundary_dims) {
      low  <- val$best$s0 <= min(s0_seq)
      high <- val$best$s0 >= max(s0_seq)
      s0_seq <- .widen_seq(s0_seq, low, high, floor = 1e-4, ceil = s0_ceiling)
    }
    # s1 widening is unreachable when s1 is coerced (bhcrr_tune_validate never
    # reports an s1 boundary then); the guard is here so the intent survives
    # any future change to the boundary test.
    if (!s1_fixed && "s1" %in% val$boundary_dims) {
      low  <- val$best$s1 <= min(s1_seq)
      high <- val$best$s1 >= max(s1_seq)
      s1_seq <- .widen_seq(s1_seq, low, high, floor = s1_floor, ceil = s1_ceiling)
    }

    # Bail out if the widen produced no new grid points. Re-tuning would fit
    # the identical grid and reproduce the identical result at full CV cost,
    # so stop here and let this attempt stand as the final one. Happens when
    # every candidate point falls outside the feasibility clamp - most often
    # widening downward from an s0 sequence already sitting near the floor.
    if (.seq_unchanged(prev_s0_seq, s0_seq) &&
        .seq_unchanged(prev_s1_seq, s1_seq)) {
      widen_stalled <- TRUE
      flags <- c(flags, sprintf(
        paste0("attempt %d: boundary max in %s, but widening added no new ",
               "points (clamped at the %s bound) -- kept this attempt instead ",
               "of re-running an identical grid."),
        attempt, paste(val$boundary_dims, collapse = "+"),
        if (any(c(val$best$s0 <= min(prev_s0_seq),
                  !s1_fixed && val$best$s1 <= min(prev_s1_seq)))) "lower" else "upper"))
      break
    }

    flags <- c(flags, sprintf(
      "attempt %d: boundary max in %s -> widened and retried.",
      attempt, paste(val$boundary_dims, collapse = "+")))
  }

  final <- history[[length(history)]]
  if (final$validation$on_boundary && !widen_stalled)
    flags <- c(flags, "Optimum still on grid boundary after widening budget exhausted.")
  if (final$validation$on_boundary && widen_stalled)
    flags <- c(flags, "Optimum still on grid boundary; the grid cannot be widened further.")
  if (!final$validation$all_pass)
    warning("Final tuning did not pass all validation checks; see $validation$checks.")

  structure(list(
    preflight  = pf,
    tuning     = final$tuning,
    validation = final$validation,
    best       = final$validation$best,
    fix_s1     = fix_s1,
    s0_seq     = final$s0_seq,
    s1_seq     = final$s1_seq,
    history    = history,
    n_attempts = length(history),
    widen_stalled = widen_stalled,
    flags      = flags,
    all_pass   = final$validation$all_pass
  ), class = "bhcrr_autotune")
}

#' @export
print.bhcrr_autotune <- function(x, ...) {
  cat(sprintf("bhCRR autotune: %d attempt(s), %s%s\n",
              x$n_attempts,
              if (x$all_pass) "validation PASSED" else "validation FLAGS",
              if (is.null(x$fix_s1)) ""
              else sprintf("  [s1 coerced to %.4g, s0 sweep only]", x$fix_s1)))
  if (!is.null(x$best))
    cat(sprintf("  best: s0=%.4f  s1=%.3f  score=%.4f  zero gap=%.4g\n",
                x$best$s0, x$best$s1, x$best$score_mean, x$validation$gap_at_opt))
  if (length(x$flags))
    cat("  flags:\n", paste0("   - ", x$flags, collapse = "\n"), "\n", sep = "")
  invisible(x)
}
