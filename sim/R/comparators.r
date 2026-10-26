# sim/R/comparators.r -- LASSO, aLASSO, SCAD, MCP comparator fits
#
# OWNS
#   fit_lasso(x, y)     cause-1 Fine-Gray LASSO via cv.crr (FHcomp / crrlasso)
#   fit_alasso(x, y)    adaptive LASSO with LASSO first-stage weights
#   fit_scad(x, y)      SCAD-penalised Fine-Gray via ncvcr or equivalent
#   fit_mcp(x, y)       MCP-penalised Fine-Gray via ncvcr or equivalent
#   fit_comparators(x, y)   run all four; return a named list of coef vectors
#
# INPUTS
#   x  matrix n x p   design matrix (same object passed to fit_bhcrr)
#   y  matrix n x 2   outcome (time, status); status in {0, 1, 2}
#
# RETURNS
#   fit_lasso / fit_alasso / fit_scad / fit_mcp:
#     list(coef = numeric p, converged = logical(1), meta = list(...))
#   fit_comparators:
#     named list with elements "lasso", "alasso", "scad", "mcp", each a
#     list(coef, converged, meta); failed fits carry coef = rep(NA_real_, p)
#     so score.r never sees a missing element
#
# NOTE ON LEGACY DESIGN
#   In sim/_legacy/code/other_comparators_cv_section.R, npredictors, x, and
#   sim_data were inherited from the calling environment rather than passed as
#   explicit arguments. Here they are explicit function arguments so
#   fit_comparators() can be called from run_one.r without environment leakage.
#
# DEPENDS ON
#   sim/R/load.r
#   bhCRR:: (none; comparators use external packages only)
#   Packages: FHcomp or crrlasso (LASSO), ncvcr or equivalent (SCAD/MCP)
#
# REPLACES
#   sim/_legacy/code/other_comparators_cv_section.R  (entire file)
#   sim/_legacy/code/batch_run_cv_wolbers.R          comparator section
#   sim/_legacy/code/batch_run_comparison.R          (entire file)
#   sim/_legacy/code/single_run_comparison.R         (entire file)
#
# NOT YET IMPLEMENTED
