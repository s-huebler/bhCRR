# sim/R/comparators.r -- LASSO, aLASSO, SCAD, MCP, and oracle comparator fits
#
# OWNS
#   fit_lasso(x, y)              cause-1 Fine-Gray LASSO via
#                                fastcmprsk::fastCrrp(penalty = "LASSO"),
#                                lambda selected by cv_comparator_lambda()
#   fit_alasso(x, y)             adaptive LASSO; LASSO first stage (ridge
#                                fallback) to get initial betas, then weights
#                                1 / pmax(abs(initial_beta_std), floor)^power,
#                                then fastCrrp(penalty = "LASSO", penalty.wt)
#   fit_scad(x, y)               SCAD-penalised Fine-Gray via
#                                fastcmprsk::fastCrrp(penalty = "SCAD")
#   fit_mcp(x, y)                MCP-penalised Fine-Gray via
#                                fastcmprsk::fastCrrp(penalty = "MCP")
#   fit_fg_oracle(x, y, beta1)   Fine-Gray fit on the TRUE active subset
#                                (columns where beta1 != 0); uses
#                                fastcmprsk::fastCrr with no penalisation;
#                                corresponds to the fg_true object in
#                                sim/_legacy/code/batch_run_complex.R:173 that
#                                was fitted there but never scored
#   cv_comparator_lambda(x, y, penalty, nfolds)
#                                fold-based lambda selection for one
#                                comparator; the legacy file
#                                (sim/_legacy/code/other_comparators_cv_section.R)
#                                carried its own CV loop independent of
#                                bhcrr_cv(); this function centralises that
#                                logic for all penalised comparators
#   fit_comparators(x, y, beta1) run all five (lasso, alasso, scad, mcp,
#                                fg_oracle); return a named list
#
# INPUTS
#   x      matrix n x p   design matrix; predictors are standardised ONCE
#                         before any comparator is called, and the same
#                         standardised matrix is passed to fit_bhcrr() --
#                         no method-to-method preprocessing differences appear
#   y      matrix n x 2   outcome (time, status); status in {0, 1, 2}
#   beta1  numeric p      full cause-1 coefficient vector (from sim_data());
#                         used by fit_fg_oracle() to identify the true active
#                         subset; also used by fit_alasso() as a fallback
#                         sanity check in tests
#   penalty character(1)  one of "LASSO", "SCAD", "MCP" (for cv_comparator_lambda)
#   nfolds  integer(1)    number of CV folds (for cv_comparator_lambda)
#
# RETURNS
#   fit_lasso / fit_alasso / fit_scad / fit_mcp / fit_fg_oracle:
#     list(coef = numeric p, converged = logical(1), fit_seconds = numeric(1),
#          meta = list(...))
#     failed fits carry coef = rep(NA_real_, p) so score.r never sees a
#     missing element; fit_seconds is numeric seconds (not difftime -- see
#     score.r timing note)
#   fit_comparators:
#     named list with elements "lasso", "alasso", "scad", "mcp", "fg_oracle",
#     each a list(coef, converged, fit_seconds, meta)
#
#   TIMING: each comparator is timed individually and returns its own
#   fit_seconds. The legacy code instrumented this wrongly: a single
#   other_time spanning all four fits
#   (sim/_legacy/code/batch_run_complex.R:492 and :589) cannot answer the
#   question of how bhCRR compares to each competitor separately.
#
# FUTURE WORK
#   The comparator set (lasso, alasso, scad, mcp, fg_oracle), the adaptive-
#   weight scheme (1/|beta|^power with ridge fallback), and the lambda
#   selection (fold-based CV within fit_comparators) are a deliberate minimal
#   first pass to get an end-to-end run working, and are expected to be
#   revisited once the simulation is running.
#
# NOTE ON LEGACY DESIGN
#   In sim/_legacy/code/other_comparators_cv_section.R, npredictors, x, and
#   sim_data were inherited from the calling environment rather than passed as
#   explicit arguments. Here they are explicit function arguments so
#   fit_comparators() can be called from run_one.r without environment leakage.
#
# DEPENDS ON
#   sim/R/load.r
#   bhCRR::                  (none; comparators use external packages only)
#   fastcmprsk::fastCrrp     (LASSO, aLASSO, SCAD, MCP penalised fits)
#   fastcmprsk::fastCrr      (fg_oracle unpenalised fit on active subset)
#
# REPLACES
#   sim/_legacy/code/other_comparators_cv_section.R  (entire file)
#   sim/_legacy/code/batch_run_cv_wolbers.R          comparator section
#   sim/_legacy/code/batch_run_comparison.R          (entire file)
#   sim/_legacy/code/single_run_comparison.R         (entire file)
#
# NOT YET IMPLEMENTED
