# sim/R/score.r -- compute all evaluation metrics for one run
#
# OWNS
#   score_run(fit_result, comparator_results, x, y, scenario, beta1, block_map)
#                       compute all metrics for bhCRR and each comparator;
#                       all metrics are computed on every run because they
#                       are cheap relative to fitting, which eliminates the
#                       need for separate "wolbers", "wolbers_cv", etc. script
#                       variants that existed in the legacy code;
#                       per-block metrics are computed here, not in parse.r,
#                       so that parse.r is a pure stacker
#
# INPUTS
#   fit_result          list    output of fit.r::fit_bhcrr()
#   comparator_results  list    output of comparators.r::fit_comparators();
#                               must contain elements "fg_oracle", "lasso",
#                               "alasso", "scad", "mcp"
#   x                   matrix n x p
#   y                   matrix n x 2  (time, status)
#   scenario            list    validated scenario spec; used for zero_gap_target
#   beta1               numeric p   data$beta1, the cause-1 true coefficient
#                                   vector (full length, zeros for inactives)
#   block_map           data.frame  data.frame from $meta$block_map; one row
#                                   per predictor with columns name, block,
#                                   block_rho, designed, active1, active2
#
# RETURNS
#   list(overall = data.frame, by_block = data.frame)
#
#   $overall  one row per model, columns:
#     model         character   see model list below
#     wolbers_c     numeric     Wolbers IPCW C-index at eval_time
#     deviance      numeric     partial log-likelihood deviance
#     tp, fp, tn, fn            integer; threshold: coef != 0
#     precision     numeric     tp / (tp + fp)
#     recall        numeric     tp / (tp + fn)   [= sensitivity]
#     specificity   numeric     tn / (tn + fp)
#     mse_coef      numeric     mean squared error of coef vs beta1
#     fit_seconds   numeric     wall-clock seconds for THIS model's fit only;
#                               stored as numeric seconds, not difftime, because
#                               difftime units switch between seconds and minutes
#                               with magnitude and would mix units silently when
#                               stacked across runs
#
#   $by_block  one row per (model x block), columns:
#     model, block, block_rho,
#     tp, fp, tn, fn, precision, recall, mse_coef
#
#   CRITICAL: wolbers_c and deviance appear in $overall ONLY. They are
#   whole-model quantities; a C-index computed over a subset of predictors is
#   not a block-stratified C-index -- it is meaningless. Only selection and
#   coefficient-recovery metrics stratify by block.
#
#   Model list (six rows, in this order):
#     "bhcrr"     the primary method
#     "fg_oracle" Fine-Gray fit on the TRUE active subset (known beta1 != 0);
#                 an upper bound on achievable performance, not a competing
#                 method; included so bhCRR's gap to oracle can be reported
#     "lasso"     penalised Fine-Gray, lambda by CV
#     "alasso"    adaptive LASSO, first-stage LASSO weights
#     "scad"      SCAD-penalised Fine-Gray
#     "mcp"       MCP-penalised Fine-Gray
#
# DEPENDS ON
#   sim/R/load.r
#   bhCRR::wolbers_c
#
# REPLACES
#   Scoring logic scattered across sim/_legacy/code/batch_run_complex.R,
#   batch_run_wolbers.R, batch_run_cv_wolbers.R, single_run.R (all variants)
#   sim/chpc/R/run_batch.R scoring section
#
# NOT YET IMPLEMENTED
