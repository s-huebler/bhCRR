# sim/R/score.r -- compute all evaluation metrics for one run
#
# OWNS
#   score_run(fit_result, comparator_results, x, y, scenario, beta_true)
#                       compute all metrics for bhCRR and each comparator;
#                       all metrics are computed on every run because they
#                       are cheap relative to fitting, which eliminates the
#                       need for separate "wolbers", "wolbers_cv", etc. script
#                       variants that existed in the legacy code
#
# INPUTS
#   fit_result          list    output of fit.r::fit_bhcrr()
#   comparator_results  list    output of comparators.r::fit_comparators()
#   x                   matrix n x p
#   y                   matrix n x 2   (time, status)
#   scenario            list    validated scenario spec; used for zero_gap_target
#                               and beta_true reconstruction
#   beta_true           numeric p   true coefficient vector for this run
#                               (drawn alongside the data in run_one.r)
#
# RETURNS
#   data.frame with one row and columns:
#     model              character   "bhcrr", "lasso", "alasso", "scad", "mcp"
#     wolbers_c          numeric     Wolbers IPCW C-index at eval_time
#     deviance           numeric     partial log-likelihood deviance
#     tp, fp, tn, fn     integer     true/false positive/negative counts
#                                    (threshold: coef != 0)
#     sensitivity        numeric     tp / (tp + fn)
#     specificity        numeric     tn / (tn + fp)
#     mse_coef           numeric     mean squared error vs beta_true
#   (one row per model; rbind into a p-row frame by run_one.r)
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
