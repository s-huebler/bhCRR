# sim/R/fit.r -- tuning + final fit for one simulated dataset
#
# OWNS
#   fit_bhcrr(x, y, scenario, run_id)   run bhcrr_autotune() then refit at
#                                       the chosen (s0, s1); return a result
#                                       list suitable for score.r
#
# INPUTS
#   x         matrix n x p   design matrix from sim/R/dgp.r
#   y         matrix n x 2   outcome from sim/R/dgp.r (time, status)
#   scenario  list           validated scenario spec; used for zero_gap_target
#                            (-> beta_min) and to construct the cv control
#   run_id    integer(1)     run index, used only for logging
#
# RETURNS
#   list with fields:
#     $coef        numeric p     final coefficient vector (named, length p)
#     $s0          numeric(1)    tuned spike scale
#     $s1          numeric(1)    tuned slab scale
#     $autotune    bhcrr_autotune  full autotune object (for diagnostics)
#     $fit         ssl_psdh      final fit at chosen (s0, s1)
#     $converged   logical(1)
#     $all_pass    logical(1)    from autotune$all_pass
#     $flags       character     from autotune$flags
#     $n_attempts  integer(1)    from autotune$n_attempts
#
# DEPENDS ON
#   sim/R/load.r
#   bhCRR::fit_ssl_psdh, bhCRR::bhcrr_autotune, bhCRR::bhcrr_cv_control
#
# STABLE SURFACE CONTRACT
#   This file is the ONLY caller of bhcrr_autotune(). It is permitted to read
#   exactly five fields from the returned bhcrr_autotune object:
#     $best$s0, $best$s1, $all_pass, $flags, $n_attempts
#   Any other field is internal to the bridge and may change without notice.
#   See design/autotune-bridge.md §4 (STABLE SURFACE) for the rationale.
#
# REPLACES
#   sim/_legacy/code/batch_run_complex.R lines ~80-180 (tuning + fit loop)
#   sim/_legacy/code/single_run_complex.R (entire file; same logic, one run)
#   sim/chpc/R/run_batch.R becomes sim/chpc/R/run_array.r, a thin wrapper
#     that reads env vars, resolves the scenario, and calls run_one() with
#     no workflow logic of its own; the fitting logic it contained is
#     absorbed by this file and run_one.r, not by the CHPC wrapper
#
# NOT YET IMPLEMENTED
