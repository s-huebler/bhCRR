# sim/scenarios/scenarioA.R -- baseline scenario definition
#
# OWNS
#   (no functions; this file is sourced by scenario_from_file() in scenario.r)
#
# WHAT IT DOES
#   Returns a plain named list assigned to `scenario` that specifies all DGP
#   parameters for Scenario A. scenario_from_file() reads this list, validates
#   it against the scenario() constructor, and passes it to run_one().
#
# CONTRACT
#   This file MUST assign a list named `scenario`.
#   It MUST NOT source, library(), or otherwise produce side effects.
#   Every key below is required by dgp.r::sim_data(); removing or renaming a
#   key is a breaking change that requires a new scenario name.
#
# DGP KEYS (transcribed from sim/chpc/config.sh defaults)
#   name            character(1)  scenario slug used in file names and results
#   nobs            integer(1)    sample size
#   zero_gap_target numeric(1)    clinically-relevant minimum treatment effect;
#                                 passed to bhcrr_autotune() via fit_bhcrr()
#   beta1_active    numeric       active cause-1 coefficient vector
#   beta2_active    numeric or "neg_beta1"
#                                 cause-2 coefficients; "neg_beta1" means
#                                 -beta1_active
#   active_block    character     block label for each active predictor,
#                                 aligned with beta1_active
#   block_props     named numeric  proportion of predictors in each block type
#   block_rho       named numeric  within-block correlation per block type
#   latent_type     character(1)  latent variable family ("gaussian")
#   latent_q        named numeric  latent quantile per block type
#   cens_rate       numeric(1)    Exp rate for random censoring (0 = none)
#   u_min           numeric(1)    lower bound for observation window
#   u_max           numeric(1)    upper bound for observation window
#
# DEPENDS ON
#   sim/R/scenario.r  (scenario_from_file sources this file)
#
# NOT YET IMPLEMENTED
