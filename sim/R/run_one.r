# sim/R/run_one.r -- scenario + run index -> one complete result list
#
# OWNS
#   run_one(scenario, p, run_id, load_mode)
#                       the single unit of work; draws data, fits bhCRR,
#                       fits comparators, scores all models, returns a result
#                       list; sim/local/batch.R and sim/chpc/R/run_batch.R
#                       both call this function and contain no workflow logic
#                       of their own -- all simulation logic lives here
#
# INPUTS
#   scenario   list or character(1)   validated scenario spec (list) or path
#                                     to a scenarios/*.R file (character);
#                                     when a path is given, scenario_from_file()
#                                     is called internally
#   p          integer(1)             number of predictors for this run
#   run_id     integer(1)             run index; used for set.seed() and for
#                                     the $run_id field of the result
#   load_mode  character(1)           passed to load_bhcrr(); "dev" or "cluster"
#
# RETURNS
#   list with fields:
#     $run_id       integer(1)
#     $p            integer(1)
#     $scenario     list          the validated scenario spec used
#     $data         list(x, y, beta_true)
#     $fit          list          output of fit.r::fit_bhcrr()
#     $comparators  list          output of comparators.r::fit_comparators()
#     $scores       data.frame    output of score.r::score_run(); one row per
#                                 model, columns as specified in score.r
#     $meta         list          timing (proc.time()), bhCRR version, R version,
#                                 git SHA, run_id, p, scenario name
#
#   This exact shape is the contract that parse.r depends on. Any change to
#   $scores columns or $meta fields must be coordinated with parse.r.
#
# DEPENDS ON
#   sim/R/load.r
#   sim/R/scenario.r     scenario_from_file(), check_manifest()
#   sim/R/dgp.r          sim_data()
#   sim/R/fit.r          fit_bhcrr()
#   sim/R/comparators.r  fit_comparators()
#   sim/R/score.r        score_run()
#
# REPLACES
#   sim/_legacy/code/single_run_complex.R     (entire file; one-run case)
#   sim/_legacy/code/batch_run_complex.R      (the inner loop body)
#   sim/chpc/R/run_batch.R                    (CHPC adaptation of same)
#
# NOT YET IMPLEMENTED
