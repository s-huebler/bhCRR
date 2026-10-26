# sim/R/run_one.r -- scenario + run index -> one complete result list
#
# OWNS
#   run_one(scenario, p, run_id, load_mode)
#                       the single unit of work; draws data, fits bhCRR,
#                       fits comparators, scores all models, returns a result
#                       list; sim/local/batch.R and sim/chpc/R/run_array.r
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
#
#     $data         list(x, y, beta1, beta2)
#                     x:     matrix n x p  design matrix
#                     y:     matrix n x 2  (time, status), status in {0,1,2}
#                     beta1: numeric p     full cause-1 coefficient vector;
#                            active predictors at their named positions,
#                            zeros elsewhere
#                     beta2: numeric p     full cause-2 coefficient vector
#                            (same layout as beta1)
#                   block_map and block_signature from sim_data() are stored
#                   in $meta, not $data; one canonical home avoids drift
#
#     $fit          list          output of fit.r::fit_bhcrr()
#     $comparators  list          output of comparators.r::fit_comparators()
#
#     $scores       list(overall = data.frame, by_block = data.frame)
#                     output of score.r::score_run(); see score.r for the
#                     full column contracts of each frame
#
#     $meta         list
#       $meta$run_id          integer(1)
#       $meta$p               integer(1)
#       $meta$scenario_name   character(1)
#       $meta$git_sha         character(1)   git SHA of HEAD at run time
#       $meta$r_version       character(1)
#       $meta$bhcrr_version   character(1)
#       $meta$timing          proc.time() snapshot (wall, user, system) for
#                             the entire run_one() call
#       $meta$block_map       data.frame    one row per predictor; columns:
#                               name       character  column name ("V1", ...)
#                               block      character  block label
#                               block_rho  numeric    within-block rho
#                               designed   logical    in designed active set
#                               active1    logical    beta1 != 0
#                               active2    logical    beta2 != 0
#       $meta$block_signature character(1)  fingerprint of the realized
#                                           column-to-block map (see below)
#
#   BLOCK MAP PROVENANCE
#   The map uses NO RNG. It is a pure function of p, block_props, block_rho,
#   active_block and layout. The RNG sets a column's values, never its block.
#   So the map is stable across replicates and can be stored once per run and
#   joined downstream by predictor NAME.
#
#   It is stored rather than recomputed so that a later edit to block_props,
#   block_rho or layout cannot silently relabel results already on disk. This
#   is a provenance guarantee about time, not about scope: all reporting
#   happens within a single scenario, and the per-p directory layout already
#   keeps p groups apart.
#
#   block_signature exists because the scenario manifest pins the twelve DGP
#   parameters but not the generator CODE. Editing dgp.r between two
#   submissions to the same scenario is otherwise undetectable, since the
#   manifest's git SHA is written once at scenario creation and never
#   re-compared.
#
#   This exact shape is the contract that parse.r depends on. Any change to
#   $scores columns or $meta fields must be coordinated with parse.r.
#
# WHAT GETS SAVED
#   For now the result keeps as much as possible -- fitted objects, the
#   autotune object, comparator fits, x and y -- so that any diagnostic can
#   be run after the fact while the method is still being troubleshooted.
#   This is a known cost in file size and parse time, accepted on purpose.
#   The revisit trigger is the point where the presented summary statistics
#   are settled: average predictive accuracy of bhCRR against each comparator
#   and against fg_oracle, per-block precision and recall, and per-model
#   computation time. After that point the result can be trimmed to what those
#   require plus whatever diagnostics are still in use.
#   Do NOT trim anything during the current migration pass.
#   x and y live in $data; fitted objects must not carry redundant second
#   copies where that is avoidable.
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
#   sim/chpc/R/run_batch.R becomes sim/chpc/R/run_array.r, a thin wrapper
#     that reads env vars, resolves the scenario, and calls run_one() with
#     no workflow logic of its own; that file is not absorbed here
#
# NOT YET IMPLEMENTED
