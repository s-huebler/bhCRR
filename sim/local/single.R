# sim/local/single.R -- interactive driver: one run, objects left in session
#
# OWNS
#   (no functions; this is a top-level script sourced interactively)
#
# WHAT IT DOES
#   Sources sim/R/load.r, reads a scenario file from sim/scenarios/, calls
#   run_one() for a single (p, run_id) pair, and assigns the result to
#   `result` in the global environment so all fields are available for
#   interactive inspection. Does NOT write any files.
#
# INPUTS (set before sourcing, or passed as variables)
#   SCENARIO_FILE  character(1)  path to a scenarios/*.R file;
#                                default: sim/scenarios/scenarioA.R
#   P              integer(1)    number of predictors; default: 25
#   RUN_ID         integer(1)    run index; default: 1
#   LOAD_MODE      character(1)  "dev" or "cluster"; default: "dev"
#
# RETURNS
#   Assigns `result` (output of run_one()) to the global environment.
#
# DEPENDS ON
#   sim/R/load.r, sim/R/run_one.r (transitively all of sim/R/)
#
# REPLACES
#   sim/_legacy/code/single_run_complex.R  (interactive use case)
#   sim/_legacy/code/single_run.R          (idem, simpler variant)
#
# NOT YET IMPLEMENTED
