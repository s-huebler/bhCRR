# sim/local/batch.R -- interactive driver: sweep over (p, run_id) pairs
#
# OWNS
#   (no functions; this is a top-level script sourced interactively)
#
# WHAT IT DOES
#   Sources sim/R/load.r, reads a scenario file from sim/scenarios/, then
#   calls run_one() for every combination of NPREDICTORS × RUN_START:RUN_END.
#   Each result is saved to
#     sim/results/<SCENARIO>/runs/result_p<p>_run<run_id>.Rdata
#   in the same layout that sim/chpc/R/run_batch.R produces, so parse.r can
#   process local and CHPC results identically. Does NOT accumulate results in
#   memory.
#
# INPUTS (set before sourcing, or passed as variables)
#   SCENARIO_FILE   character(1)   path to a scenarios/*.R file;
#                                  default: sim/scenarios/scenarioA.R
#   NPREDICTORS     integer(1+)    vector of p values; default: c(25)
#   RUN_START       integer(1)     first run index (inclusive); default: 1
#   RUN_END         integer(1)     last run index (inclusive); default: 10
#   RESULTS_DIR     character(1)   root of the results tree;
#                                  default: sim/results/
#   LOAD_MODE       character(1)   "dev" or "cluster"; default: "dev"
#
# OUTPUTS
#   sim/results/<SCENARIO>/runs/result_p<p>_run<run_id>.Rdata  (one per run)
#   Each file contains a single object `result` (output of run_one()).
#
# DEPENDS ON
#   sim/R/load.r, sim/R/run_one.r (transitively all of sim/R/)
#
# REPLACES
#   sim/_legacy/code/batch_run_complex.R   (local batch use case)
#   sim/_legacy/code/batch_run.R           (idem, simpler variant)
#   sim/_legacy/code/batch_run_wolbers.R   (idem, Wolbers scoring variant)
#
# NOT YET IMPLEMENTED
