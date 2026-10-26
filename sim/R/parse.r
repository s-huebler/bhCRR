# sim/R/parse.r -- many result files -> CSV and .Rdata summary tables
#
# OWNS
#   parse_results(results_dir, scenario_name, out_dir)
#                       find all .Rdata run files under
#                       results_dir/<scenario_name>/runs/, apply two hard
#                       stops (see below), then stack and average the
#                       precomputed metrics from each run's $scores$overall
#                       and $scores$by_block; write summary CSVs and .Rdata
#
# INPUTS
#   results_dir    character(1)  path to the results tree, e.g. sim/results/
#   scenario_name  character(1)  scenario slug, e.g. "scenarioA"
#   out_dir        character(1)  where to write the summary files; defaults to
#                                results_dir/<scenario_name>/parsed/
#
# RETURNS
#   invisible(list) with paths to all files written; called for its side
#   effect of writing CSV and .Rdata outputs
#
#   Output files (all written to out_dir):
#     overall_p<P>_<scenario>.csv        stacked $scores$overall for one p
#     by_block_p<P>_<scenario>.csv       stacked $scores$by_block for one p
#     overall_<scenario>.csv             overall stacked across all p values
#     by_block_<scenario>.csv            by_block stacked across all p values
#
# DESIGN: PURE STACKER
#   All metrics arrive precomputed in each run's $scores$overall and
#   $scores$by_block; parse_results() stacks and averages them and derives
#   nothing new. Moving metric computation to score.r and out of parse.r
#   means that adding a metric never requires touching the parser.
#
#   Averaging per-block rates across runs is legitimate because the block map
#   is deterministic and the beta pattern is fixed: the number of true zeros
#   per block is identical in every run of a (scenario, p) group, which makes
#   the pooled rate exactly equal to the mean of per-run rates.
#
#   Block and BlockRho in the output tables come from each run's
#   $meta$block_map, joined by predictor NAME; they are never recomputed
#   from scenario parameters.
#
# HARD STOP 1 -- scenario isolation
#   If run files from more than one scenario name are found under the directory
#   being parsed, stop with an error. Scenario isolation is checked, not merely
#   assumed from the directory layout.
#
# HARD STOP 2 -- block signature consistency
#   If two runs in the same (scenario, p) group report different
#   $meta$block_signature values, stop with an error. A signature mismatch
#   means dgp.r was edited between submissions to the same scenario; mixing
#   those runs would produce silently incorrect per-block summaries.
#
# NOTE ON CONSOLIDATION
#   This file absorbs two hand-maintained copies of the same logic:
#     sim/_legacy/code/parse_batch_run.R  (local version)
#     sim/chpc/R/parse_batch.R            (CHPC version, diverged slightly)
#   Both are replaced by a single implementation here. sim/chpc/R/parse_batch.R
#   becomes a thin wrapper: source load.r, call parse_results(), done.
#
# DEPENDS ON
#   sim/R/load.r
#   bhCRR:: (none; parsing is pure R)
#
# REPLACES
#   sim/_legacy/code/parse_batch_run.R   (entire file)
#   sim/chpc/R/parse_batch.R             (entire file, becomes a wrapper)
#
# NOT YET IMPLEMENTED
