# sim/R/parse.r -- many result files -> CSV and .Rdata summary tables
#
# OWNS
#   parse_results(results_dir, scenario_name, out_dir)
#                       find all .Rdata run files under
#                       results_dir/<scenario_name>/runs/, validate their
#                       $meta and $scores shapes, stack into summary tables,
#                       write betas_p<P>_<scenario>.csv/.Rdata,
#                       model_performance_p<P>_<scenario>.csv,
#                       false_positives_p<P>_<scenario>.csv, and the
#                       stacked model_performance_<scenario>.csv when >1 p
#
# INPUTS
#   results_dir    character(1)  path to the results tree, e.g.
#                                sim/results/ (local) or the REPO_ROOT-relative
#                                equivalent on CHPC
#   scenario_name  character(1)  scenario slug, e.g. "scenarioA"
#   out_dir        character(1)  where to write the summary files; defaults to
#                                results_dir/<scenario_name>/parsed/
#
# RETURNS
#   invisible(list) with paths to all files written; called for its side
#   effect of writing CSV and .Rdata outputs
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
