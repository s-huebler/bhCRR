# sim/chpc/R/parse.r -- CHPC parse-job wrapper
#
# THIS FILE IS A THIN WRAPPER. All parsing logic lives in sim/R/parse.r.
# This file reads env vars and calls parse_results(). Nothing else.
#
# REPLACES
#   sim/chpc/R/parse_batch.R  (525 lines; entire body gone)
#   The old out_dir was file.path(repo_root, "Sim", "chpc_results", scenario).
#   The new out_dir is file.path(REPO_ROOT, "sim", "results", SCENARIO, "parsed"),
#   matching the layout sim/R/parse.r documents.
#
# CALL SITE
#   parse.slurm  (invoked as: Rscript $REPO_ROOT/sim/chpc/R/parse.r)
#
# ENVIRONMENT VARIABLES
#   REPO_ROOT        required; path to the bhCRR repo checkout
#   SCRATCH_RUN_DIR  required; directory containing p<P>/ subdirs (runs_root)
#   SCENARIO         required; scenario slug, e.g. "scenarioA"
#
# LOAD MODE
#   Always calls load_sim("cluster") -- never reads BHCRR_LOAD. Same rationale
#   as run_array.r: deliberate defense in depth against a missing export.
#
# DEPENDS ON
#   sim/R/load.r    repo_root(), load_sim()
#   sim/R/parse.r   parse_results()  (sourced via load_sim)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

get_env <- function(name) {
  v <- Sys.getenv(name, unset = NA_character_)
  if (is.na(v) || !nzchar(v))
    stop("Required environment variable not set: ", name, call. = FALSE)
  v
}

# ---------------------------------------------------------------------------
# Read required env vars
# ---------------------------------------------------------------------------

repo_root_val   <- get_env("REPO_ROOT")
scratch_run_dir <- get_env("SCRATCH_RUN_DIR")
scenario_name   <- get_env("SCENARIO")

out_dir <- file.path(repo_root_val, "sim", "results", scenario_name, "parsed")

# ---------------------------------------------------------------------------
# Load bhCRR and sim modules (always cluster mode -- never read BHCRR_LOAD)
# ---------------------------------------------------------------------------

source(file.path(repo_root_val, "sim", "R", "load.r"))
load_sim("cluster")

# ---------------------------------------------------------------------------
# Parse
# ---------------------------------------------------------------------------

message("=== bhCRR parse.r ===")
message(sprintf("  repo_root    : %s", repo_root_val))
message(sprintf("  scenario     : %s", scenario_name))
message(sprintf("  runs_root    : %s", scratch_run_dir))
message(sprintf("  out_dir      : %s", out_dir))

parse_results(
  runs_root     = scratch_run_dir,
  scenario_name = scenario_name,
  out_dir       = out_dir,
  manifest_dir  = scratch_run_dir
)
