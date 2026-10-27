# sim/chpc/R/run_array.r -- CHPC array-task wrapper
#
# THIS FILE IS A THIN WRAPPER. All fitting, seeding, failure handling, and
# manifest logic live in sim/R/run_one.r and sim/local/batch.R. This file
# reads env vars, validates them, and calls batch_sweep(). Nothing else.
#
# REPLACES
#   sim/chpc/R/run_batch.R  (654 lines; entire body gone)
#
# CALL SITE
#   run_array.slurm  (invoked as: Rscript $REPO_ROOT/sim/chpc/R/run_array.r)
#
# ENVIRONMENT VARIABLES (required; no defaults for the first four)
#   REPO_ROOT         path to the bhCRR repo checkout
#   SCENARIO          scenario slug, e.g. "scenarioA"
#   SCRATCH_RUN_DIR   parent of p<P>/ subdirs; also used as manifest_dir
#   NPREDICTORS       comma-separated integer list, e.g. "25,50"
#   RUN_SUBSTART      first run index this task handles
#                     falls back to RUN_START if RUN_SUBSTART is unset
#   RUN_SUBEND        last run index this task handles
#                     falls back to RUN_END if RUN_SUBEND is unset
#
# DGP ENVIRONMENT VARIABLES (MUST NOT be set; scenario file is authoritative)
#   NOBS, ZERO_GAP_TARGET, BETA1_ACTIVE, BETA2_ACTIVE, ACTIVE_BLOCK,
#   BLOCK_PROPS, BLOCK_RHO, LATENT_TYPE, LATENT_Q, CENS_RATE, U_MIN, U_MAX
#
#   If any of these twelve are set in the environment, this script errors
#   immediately, naming them and the authoritative scenario file. A stale export
#   in the shell must not silently override what the scenario file says.
#
# LOAD MODE
#   Always calls load_sim("cluster") -- never reads BHCRR_LOAD. This is
#   deliberate defense in depth: if config.sh ever fails to export BHCRR_LOAD,
#   reading the env var would fall back to "dev", which compiles, and concurrent
#   array tasks would race on src/*.o and src/bhCRR.so. The wrappers never
#   trust the env var; they always name "cluster" explicitly.
#
# DEPENDS ON
#   sim/R/load.r        repo_root(), load_sim()
#   sim/local/batch.R   batch_sweep()
#   sim/R/scenario.r    (transitively, via batch_sweep)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

get_env <- function(name, default = NULL) {
  v <- Sys.getenv(name, unset = NA_character_)
  if (is.na(v) || !nzchar(v)) {
    if (is.null(default))
      stop("Required environment variable not set: ", name, call. = FALSE)
    return(default)
  }
  v
}

int_list <- function(s) as.integer(strsplit(s, ",", fixed = TRUE)[[1]])

# ---------------------------------------------------------------------------
# Read required env vars
# ---------------------------------------------------------------------------

repo_root_val   <- get_env("REPO_ROOT")
scenario_name   <- get_env("SCENARIO")
scratch_run_dir <- get_env("SCRATCH_RUN_DIR")
npredictors_str <- get_env("NPREDICTORS")

npredictors_grid <- int_list(npredictors_str)
run_start <- as.integer(get_env("RUN_SUBSTART", get_env("RUN_START", "1")))
run_end   <- as.integer(get_env("RUN_SUBEND",   get_env("RUN_END",   "10")))

# ---------------------------------------------------------------------------
# Guard: DGP env vars must NOT be set; the scenario file is authoritative
# ---------------------------------------------------------------------------

.dgp_vars <- c("NOBS", "ZERO_GAP_TARGET", "BETA1_ACTIVE", "BETA2_ACTIVE",
               "ACTIVE_BLOCK", "BLOCK_PROPS", "BLOCK_RHO", "LATENT_TYPE",
               "LATENT_Q", "CENS_RATE", "U_MIN", "U_MAX")

.set_dgp <- .dgp_vars[nzchar(vapply(.dgp_vars, Sys.getenv, character(1)))]
if (length(.set_dgp) > 0L)
  stop(
    "The following DGP environment variable(s) are set: ",
    paste(.set_dgp, collapse = ", "), "\n",
    "The scenario file is authoritative for all DGP parameters. ",
    "Unset them to prevent silent override:\n",
    "  unset ", paste(.set_dgp, collapse = " "), "\n",
    "Authoritative scenario file: ",
    file.path(repo_root_val, "sim", "scenarios", paste0(scenario_name, ".R")),
    call. = FALSE
  )

# ---------------------------------------------------------------------------
# Load bhCRR and sim modules (always cluster mode -- never read BHCRR_LOAD)
# ---------------------------------------------------------------------------

source(file.path(repo_root_val, "sim", "R", "load.r"))
load_sim("cluster")

# Batch driver is in sim/local/batch.R, not autoloaded by load_sim
source(file.path(repo_root_val, "sim", "local", "batch.R"))

# ---------------------------------------------------------------------------
# Resolve and validate scenario file
# ---------------------------------------------------------------------------

scenario_file <- file.path(repo_root_val, "sim", "scenarios",
                            paste0(scenario_name, ".R"))
if (!file.exists(scenario_file))
  stop("Scenario file not found: ", scenario_file, "\n",
       "Create it with: new_scenario(\"", scenario_name,
       "\", dir = \"sim/scenarios\")", call. = FALSE)

# ---------------------------------------------------------------------------
# Startup banner (goes to the SLURM log)
# ---------------------------------------------------------------------------

message("=== bhCRR run_array.r ===")
message(sprintf("  repo_root    : %s", repo_root_val))
message(sprintf("  scenario     : %s", scenario_name))
message(sprintf("  scenario_file: %s", scenario_file))
message(sprintf("  scratch_dir  : %s", scratch_run_dir))
message(sprintf("  npredictors  : %s", paste(npredictors_grid, collapse = ", ")))
message(sprintf("  runs         : %d .. %d", run_start, run_end))

# ---------------------------------------------------------------------------
# Run the sweep (failure records, collision guard, manifest all in batch_sweep)
# ---------------------------------------------------------------------------

batch_sweep(
  scenario_file = scenario_file,
  npredictors   = npredictors_grid,
  run_start     = run_start,
  run_end       = run_end,
  runs_root     = scratch_run_dir,
  manifest_dir  = scratch_run_dir,
  manifest      = "check_only",
  load_mode     = "cluster",
  overwrite     = FALSE,
  dryrun        = FALSE
)
