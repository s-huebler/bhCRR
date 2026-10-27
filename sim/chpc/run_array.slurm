#!/bin/bash
# =============================================================================
# sim/chpc/run_array.slurm  --  one array task = a chunk of run indices.
#
# Submitted by run_sim.sh with --export=ALL, so every config.sh variable
# (SCENARIO, NPREDICTORS, RUN_START, RUN_END, RUNS_PER_TASK,
# SCRATCH_RUN_DIR, REPO_ROOT, R_MODULE, BHCRR_LOAD)
# is already in the environment. DGP parameters are NOT exported -- they live
# solely in sim/scenarios/$SCENARIO.R (run_array.r errors if any DGP env var
# is set). sbatch resource flags (--array, --account, --time, ...) are set on
# the run_sim.sh command line, not here.
#
# renv::restore() and compile_dll() are NOT run here -- the setup job (which
# this array depends on via afterok) already did both. Running them again here
# would (a) re-race on src/*.o across concurrent tasks and (b) waste time.
# =============================================================================
set -euo pipefail

echo "[$(date '+%F %T')] host=$(hostname) job=${SLURM_JOB_ID:-NA} task=${SLURM_ARRAY_TASK_ID:-NA}"

# ---- This task's slice of the run range -------------------------------------
task="${SLURM_ARRAY_TASK_ID:-0}"
sub_start=$(( RUN_START + task * RUNS_PER_TASK ))
sub_end=$(( sub_start + RUNS_PER_TASK - 1 ))
if (( sub_end > RUN_END )); then sub_end="$RUN_END"; fi
if (( sub_start > RUN_END )); then
  echo "task $task has no runs (sub_start=$sub_start > RUN_END=$RUN_END); exiting cleanly."
  exit 0
fi
export RUN_SUBSTART="$sub_start" RUN_SUBEND="$sub_end"
echo "task $task -> runs ${RUN_SUBSTART}..${RUN_SUBEND}  (npredictors=$NPREDICTORS)"

# ---- R environment ----------------------------------------------------------
module purge 2>/dev/null || true
module load "$R_MODULE"

mkdir -p "$SCRATCH_RUN_DIR"
cd "$REPO_ROOT"   # so .Rprofile activates renv for this project

# ---- Run the models ----------------------------------------------------------
Rscript "$REPO_ROOT/sim/chpc/R/run_array.r"

echo "[$(date '+%F %T')] task ${task} done."
