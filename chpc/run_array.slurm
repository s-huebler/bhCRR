#!/bin/bash
# =============================================================================
# chpc/run_array.slurm  --  one array task = a chunk of run indices.
#
# Submitted by run_sim.sh with --export=ALL, so every config.sh variable
# (SCENARIO, NOBS, NPREDICTORS, RUN_START, RUN_END, RUNS_PER_TASK,
# SCRATCH_RUN_DIR, REPO_ROOT, R_MODULE, USE_RENV, ZERO_GAP_TARGET,
# BETA1_ACTIVE, BETA2_ACTIVE)
# is already in the environment. sbatch resource flags (--array, --account,
# --time, ...) are set on the run_sim.sh command line, not here.
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
echo "task $task -> runs ${RUN_SUBSTART}..${RUN_SUBEND}  (nobs=$NOBS, npredictors=$NPREDICTORS)"

# ---- R environment ----------------------------------------------------------
module purge 2>/dev/null || true
module load "$R_MODULE"

mkdir -p "$SCRATCH_RUN_DIR"
cd "$REPO_ROOT"   # so .Rprofile activates renv for this project

# Restore the renv library the first time (subsequent tasks are near-instant).
if [[ "${USE_RENV:-1}" == "1" ]]; then
  Rscript -e 'renv::restore(prompt = FALSE)'
fi

# ---- Run the models ----------------------------------------------------------
Rscript "$REPO_ROOT/chpc/R/run_batch.R"

echo "[$(date '+%F %T')] task ${task} done."
