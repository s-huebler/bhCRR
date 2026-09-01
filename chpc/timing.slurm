#!/bin/bash
# =============================================================================
# chpc/timing.slurm  --  one task = one bhCRR timing run on a real data config.
#
# Submitted by run_timing.sh with --export=ALL, so every timing_config.sh
# variable is already in the environment. Resource flags (--time, --mem, etc.)
# are passed on the sbatch command line, not here.
#
# TIMING_TAG is exported by run_timing.sh before submission. If not set (e.g.,
# when submitting this script directly as an array job), it is derived from
# SLURM_ARRAY_TASK_ID using the canonical index-to-tag mapping below:
#
#   0 -> n1309-p300     3 -> n1309-p10000
#   1 -> n1309-p1000    4 -> n1309-pall
#   2 -> n1309-p5000
# =============================================================================
set -euo pipefail

# ---- Resolve the config tag -------------------------------------------------
_TIMING_TAGS=(
    "n1309-p300"
    "n1309-p1000"
    "n1309-p5000"
    "n1309-p10000"
    "n1309-pall"
)

if [[ -z "${TIMING_TAG:-}" ]]; then
    _task="${SLURM_ARRAY_TASK_ID:-0}"
    if (( _task >= ${#_TIMING_TAGS[@]} )); then
        echo "ERROR: SLURM_ARRAY_TASK_ID=$_task is out of range (max $((${#_TIMING_TAGS[@]}-1)))." >&2
        exit 1
    fi
    export TIMING_TAG="${_TIMING_TAGS[$_task]}"
fi

echo "[$(date '+%F %T')] START host=$(hostname) job=${SLURM_JOB_ID:-NA} task=${SLURM_ARRAY_TASK_ID:-NA} tag=${TIMING_TAG}"

# ---- R environment ----------------------------------------------------------
module purge 2>/dev/null || true
module load "$R_MODULE"

mkdir -p "$TIMING_OUT_DIR"
cd "$REPO_ROOT"   # so .Rprofile activates renv for this project

if [[ "${USE_RENV:-0}" == "1" ]]; then
    Rscript -e 'renv::restore(prompt = FALSE)'
fi

# ---- Run the model ----------------------------------------------------------
Rscript "$REPO_ROOT/chpc/R/run_timing.R"

echo "[$(date '+%F %T')] END   host=$(hostname) job=${SLURM_JOB_ID:-NA} task=${SLURM_ARRAY_TASK_ID:-NA} tag=${TIMING_TAG}"
