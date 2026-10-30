#!/bin/bash
# =============================================================================
# sim/timing/collect.slurm  --  collect all timing result files into one CSV
#
# Submitted by submit_timing.sh with --dependency=afterany:<array_job_id>
# and --export=ALL,TIMING_TAG=...,TIMING_OUT_DIR=...,GRID_FILE=...,
#                REPS=...,ARRAY_JOB_ID=...
#
# This job runs even if some array tasks failed (afterany dependency).
# It writes sacct output for the array job to a temp file, then calls
# collect_timing.r which reads it to fill in missing (p,rep) cells.
# =============================================================================
#SBATCH --job-name=bhcrr_collect
#SBATCH --output=sim/chpc/logs/collect_%j.out

set -euo pipefail

echo "========================================"
echo "bhCRR timing collect"
echo "[$(date '+%F %T')] host=$(hostname) job=${SLURM_JOB_ID:-NA}"
echo "  TIMING_TAG     : ${TIMING_TAG:-<unset>}"
echo "  TIMING_OUT_DIR : ${TIMING_OUT_DIR:-<unset>}"
echo "  ARRAY_JOB_ID   : ${ARRAY_JOB_ID:-<unset>}"
echo "========================================"

module purge 2>/dev/null || true
module load "$R_MODULE"

cd "$REPO_ROOT"   # activates renv via .Rprofile

# Dump sacct output for the array job so R code can read it
SACCT_FILE="${TIMING_OUT_DIR}/sacct_${ARRAY_JOB_ID}.txt"
if [[ -n "${ARRAY_JOB_ID:-}" ]]; then
  SACCT_CLUSTER_FLAG=()
  [[ -n "${CHPC_CLUSTER:-}" ]] && SACCT_CLUSTER_FLAG=(--clusters="$CHPC_CLUSTER")
  sacct "${SACCT_CLUSTER_FLAG[@]}" \
    -j "$ARRAY_JOB_ID" \
    --format="JobID,JobName,State,Elapsed,Timelimit,MaxRSS,ExitCode" \
    --parsable2 --noheader=false \
    > "$SACCT_FILE" 2>/dev/null || true
  echo "[collect] sacct written to: $SACCT_FILE"
else
  SACCT_FILE=""
  echo "[collect] ARRAY_JOB_ID not set; skipping sacct."
fi

export SACCT_FILE

Rscript sim/timing/collect_timing.r

echo "========================================"
echo "[$(date '+%F %T')] collect job complete."
echo "========================================"
