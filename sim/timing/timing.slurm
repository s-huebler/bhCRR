#!/bin/bash
# =============================================================================
# sim/timing/timing.slurm  --  one array task = one fit_ssl_psdh call
#
# Resource flags (--time, --mem, --array, etc.) are NOT set here.
# Pass them on the sbatch line from submit_timing.sh.
#
# Submitted by submit_timing.sh with --export=ALL.
# Required exports: TIMING_INPUT, TIMING_OUT_DIR, TIMING_TAG, REPO_ROOT,
#                   R_MODULE, BHCRR_LOAD, (optionally) REPS, INIT_METHOD,
#                   SS0, SS1.
# =============================================================================
#SBATCH --job-name=bhcrr_timing
#SBATCH --output=sim/chpc/logs/timing_%A_%a.out

set -euo pipefail

# Single-thread BLAS/LAPACK so timing measures wall time, not core count.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "========================================"
echo "bhCRR timing task"
echo "[$(date '+%F %T')] host=$(hostname) job=${SLURM_JOB_ID:-NA} task=${SLURM_ARRAY_TASK_ID:-NA}"
echo "  TIMING_TAG   : ${TIMING_TAG:-<unset>}"
echo "  TIMING_INPUT : ${TIMING_INPUT:-<unset>}"
echo "  TIMING_OUT   : ${TIMING_OUT_DIR:-<unset>}"
echo "  REPO_ROOT    : ${REPO_ROOT:-<unset>}"
echo "  BHCRR_LOAD   : ${BHCRR_LOAD:-cluster}"
echo "  INIT_METHOD  : ${INIT_METHOD:-LASSO_bic}"
echo "========================================"

module purge 2>/dev/null || true
module load "$R_MODULE"

cd "$REPO_ROOT"   # activates renv via .Rprofile

Rscript sim/timing/fit_one_timing.r
