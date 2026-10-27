#!/bin/bash
# =============================================================================
# sim/chpc/parse.slurm  --  runs once, after the whole run array succeeds
# (--dependency=afterok). Collects every .Rdata in $SCRATCH_RUN_DIR, groups
# by npredictors, and writes summary CSVs to sim/results/$SCENARIO/parsed/.
#
# Inherits the environment from run_sim.sh via --export=ALL
# (SCENARIO, SCRATCH_RUN_DIR, REPO_ROOT, R_MODULE, BHCRR_LOAD).
# DGP parameters are NOT exported -- they live solely in the scenario file.
#
# renv::restore() and compile_dll() are NOT run here -- the setup job already
# did both.
# =============================================================================
set -euo pipefail

echo "[$(date '+%F %T')] parse on host=$(hostname) job=${SLURM_JOB_ID:-NA}"
echo "scratch : $SCRATCH_RUN_DIR"
echo "results : $REPO_ROOT/sim/results/$SCENARIO/parsed"

module purge 2>/dev/null || true
module load "$R_MODULE"

cd "$REPO_ROOT"

Rscript "$REPO_ROOT/sim/chpc/R/parse.r"

echo "[$(date '+%F %T')] parse done."
