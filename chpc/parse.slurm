#!/bin/bash
# =============================================================================
# chpc/parse.slurm  --  runs once, after the whole run array succeeds
# (--dependency=afterok). Collects every RData in $SCRATCH_RUN_DIR, groups by
# npredictors, and writes betas_p<P>_n<NOBS>.{csv,Rdata} to Sim/chpc_results/.
#
# Inherits the environment from run_sim.sh via --export=ALL
# (SCENARIO, SCRATCH_RUN_DIR, REPO_ROOT, R_MODULE, USE_RENV, and all DGP vars).
# =============================================================================
set -euo pipefail

echo "[$(date '+%F %T')] parse on host=$(hostname) job=${SLURM_JOB_ID:-NA}"
echo "scratch : $SCRATCH_RUN_DIR"
echo "results : $REPO_ROOT/Sim/chpc_results"

module purge 2>/dev/null || true
module load "$R_MODULE"

cd "$REPO_ROOT"
if [[ "${USE_RENV:-1}" == "1" ]]; then
  Rscript -e 'renv::restore(prompt = FALSE)'
fi

Rscript "$REPO_ROOT/chpc/R/parse_batch.R"

echo "[$(date '+%F %T')] parse done."
