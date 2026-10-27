#!/bin/bash
# =============================================================================
# sim/_legacy/chpc-timing/prep_timing.slurm  --  one-shot job: materialize the five bhCRR timing
# subsets (+ manifest) from the ClinTALL data. No model fitting here.
#
# Submitted directly with sbatch; resource flags (--time, --mem, etc.) are
# passed on the sbatch command line, matching timing.slurm's convention
# rather than being embedded as #SBATCH directives here.
#
# CLINTALL_DIR is exported to the repo-local inst/ClinTall BELOW, before
# sourcing timing_config.sh, so timing_config.sh's own
# ": ${CLINTALL_DIR:=$TIMING_SCRATCH/raw}" default does not win. The uploaded
# ClinTALL files live in the repo (inst/ClinTall/), not in scratch.
# =============================================================================
set -euo pipefail

echo "[$(date '+%F %T')] START host=$(hostname) job=${SLURM_JOB_ID:-NA}"

# ---- Resolve REPO_ROOT and force CLINTALL_DIR before config defaults apply --
# NOTE: under sbatch, SLURM copies this script into a per-job spool dir and
# runs it from there, so ${BASH_SOURCE[0]} does NOT point at the real repo
# checkout. Use $SLURM_SUBMIT_DIR (the directory `sbatch` was invoked from)
# instead, falling back to BASH_SOURCE only for non-batch invocations.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    export REPO_ROOT="$SLURM_SUBMIT_DIR"
else
    export REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
fi
HERE="$REPO_ROOT/sim/_legacy/chpc-timing"
export CLINTALL_DIR="$REPO_ROOT/inst/ClinTall"

source "$HERE/timing_config.sh"
echo "CLINTALL_DIR (resolved)    = $CLINTALL_DIR"
echo "TIMING_DATA_DIR (resolved) = $TIMING_DATA_DIR"

# This is a real prep run, not a smoke test -- make sure no stray env var
# (e.g. inherited via --export=ALL) turns on the smoke config.
unset TIMING_SMOKE || true

mkdir -p "$TIMING_DATA_DIR"

# ---- R environment ----------------------------------------------------------
module purge 2>/dev/null || true
module load "$R_MODULE"
cd "$REPO_ROOT"   # so .Rprofile activates renv for this project

if [[ "${USE_RENV:-0}" == "1" ]]; then
    Rscript -e 'renv::restore(prompt = FALSE)'
fi

# ---- Library / environment diagnostics --------------------------------------
# Deliberately run here, in the job, not on the login node: renv 1.1+ keeps
# the project library outside the repo and .Rprofile activation differs
# between login-node and submitted-job contexts, so a login-node library
# check is not just discouraged, it is actively misleading.
Rscript -e 'cat("R:", R.version.string, "\n")' \
        -e 'cat("libPaths:\n"); cat(.libPaths(), sep="\n"); cat("\n")' \
        -e 'cat("fastcmprsk:", system.file(package="fastcmprsk"), "\n")' \
        -e 'cat("tidyverse:",  system.file(package="tidyverse"),  "\n")'

# ---- Materialize the five timing subsets ------------------------------------
Rscript "$REPO_ROOT/sim/_legacy/chpc-timing/prep_timing_data.R"

# ---- Post-check: confirm the largest subset loads with expected dimensions --
# A real memory allocation -- deliberately done here, inside the job, not on
# the login node.
Rscript -e '
  f <- file.path(Sys.getenv("TIMING_DATA_DIR"), "timing_n1309-pall.rds")
  obj <- readRDS(f)
  cat("n1309-pall loaded dims:", nrow(obj$x), "x", ncol(obj$x), "\n")
  stopifnot(nrow(obj$x) == 1309, ncol(obj$x) == 24619)
  cat("n1309-pall dimension check: OK\n")
'

echo "[$(date '+%F %T')] END   host=$(hostname) job=${SLURM_JOB_ID:-NA}"
