#!/bin/bash
# =============================================================================
# sim/chpc/setup.slurm  --  one-shot prepare job: restore renv + compile bhCRR
#
# WHY THIS JOB EXISTS:
#   (1) CHPC Arbiter: login and interactive nodes are resource-governed (4 cores
#       / 8 GB cap). A full renv::restore() of 166 source packages and a
#       compile_dll() routinely exceed that budget and will be killed. Both MUST
#       run in a submitted batch job.
#   (2) Concurrent safety: array tasks share one REPO_ROOT on a shared
#       filesystem. If they each ran compile_dll(), they would race on
#       src/*.o and src/bhCRR.so, producing corrupt or mismatched objects.
#       load_sim("cluster") therefore uses compile=FALSE. The setup job compiles
#       once; the array tasks never compile.
#
# The renv library lives outside the repo under ~/.cache/R/renv/library/ and is
# visible from every CHPC resource, so one setup job covers all clusters.
#
# Submitted by run_sim.sh with --export=ALL. Resource requests (#SBATCH lines
# below are DEFAULTS; run_sim.sh overrides --time/--mem/--cpus-per-task on the
# sbatch command line from SB_SETUP_TIME/MEM/CPUS).
# =============================================================================
#SBATCH --job-name=bhcrr_setup
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=sim/chpc/logs/setup_%j.out

set -euo pipefail

echo "========================================"
echo "bhCRR setup job"
echo "[$(date '+%F %T')] host=$(hostname) job=${SLURM_JOB_ID:-NA}"
echo "  repo_root  : $REPO_ROOT"
echo "  R_MODULE   : $R_MODULE"
echo "  USE_RENV   : $USE_RENV"
echo "  CPUS       : ${SB_SETUP_CPUS:-4}"
echo "========================================"

# ---- R environment ----------------------------------------------------------
module purge 2>/dev/null || true
module load "$R_MODULE"

cd "$REPO_ROOT"   # activates renv via .Rprofile

# ---- Diagnostics ------------------------------------------------------------
Rscript - <<'REOF'
cat("R version:", R.version$version.string, "\n")
cat(".libPaths():\n"); cat(paste0("  ", .libPaths()), sep = "\n"); cat("\n")
for (pkg in c("fastcmprsk", "pkgload", "pkgbuild")) {
  loc <- tryCatch(find.package(pkg), error = function(e) NA_character_)
  cat(sprintf("  %-15s : %s\n", pkg, if (is.na(loc)) "NOT FOUND" else loc))
}
REOF

# ---- renv restore -----------------------------------------------------------
use_renv="${USE_RENV:-auto}"

case "$use_renv" in
  0)
    echo "[setup] USE_RENV=0: skipping renv::restore(). Ensure library is current."
    ;;
  1)
    echo "[setup] USE_RENV=1: running renv::restore() unconditionally."
    Rscript - <<REOF
options(Ncpus = as.integer(Sys.getenv("SB_SETUP_CPUS", "${SB_SETUP_CPUS:-4}")))
renv::restore(prompt = FALSE)
REOF
    ;;
  auto|*)
    echo "[setup] USE_RENV=auto: checking synchronization status."
    Rscript - <<REOF
options(Ncpus = as.integer(Sys.getenv("SB_SETUP_CPUS", "${SB_SETUP_CPUS:-4}")))
st <- renv::status()
# renv::status() returns a list; \$synchronized may be absent in older renv
synchronized <- if (!is.null(st\$synchronized)) isTRUE(st\$synchronized) else FALSE
if (synchronized) {
  cat("[setup] renv is synchronized -- skipping restore.\n")
} else {
  cat("[setup] renv is NOT synchronized -- restoring.\n")
  renv::restore(prompt = FALSE)
}
REOF
    ;;
esac

# ---- Compile the package ----------------------------------------------------
echo "[setup] Compiling bhCRR (pkgbuild::compile_dll)..."
Rscript -e 'pkgbuild::compile_dll()'
echo "[setup] Compilation done."

# ---- Verification suite -----------------------------------------------------
echo "[setup] Running verification checks..."

# 1. Shared object exists
SO_PATH="$REPO_ROOT/src/bhCRR.so"
if [[ ! -f "$SO_PATH" ]]; then
  echo "ERROR [setup]: src/bhCRR.so not found after compile_dll(). Aborting." >&2
  exit 1
fi
echo "  [ok] src/bhCRR.so exists"

# 2. Key packages loadable
Rscript - <<'REOF' || { echo "ERROR [setup]: package loadability check failed." >&2; exit 1; }
for (pkg in c("pkgload", "pkgbuild", "fastcmprsk")) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  cat(sprintf("  [%s] requireNamespace('%s')\n", if (ok) "ok" else "FAIL", pkg))
  if (!ok) stop(pkg, " is not available.", call. = FALSE)
}
REOF

# 3. bhCRR loads via load_sim("cluster")
Rscript - <<'REOF' || { echo "ERROR [setup]: load_sim('cluster') failed." >&2; exit 1; }
source(file.path(Sys.getenv("REPO_ROOT"), "sim", "R", "load.r"))
load_sim("cluster")
cat("  [ok] load_sim('cluster') succeeded\n")
REOF

# 4. fit_ssl_psdh callable
Rscript - <<'REOF' || { echo "ERROR [setup]: fit_ssl_psdh not found after load_sim." >&2; exit 1; }
source(file.path(Sys.getenv("REPO_ROOT"), "sim", "R", "load.r"))
load_sim("cluster")
stopifnot(existsFunction("fit_ssl_psdh"))
cat("  [ok] fit_ssl_psdh is callable\n")
REOF

echo "========================================"
echo "[$(date '+%F %T')] setup job complete -- all checks passed."
echo "========================================"
