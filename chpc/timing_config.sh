# =============================================================================
# chpc/timing_config.sh  --  the ONE file you source before submitting timing jobs.
#
#   ssh <you>@notchpeak.chpc.utah.edu
#   cd ~/bhCRR
#   source chpc/timing_config.sh      # (optional; run_timing.sh sources it too)
#   ./chpc/run_timing.sh              # all seven configs at default resources
#   ./chpc/run_timing.sh tags=n200-p300,n200-p5000 init_mode=zero
#
# Every value is set only if not already in the environment -- so any key=value
# override parsed by run_timing.sh (or an env var you export by hand) WINS over
# the default here.
# =============================================================================

# ---- Paths ------------------------------------------------------------------
: "${REPO_ROOT:=$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/.." && pwd)}"

# Scratch base for timing runs. Separate from SCRATCH_BASE (simulation runs)
# to avoid accidental mixing.
: "${TIMING_SCRATCH:=/scratch/general/vast/$USER/bhCRR_timing}"

# Where prep_timing_data.R wrote the subsets (must exist before run_timing.sh).
: "${TIMING_DATA_DIR:=$TIMING_SCRATCH/data}"

# Where run_timing.R writes its result .rds files.
: "${TIMING_OUT_DIR:=$TIMING_SCRATCH/results}"

# Where dat_cr.rds and dat_design.rds live on CHPC (copy or symlink them here).
: "${CLINTALL_DIR:=$TIMING_SCRATCH/raw}"

# ---- R environment ----------------------------------------------------------
: "${R_MODULE:=R/4.5.2}"          # confirm with: module spider R
: "${USE_RENV:=0}"
# PPM binary endpoint disabled: CHPC's R has no libRlapack.so and PPM binaries
# built against standard R fail to load. Force source installs instead.
: "${RENV_CONFIG_PPM_ENABLED:=FALSE}"

# ---- SLURM allocation -------------------------------------------------------
export CHPC_ACCOUNT="${CHPC_ACCOUNT:-biostat-division-np}"
export CHPC_PARTITION="${CHPC_PARTITION:-biostat-division-np}"
# Each CHPC cluster has its own Slurm controller; empty = login cluster's scheduler.
export CHPC_CLUSTER="${CHPC_CLUSTER-notchpeak}"

# ---- Seed -------------------------------------------------------------------
: "${TIMING_SEED:=20260901}"

# ---- Model hyperparameters --------------------------------------------------
# These match run_batch.R defaults; override any of them on the run_timing.sh line.

: "${SS0:=0.001}"
: "${SS1:=0.4}"
: "${THETA_A:=1}"

# THETA_B default matches run_batch.R (theta_b=1), NOT fit_ssl_psdh's own default
# (theta_b=ncol(x)). At p≈24619 this choice can substantially change the EM
# iteration count and therefore runtime -- override to explore the tradeoff.
: "${THETA_B:=1}"

: "${INITIAL_SPARSITY:=0.05}"
: "${MAXIT:=50}"
: "${EPSILON:=1e-04}"

# INIT_MODE controls how fit_ssl_psdh is initialized.
#   "cv"   -> init=NULL  (5-fold CV LASSO over init_lam_path; what you'd actually run)
#   "zero" -> init=rep(0,p)  (skips CV entirely; isolates EM cost)
# Use "zero" as a fallback if the CV initialization is the walltime bottleneck at
# high p, or to attribute cost between initialization and the EM loop.
: "${INIT_MODE:=cv}"

: "${INNER_MAXIT_START:=1000}"

# ---- Export everything -------------------------------------------------------
export REPO_ROOT TIMING_SCRATCH TIMING_DATA_DIR TIMING_OUT_DIR CLINTALL_DIR \
       R_MODULE USE_RENV RENV_CONFIG_PPM_ENABLED \
       TIMING_SEED \
       SS0 SS1 THETA_A THETA_B INITIAL_SPARSITY MAXIT EPSILON \
       INIT_MODE INNER_MAXIT_START

# Fail fast if the allocation was left unset.
check_allocation() {
    if [[ -z "$CHPC_ACCOUNT" || "$CHPC_ACCOUNT" == "CHANGE_ME" ]]; then
        echo "ERROR: set CHPC_ACCOUNT in chpc/timing_config.sh (or export it) before submitting." >&2
        return 1
    fi
    if [[ -z "$CHPC_PARTITION" ]]; then
        echo "ERROR: set CHPC_PARTITION in chpc/timing_config.sh (or export it) before submitting." >&2
        return 1
    fi
}
