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
#
# ORDERING NOTE: INIT_MODE must appear before TIMING_OUT_DIR in this file,
# because TIMING_OUT_DIR's default is derived from INIT_MODE via the
# `: "${VAR:=default}"` idiom, which evaluates the default at the line where
# it appears. Do not move TIMING_OUT_DIR above INIT_MODE.
# =============================================================================

# ---- Paths ------------------------------------------------------------------
: "${REPO_ROOT:=$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/.." && pwd)}"

# Scratch base for timing runs. Separate from SCRATCH_BASE (simulation runs)
# to avoid accidental mixing.
: "${TIMING_SCRATCH:=/scratch/general/vast/$USER/bhCRR_timing}"

# Where prep_timing_data.R wrote the subsets (must exist before run_timing.sh).
: "${TIMING_DATA_DIR:=$TIMING_SCRATCH/data}"

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

# INIT_MODE controls how fit_ssl_psdh is initialized. Three options:
#
#   "bic"  -- (default) Fit one fastCrrp LASSO path over init_lam_path, select
#              the lambda with the lowest BIC, and pass those coefficients as
#              `init`. No fold splitting; stable even with few events. Uses the
#              same lambda grid as "cv" but fits the data once, not k+1 times.
#
#   "cv"   -- pass init=NULL so fit_ssl_psdh runs its own internal 5-fold
#              cv_fastCrrp_cpp search. Realistic (this is what a real analysis
#              would run), but with <4 cause-1 events per fold in n200 configs
#              the CV can be noisy or degenerate. Fit time includes the CV.
#
#   "zero" -- pass init=rep(0,p), skipping any initial model entirely. Fastest;
#              isolates pure EM cost and is a useful lower bound on wall time.
#
# TIMING_OUT_DIR below is derived from this variable -- INIT_MODE MUST be
# defined here before TIMING_OUT_DIR is evaluated.
: "${INIT_MODE:=bic}"

# ---- Output directory (derived from TIMING_SCRATCH and INIT_MODE) -----------
# Each init mode writes to its own subdirectory so runs never overwrite each
# other. An explicit timing_out_dir= override on the run_timing.sh command line
# still wins (exported before this file is sourced, so the := skips the default).
#
# Stale-source caveat: if you `source timing_config.sh` twice in the SAME shell
# with different INIT_MODE values, the second source does NOT update
# TIMING_OUT_DIR because the variable is already exported from the first source.
# This is harmless for our usage (every Bash call is a fresh shell) but worth
# knowing if you source this file interactively.
: "${TIMING_OUT_DIR:=$TIMING_SCRATCH/results_$INIT_MODE}"

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
