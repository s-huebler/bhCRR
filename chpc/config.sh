# =============================================================================
# chpc/config.sh  --  the ONE file you source before submitting.
#
#   ssh <you>@notchpeak.chpc.utah.edu
#   cd ~/bhCRR                       # your CHPC checkout of the repo
#   source chpc/config.sh            # (optional; run_sim.sh sources it too)
#   ./chpc/run_sim.sh run_start=10 run_end=40
#
# Every value below is set only if not already in the environment -- so any
# key=value override parsed by run_sim.sh (or an env var you export by hand) WINS
# over the default here. So `./chpc/run_sim.sh nobs=140` overrides NOBS, and
# `CHPC_PARTITION=notchpeak-shared-short ./chpc/run_sim.sh` overrides the
# partition, without you touching this file.
#
# Allocation is prefilled for the qiaox / lonepeak setup (carried over from the
# ODSiData CHPC config). The only thing to confirm once is R_MODULE.
# =============================================================================

# ---- Paths ------------------------------------------------------------------
# Where the repo lives on CHPC (this file's grandparent dir by default).
: "${REPO_ROOT:=$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/.." && pwd)}"

# Shared scratch base -- must be visible from every compute node so the job
# array and the parse job all see the same files. CHPC general scratch:
#   /scratch/general/vast/$USER   or   /scratch/general/nfs1/$USER
: "${SCRATCH_BASE:=/scratch/general/vast/$USER/bhCRR}"

# ---- R environment ----------------------------------------------------------
# `module load` name. renv.lock pins R 4.5.1, so match the minor version.
: "${R_MODULE:=R/4.5.2}"          # <-- confirm this module exists: `module spider R`
# 1 = run renv::restore() before the models (recommended, matches renv.lock).
: "${USE_RENV:=0}"
# CHPC's R is built against an external LAPACK, so it has no libRlapack.so.
# Posit Package Manager's Linux binaries are built against an R that does, and
# fail to load here ("libRlapack.so: cannot open shared object file"). Disabling
# the PPM binary endpoint forces source installs that link against this R.
: "${RENV_CONFIG_PPM_ENABLED:=FALSE}"

# ---- SLURM allocation (logic carried from the ODSiData CHPC config) ----------
# Set these to your CHPC allocation. Override at submit time with env vars, e.g.
#   CHPC_ACCOUNT=my-alloc CHPC_PARTITION=notchpeak ./chpc/run_sim.sh nobs=140
export CHPC_ACCOUNT="${CHPC_ACCOUNT:-biostat-division-np}"          # sbatch -A / --account
# Default to a SHARED partition, not a whole-node one. These jobs (each array
# task = a few cores + tens of GB) do not need a full node; a whole-node
# partition makes the job wait for an entire free node, the main cause of long
# queue waits. lonepeak-shared (qiaox allocation) is shared, has generous
# walltime, and usually has free nodes. It lives on the lonepeak scheduler, so
# CHPC_CLUSTER below routes sbatch there.
#
# Common overrides (set on the run_sim.sh command line or as env vars):
#   whole node:        CHPC_CLUSTER=lonepeak  CHPC_PARTITION=lonepeak
#   back to kingspeak: CHPC_CLUSTER=""        CHPC_PARTITION=kingspeak-shared
#   free short (<=8h): CHPC_CLUSTER=notchpeak CHPC_ACCOUNT=notchpeak-shared-short \
#                      CHPC_PARTITION=notchpeak-shared-short
export CHPC_PARTITION="${CHPC_PARTITION:-biostat-division-np}"  # sbatch -p / --partition
# Each CHPC cluster runs its OWN Slurm controller, so a partition is only valid
# on its own cluster -- submitting a lonepeak partition from a kingspeak context
# fails with "invalid partition specified". Empty = your login cluster's
# scheduler; set this (paired with a matching account/partition) to reach another
# cluster. run_sim.sh turns this into `sbatch --clusters=$CHPC_CLUSTER`.
# Note the `-` (not `:-`): an explicitly-set empty CHPC_CLUSTER="" is preserved
# (routes to your login cluster), while leaving it unset defaults to lonepeak.
export CHPC_CLUSTER="${CHPC_CLUSTER-notchpeak}"        # sbatch -M / --clusters ("" = login cluster)

# ---- SLURM resource requests (per job; tune freely) -------------------------
: "${SB_TIME:=00:15:00}"           # walltime per array task (rough guess for p=221 + autotune; trim after first timing run)
: "${SB_MEM:=4G}"                # memory per array task
: "${SB_CPUS:=1}"                 # cpus-per-task
: "${SB_PARSE_TIME:=00:20:00}"    # walltime for the (single) parse job

# ---- Simulation defaults (override any of these on the run_sim.sh line) ------
#
# SCENARIO names a user-maintained set of data-generating settings. All DGP
# parameters below (NOBS through U_MAX) are considered part of the scenario
# definition -- changing any of them requires a new scenario name. n (NOBS) is
# part of the scenario and therefore drops out of individual file naming.
#
# NOT part of the scenario: NPREDICTORS, RUN_START, RUN_END, RUNS_PER_TASK
# (these control how you sweep or extend a scenario, not what the scenario is),
# and all SLURM / account knobs.
#
# There is no default -- you must pass scenario=<name> on the run_sim.sh line.
: "${SCENARIO:=}"

: "${NOBS:=200}"                                   # sample size
: "${NPREDICTORS:=25}"                            # single value OR comma list, e.g. 25,50,221
: "${RUN_START:=1}"                                # first run index (inclusive)
: "${RUN_END:=10}"                                 # last  run index (inclusive)
: "${RUNS_PER_TASK:=5}"                            # runs handled by each array task
: "${ZERO_GAP_TARGET:=0.1}"                        # clinically-relevant min treatment effect
: "${BETA1_ACTIVE:=0.40,-0.50,0.60,0.75,-0.80}"    # active cause-1 coefficients
: "${BETA2_ACTIVE:=neg_beta1}"                      # cause-2 coefficients; "neg_beta1" = -BETA1_ACTIVE
#
# Named-vector variables below use "name=val,name=val" format. run_sim.sh parses
# overrides with key="${arg%%=*}" / val="${arg#*=}", which splits only on the FIRST
# '=', so extra '=' inside the value are safe. Example:
#   ./chpc/run_sim.sh block_rho=indep=0.0,weak=0.5,strong=0.8
: "${ACTIVE_BLOCK:=strong,strong,weak,weak,indep}"   # block labels aligned with BETA1_ACTIVE
: "${BLOCK_PROPS:=indep=0.50,weak=0.25,strong=0.25}" # block proportions (name=val,...)
: "${BLOCK_RHO:=indep=0.00,weak=0.30,strong=0.60}"   # within-block correlations (name=val,...)
: "${LATENT_TYPE:=gaussian}"                          # latent variable distribution
: "${LATENT_Q:=indep=0.5,weak=0.4,strong=0.5}"       # latent quantile per block (name=val,...)
: "${CENS_RATE:=0.05}"                                # Exp rate for random censoring (0 = none)
: "${U_MIN:=100}"                                     # lower bound for observation window
: "${U_MAX:=100}"                                     # upper bound for observation window

# CHPC_ACCOUNT / CHPC_PARTITION / CHPC_CLUSTER are already exported inline above.
export REPO_ROOT SCRATCH_BASE R_MODULE USE_RENV RENV_CONFIG_PPM_ENABLED \
       SB_TIME SB_MEM SB_CPUS SB_PARSE_TIME \
       SCENARIO \
       NOBS NPREDICTORS RUN_START RUN_END RUNS_PER_TASK ZERO_GAP_TARGET \
       BETA1_ACTIVE BETA2_ACTIVE \
       ACTIVE_BLOCK BLOCK_PROPS BLOCK_RHO LATENT_TYPE LATENT_Q \
       CENS_RATE U_MIN U_MAX

# Fail fast if the allocation was left unset (mirrors ODSiData's check).
check_allocation() {
    if [[ -z "$CHPC_ACCOUNT" || "$CHPC_ACCOUNT" == "CHANGE_ME" ]]; then
        echo "ERROR: set CHPC_ACCOUNT in chpc/config.sh (or export it) before submitting." >&2
        return 1
    fi
    if [[ -z "$CHPC_PARTITION" ]]; then
        echo "ERROR: set CHPC_PARTITION in chpc/config.sh (or export it) before submitting." >&2
        return 1
    fi
}
