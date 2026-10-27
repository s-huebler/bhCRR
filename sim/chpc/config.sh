# =============================================================================
# sim/chpc/config.sh  --  the ONE file you source before submitting.
#
# What lives here: cluster identity (allocation, partition), paths, sweep scope
# (npredictors, run range), and renv policy. NOTHING about the data-generating
# process -- DGP parameters live solely in sim/scenarios/<slug>.R.
#
#   ssh <you>@notchpeak.chpc.utah.edu
#   cd ~/Documents/bhCRR           # your CHPC checkout of the repo
#   source sim/chpc/config.sh      # (optional; run_sim.sh sources it too)
#   ./sim/chpc/run_sim.sh scenario=<name> run_start=1 run_end=40
#
# Every value below is set only if not already in the environment -- so any
# key=value override parsed by run_sim.sh (or an env var you export by hand)
# WINS over the default here. So `./sim/chpc/run_sim.sh R_MODULE=R/4.6.0`
# overrides R_MODULE without touching this file.
#
# Allocation is prefilled for the qiaox / lonepeak setup (carried over from the
# ODSiData CHPC config). The only thing to confirm once is R_MODULE.
# =============================================================================

# ---- Paths ------------------------------------------------------------------
# Where the repo lives on CHPC (this file's great-grandparent dir: sim/chpc/ -> repo root).
: "${REPO_ROOT:=$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/../.." && pwd)}"

# Shared scratch base -- must be visible from every compute node so the job
# array and the parse job all see the same files. CHPC general scratch:
#   /scratch/general/vast/$USER   or   /scratch/general/nfs1/$USER
: "${SCRATCH_BASE:=/scratch/general/vast/$USER/bhCRR}"

# ---- R environment ----------------------------------------------------------
# `module load` name. renv.lock pins R 4.5.1, so match the minor version.
: "${R_MODULE:=R/4.5.2}"          # <-- confirm this module exists: `module spider R`

# BHCRR_LOAD=cluster tells load_sim() to use pkgload::load_all(compile=FALSE).
# "dev" would recompile on every task and concurrent array tasks would race on
# src/*.o and src/bhCRR.so. The setup job compiles once; tasks never compile.
: "${BHCRR_LOAD:=cluster}"

# USE_RENV -- three-state renv restore policy for the setup job:
#   auto  restore only if renv reports the project is NOT synchronized with the
#         lockfile (detected via renv::status()$synchronized).  This is the
#         right default: it fixes a missing or outdated library without forcing
#         a full reinstall every time.  It also resolves the 2026-10-27 state
#         where the library under ~/.cache was incomplete but USE_RENV=0 was set
#         deliberately to avoid a broken restore -- "auto" bridges that gap.
#   1     always run renv::restore(prompt=FALSE), even if in sync.
#   0     skip renv::restore() entirely. Use this only when you know the library
#         is current and want to save the restore overhead. Note: as of
#         2026-10-27 the CHPC renv library was incomplete due to a PPM binary
#         incompatibility (see RENV_CONFIG_PPM_ENABLED below). Setting 0 with
#         an incomplete library will cause load failures in the array tasks.
: "${USE_RENV:=auto}"

# CHPC's R is built against an external LAPACK, so it has no libRlapack.so.
# Posit Package Manager's Linux binaries are built against an R that does, and
# fail to load here ("libRlapack.so: cannot open shared object file"). Disabling
# the PPM binary endpoint forces source installs that link against this R.
# The renv cellar (renv/cellar/fastcmprsk_1.24.10.tar.gz) is vendored in the
# repo and reaches CHPC on `git pull`, so renv uses it instead of fetching.
: "${RENV_CONFIG_PPM_ENABLED:=FALSE}"

# ---- SLURM allocation (logic carried from the ODSiData CHPC config) ----------
# Set these to your CHPC allocation. Override at submit time with env vars, e.g.
#   CHPC_ACCOUNT=my-alloc CHPC_PARTITION=notchpeak ./sim/chpc/run_sim.sh scenario=x
export CHPC_ACCOUNT="${CHPC_ACCOUNT:-biostat-division-np}"          # sbatch -A / --account
# Default to a SHARED partition, not a whole-node one. These jobs (each array
# task = a few cores + tens of GB) do not need a full node; a whole-node
# partition makes the job wait for an entire free node, the main cause of long
# queue waits. lonepeak-shared (qiaox allocation) is shared, has generous
# walltime, and usually has free nodes. It lives on the lonepeak scheduler, so
# CHPC_CLUSTER below routes sbatch there.
#
# Common overrides (set on the run_sim.sh command line or as env vars):
#   whole node:          CHPC_CLUSTER=lonepeak  CHPC_PARTITION=lonepeak
#   back to kingspeak:   CHPC_CLUSTER=""        CHPC_PARTITION=kingspeak-shared
#   free short (<=8h):   CHPC_CLUSTER=notchpeak CHPC_ACCOUNT=notchpeak-shared-short \
#                        CHPC_PARTITION=notchpeak-shared-short
export CHPC_PARTITION="${CHPC_PARTITION:-biostat-division-np}"  # sbatch -p / --partition
# Each CHPC cluster runs its OWN Slurm controller, so a partition is only valid
# on its own cluster -- submitting a lonepeak partition from a kingspeak context
# fails with "invalid partition specified". Empty = your login cluster's
# scheduler; set this (paired with a matching account/partition) to reach another
# cluster. Note the `-` (not `:-`): an explicitly-set empty CHPC_CLUSTER="" is
# preserved (routes to your login cluster), while leaving it unset defaults to
# notchpeak.
export CHPC_CLUSTER="${CHPC_CLUSTER-notchpeak}"        # sbatch -M / --clusters ("" = login cluster)

# ---- SLURM resource requests -- setup job -----------------------------------
# Measured 2026-10-27 on notchpeak: restoring 46 missing packages from source
# on 4 cores with PPM disabled took 8.1 min, peaked at 2.0 GB; compile_dll()
# took 20 s; whole job was 8:41. A cold restore into an empty library would be
# several times longer but still well under 1 h. Keeping the request tight also
# queues faster on notchpeak-shared-short (QOS cap: 16 CPUs, 128 GB, 8 h MaxTime).
: "${SB_SETUP_TIME:=01:00:00}"     # measured worst-case ~9 min; 1 h is a generous buffer
: "${SB_SETUP_MEM:=8G}"            # measured peak 2.0 GB; 8 G gives 4x headroom
: "${SB_SETUP_CPUS:=4}"            # passed to options(Ncpus=) inside setup.slurm

# ---- SLURM resource requests -- array and parse jobs ------------------------
: "${SB_TIME:=00:15:00}"           # walltime per array task (rough guess for p=221 + autotune; trim after first timing run)
: "${SB_MEM:=4G}"                  # memory per array task
: "${SB_CPUS:=1}"                  # cpus-per-task
: "${SB_PARSE_TIME:=00:20:00}"     # walltime for the (single) parse job

# ---- Simulation sweep scope --------------------------------------------------
# SCENARIO names a file in sim/scenarios/. DGP parameters are defined there --
# not here. What goes here: the sweep dimensions (which p values, which run range)
# and task geometry for SLURM.
#
# There is no default for SCENARIO -- you must pass scenario=<name> on the
# run_sim.sh command line.
: "${SCENARIO:=}"

: "${NPREDICTORS:=25}"             # single value OR comma list, e.g. 25,50,221
: "${RUN_START:=1}"                # first run index (inclusive)
: "${RUN_END:=10}"                 # last  run index (inclusive)
: "${RUNS_PER_TASK:=5}"            # runs handled by each array task

# CHPC_ACCOUNT / CHPC_PARTITION / CHPC_CLUSTER are already exported inline above.
export REPO_ROOT SCRATCH_BASE R_MODULE BHCRR_LOAD USE_RENV RENV_CONFIG_PPM_ENABLED \
       SB_SETUP_TIME SB_SETUP_MEM SB_SETUP_CPUS \
       SB_TIME SB_MEM SB_CPUS SB_PARSE_TIME \
       SCENARIO NPREDICTORS RUN_START RUN_END RUNS_PER_TASK

# Fail fast if the allocation was left unset (mirrors ODSiData's check).
check_allocation() {
    if [[ -z "$CHPC_ACCOUNT" || "$CHPC_ACCOUNT" == "CHANGE_ME" ]]; then
        echo "ERROR: set CHPC_ACCOUNT in sim/chpc/config.sh (or export it) before submitting." >&2
        return 1
    fi
    if [[ -z "$CHPC_PARTITION" ]]; then
        echo "ERROR: set CHPC_PARTITION in sim/chpc/config.sh (or export it) before submitting." >&2
        return 1
    fi
}
