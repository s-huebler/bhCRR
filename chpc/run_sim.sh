#!/usr/bin/env bash
# =============================================================================
# chpc/run_sim.sh  --  submit the bhCRR simulation as a SLURM job array, then a
# dependent parse job that collects betas_pXX into Sim/chpc_results/.
#
# Usage (from the repo root on CHPC):
#   ./chpc/run_sim.sh                          # all defaults from config.sh
#   ./chpc/run_sim.sh run_start=10 run_end=40  # override the run range
#   ./chpc/run_sim.sh nobs=140                 # override sample size
#   ./chpc/run_sim.sh npredictors=25,50,206 run_end=100 runs_per_task=5
#
# Notes:
#   * Overrides are key=value with NO spaces around '=' (nobs=140, not nobs = 140).
#   * Keys are case-insensitive; they map to the UPPER_CASE vars in config.sh.
#   * Add  dryrun=1  to print the sbatch commands without submitting.
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---- 1. Parse key=value overrides into the environment (before config.sh) ----
DRYRUN=0
for arg in "$@"; do
  case "$arg" in
    *=*)
      key="${arg%%=*}"; val="${arg#*=}"
      key="$(printf '%s' "$key" | tr '[:lower:]' '[:upper:]' | tr -d '[:space:]')"
      if [[ "$key" == "DRYRUN" ]]; then DRYRUN="$val"; continue; fi
      export "$key=$val"
      ;;
    *)
      echo "run_sim.sh: ignoring unrecognized argument '$arg' (expected key=value)" >&2
      ;;
  esac
done

# ---- 2. Fill in anything not overridden from the single config file ----------
# shellcheck source=/dev/null
source "$HERE/config.sh"

# ---- 3. Sanity checks --------------------------------------------------------
die() { echo "run_sim.sh: $*" >&2; exit 1; }
check_allocation || exit 1     # from config.sh: CHPC_ACCOUNT / CHPC_PARTITION set
[[ -d "$REPO_ROOT" ]]              || die "REPO_ROOT does not exist: $REPO_ROOT"
[[ "$RUN_END" -ge "$RUN_START" ]]  || die "run_end ($RUN_END) must be >= run_start ($RUN_START)."
[[ "$RUNS_PER_TASK" -ge 1 ]]       || die "runs_per_task must be >= 1."

# Optional cross-cluster routing: each CHPC cluster has its own Slurm controller.
# Empty CHPC_CLUSTER = submit to the login cluster's scheduler.
CLUSTER_FLAG=()
[[ -n "${CHPC_CLUSTER:-}" ]] && CLUSTER_FLAG=(--clusters="$CHPC_CLUSTER")

# Per-nobs scratch PARENT dir. Each npredictors value then gets its own
# subdirectory beneath it, created by run_batch.R as it sweeps the grid:
#
#   $SCRATCH_BASE/runs_n200/p15/n200_p15_run1.Rdata
#   $SCRATCH_BASE/runs_n200/p221/n200_p221_run1.Rdata
#
# The block map is a function of npredictors, so runs at different p must never
# be pooled into one summary. Separate directories make that structural rather
# than a convention parse_batch.R has to enforce after the fact.
#
# parse_batch.R lists recursively, so this variable doubles as the parse scope:
# point it at the parent to parse every p present, or at a single p<P> subdir to
# parse just that one. Runs sharing the same n and p accumulate in their subdir,
# so run_start=10 run_end=40 still adds to an earlier 1..10 batch.
SCRATCH_RUN_DIR="$SCRATCH_BASE/runs_n${NOBS}"
export SCRATCH_RUN_DIR
mkdir -p "$REPO_ROOT/chpc/logs" "$REPO_ROOT/Sim/chpc_results"
if [[ "$DRYRUN" != "0" ]]; then
  echo "  [dryrun] scratch dir NOT created: $SCRATCH_RUN_DIR"
else
  mkdir -p "$SCRATCH_RUN_DIR"
fi

# ---- 4. Array geometry -------------------------------------------------------
total_runs=$(( RUN_END - RUN_START + 1 ))
ntasks=$(( (total_runs + RUNS_PER_TASK - 1) / RUNS_PER_TASK ))   # ceil
array_max=$(( ntasks - 1 ))

cat <<EOF
------------------------------------------------------------
bhCRR CHPC submission
  repo_root      : $REPO_ROOT
  scratch (runs) : $SCRATCH_RUN_DIR
  results        : $REPO_ROOT/Sim/chpc_results
  R module       : $R_MODULE   (renv restore: $USE_RENV)
  account/part.  : $CHPC_ACCOUNT / $CHPC_PARTITION  (cluster: ${CHPC_CLUSTER:-<login>})
  nobs           : $NOBS
  npredictors    : $NPREDICTORS
  runs           : $RUN_START .. $RUN_END  ($total_runs runs)
  runs/task      : $RUNS_PER_TASK  ->  array 0-$array_max ($ntasks tasks)
  -- scenario --
  beta1_active   : $BETA1_ACTIVE
  beta2_active   : $BETA2_ACTIVE
  active_block   : $ACTIVE_BLOCK
  block_props    : $BLOCK_PROPS
  block_rho      : $BLOCK_RHO
  latent_type    : $LATENT_TYPE
  latent_q       : $LATENT_Q
  cens_rate      : $CENS_RATE  u_min: $U_MIN  u_max: $U_MAX
------------------------------------------------------------
EOF

submit() {  # echo + run (or just echo under dryrun); returns job id on stdout
  if [[ "$DRYRUN" != "0" ]]; then
    echo "[dryrun] $*" >&2
    echo "DRYRUN_JOBID"
  else
    "$@"
  fi
}

# ---- 5. Submit the run array -------------------------------------------------
run_jid=$(submit sbatch --parsable \
  --account="$CHPC_ACCOUNT" --partition="$CHPC_PARTITION" "${CLUSTER_FLAG[@]}" \
  --time="$SB_TIME" --mem="$SB_MEM" --cpus-per-task="$SB_CPUS" \
  --array="0-${array_max}" \
  --job-name="bhcrr_run_n${NOBS}" \
  --output="$REPO_ROOT/chpc/logs/run_%A_%a.out" \
  --export=ALL \
  "$HERE/run_array.slurm")
# With --clusters, `sbatch --parsable` returns "jobid;cluster"; keep just the id
# so the dependency below is a clean numeric job ID.
run_jid="${run_jid%%;*}"
echo "Submitted run array : job $run_jid  (array 0-$array_max)"

# ---- 6. Submit the parse job, gated on the whole array succeeding -------------
parse_jid=$(submit sbatch --parsable \
  --account="$CHPC_ACCOUNT" --partition="$CHPC_PARTITION" "${CLUSTER_FLAG[@]}" \
  --time="$SB_PARSE_TIME" --mem="$SB_MEM" --cpus-per-task=1 \
  --dependency="afterok:${run_jid}" \
  --job-name="bhcrr_parse_n${NOBS}" \
  --output="$REPO_ROOT/chpc/logs/parse_%j.out" \
  --export=ALL \
  "$HERE/parse.slurm")
parse_jid="${parse_jid%%;*}"
echo "Submitted parse job : job $parse_jid  (runs after array $run_jid succeeds)"

echo
squeue_cmd="squeue -u \$USER"
[[ -n "${CHPC_CLUSTER:-}" ]] && squeue_cmd="squeue -M $CHPC_CLUSTER -u \$USER"
echo "Track with:  $squeue_cmd"
echo "Results will appear in: $REPO_ROOT/Sim/chpc_results/betas_p<P>_n${NOBS}.{csv,Rdata}"
