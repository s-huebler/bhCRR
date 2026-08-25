#!/usr/bin/env bash
# =============================================================================
# chpc/run_sim.sh  --  submit the bhCRR simulation as a SLURM job array, then a
# dependent parse job that collects betas_pXX into Sim/chpc_results/.
#
# Usage (from the repo root on CHPC):
#   ./chpc/run_sim.sh scenario=n200_baseline              # required; all other defaults from config.sh
#   ./chpc/run_sim.sh scenario=n200_baseline run_start=10 run_end=40
#   ./chpc/run_sim.sh scenario=n200_baseline npredictors=25,50,206 run_end=100 runs_per_task=5
#
# Notes:
#   * scenario=<name> is REQUIRED. The name must match ^[A-Za-z0-9][A-Za-z0-9.-]*$
#     (no underscores -- they delimit fields in p<P>_<scenario>_run<N> filenames).
#   * Overrides are key=value with NO spaces around '=' (nobs=140, not nobs = 140).
#   * Keys are case-insensitive; they map to the UPPER_CASE vars in config.sh.
#   * Add  dryrun=1    to print the sbatch commands without submitting.
#   * Add  overwrite=1 to re-submit run indices that already have .Rdata output.
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---- 1. Parse key=value overrides into the environment (before config.sh) ----
DRYRUN=0
OVERWRITE=0
for arg in "$@"; do
  case "$arg" in
    *=*)
      key="${arg%%=*}"; val="${arg#*=}"
      key="$(printf '%s' "$key" | tr '[:lower:]' '[:upper:]' | tr -d '[:space:]')"
      if [[ "$key" == "DRYRUN" ]];    then DRYRUN="$val";    continue; fi
      if [[ "$key" == "OVERWRITE" ]]; then OVERWRITE="$val"; continue; fi
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

# SCENARIO is required and must be filename-safe with no underscores.
# Underscores are reserved as field delimiters in the p<P>_<scenario>_run<N>
# filename format; allowing them in the scenario name would make that regex
# ambiguous when parsing results.
[[ -n "${SCENARIO:-}" ]] || \
  die "scenario is required -- pass scenario=<name> on the command line."
[[ "$SCENARIO" =~ ^[A-Za-z0-9][A-Za-z0-9.-]*$ ]] || \
  die "SCENARIO '$SCENARIO' is invalid. Use only letters, digits, dots, and hyphens (no underscores -- they are delimiters in p<P>_<scenario>_run<N> filenames)."

# Optional cross-cluster routing: each CHPC cluster has its own Slurm controller.
# Empty CHPC_CLUSTER = submit to the login cluster's scheduler.
CLUSTER_FLAG=()
[[ -n "${CHPC_CLUSTER:-}" ]] && CLUSTER_FLAG=(--clusters="$CHPC_CLUSTER")

# Per-scenario scratch PARENT dir. Each npredictors value then gets its own
# subdirectory beneath it, created by run_batch.R as it sweeps the grid:
#
#   $SCRATCH_BASE/<scenario>/p15/<scenario>_p15_run1.Rdata
#   $SCRATCH_BASE/<scenario>/p221/<scenario>_p221_run1.Rdata
#
# The block map is a function of npredictors, so runs at different p must never
# be pooled into one summary. Separate directories make that structural rather
# than a convention parse_batch.R has to enforce after the fact.
#
# parse_batch.R lists recursively, so this variable doubles as the parse scope:
# point it at the scenario root to parse every p present, or at a single p<P>
# subdir to parse just that one. Runs sharing the same scenario and p accumulate
# in their subdir, so run_start=10 run_end=40 still adds to an earlier 1..10
# batch.
SCRATCH_RUN_DIR="$SCRATCH_BASE/${SCENARIO}"
export SCRATCH_RUN_DIR

mkdir -p "$REPO_ROOT/chpc/logs" "$REPO_ROOT/Sim/chpc_results/$SCENARIO"
if [[ "$DRYRUN" != "0" ]]; then
  echo "  [dryrun] scratch dir NOT created: $SCRATCH_RUN_DIR"
else
  mkdir -p "$SCRATCH_RUN_DIR"
fi

# ---- 3b. Scenario manifest ---------------------------------------------------
# Signature = the 12 DGP keys (sorted alphabetically). Deliberately excluded:
# NPREDICTORS, RUN_START, RUN_END, RUNS_PER_TASK, and all SLURM/account knobs.
_sig_keys=(ACTIVE_BLOCK BETA1_ACTIVE BETA2_ACTIVE BLOCK_PROPS BLOCK_RHO CENS_RATE
           LATENT_Q LATENT_TYPE NOBS U_MAX U_MIN ZERO_GAP_TARGET)

scenario_env="$SCRATCH_RUN_DIR/scenario.env"

if [[ -f "$scenario_env" ]]; then
  # Compare the recorded signature line by line against current values.
  mismatch=0
  for k in "${_sig_keys[@]}"; do
    recorded=$(grep "^${k}=" "$scenario_env" | cut -d= -f2-)
    current="${!k}"
    if [[ "$recorded" != "$current" ]]; then
      printf 'run_sim.sh: scenario mismatch -- %s:\n  recorded : %s\n  current  : %s\n' \
        "$k" "$recorded" "$current" >&2
      mismatch=1
    fi
  done
  [[ "$mismatch" -eq 0 ]] || \
    die "Scenario '$SCENARIO' has recorded settings that differ from what you passed (see above). Fix the mismatch or choose a new scenario name with scenario=<new>."
elif [[ "$DRYRUN" == "0" ]]; then
  {
    printf '# scenario:  %s\n' "$SCENARIO"
    printf '# created:   %s\n' "$(date '+%Y-%m-%dT%H:%M:%S')"
    printf '# git:       %s\n' "$(git -C "$REPO_ROOT" rev-parse --short HEAD 2>/dev/null || echo 'unknown')"
    for k in "${_sig_keys[@]}"; do printf '%s=%s\n' "$k" "${!k}"; done
  } > "$scenario_env"
  echo "  wrote scenario manifest: $scenario_env"
else
  echo "  [dryrun] scenario.env would be created at: $scenario_env"
fi

# ---- 3c. Collision guard -----------------------------------------------------
# Refuse to submit if any requested run indices already have .Rdata output under
# $SCRATCH_RUN_DIR/p<P>/, unless overwrite=1 was passed.
if [[ "$OVERWRITE" != "1" ]]; then
  any_collision=0
  for P in $(echo "$NPREDICTORS" | tr ',' ' '); do
    pdir="$SCRATCH_RUN_DIR/p${P}"
    colliding=()
    max_collide=0
    for i in $(seq "$RUN_START" "$RUN_END"); do
      files=( "${pdir}/"*"_run${i}.Rdata" )
      if [[ -e "${files[0]}" ]]; then
        colliding+=("$i")
        [[ "$i" -gt "$max_collide" ]] && max_collide=$i
      fi
    done
    if [[ ${#colliding[@]} -gt 0 ]]; then
      printf 'run_sim.sh: collision for p=%s -- existing runs: %s\n' \
        "$P" "${colliding[*]}" >&2
      printf '  Pass overwrite=1 to overwrite, or run_start=%d to append.\n' \
        $(( max_collide + 1 )) >&2
      any_collision=1
    fi
  done
  [[ "$any_collision" -eq 0 ]] || \
    die "Refusing to submit due to run collisions (see above). Pass overwrite=1 or adjust run_start."
fi

# ---- 4. Array geometry -------------------------------------------------------
total_runs=$(( RUN_END - RUN_START + 1 ))
ntasks=$(( (total_runs + RUNS_PER_TASK - 1) / RUNS_PER_TASK ))   # ceil
array_max=$(( ntasks - 1 ))

cat <<EOF
------------------------------------------------------------
bhCRR CHPC submission
  scenario       : $SCENARIO
  repo_root      : $REPO_ROOT
  scratch (runs) : $SCRATCH_RUN_DIR
  results        : $REPO_ROOT/Sim/chpc_results/$SCENARIO
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
  --job-name="bhcrr_run_${SCENARIO}" \
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
  --job-name="bhcrr_parse_${SCENARIO}" \
  --output="$REPO_ROOT/chpc/logs/parse_%j.out" \
  --export=ALL \
  "$HERE/parse.slurm")
parse_jid="${parse_jid%%;*}"
echo "Submitted parse job : job $parse_jid  (runs after array $run_jid succeeds)"

echo
squeue_cmd="squeue -u \$USER"
[[ -n "${CHPC_CLUSTER:-}" ]] && squeue_cmd="squeue -M $CHPC_CLUSTER -u \$USER"
echo "Track with:  $squeue_cmd"
echo "Results will appear in: $REPO_ROOT/Sim/chpc_results/$SCENARIO/betas_p<P>_${SCENARIO}.{csv,Rdata}"
