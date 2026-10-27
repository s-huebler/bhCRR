#!/usr/bin/env bash
# =============================================================================
# sim/chpc/run_sim.sh  --  submit the bhCRR simulation as a SLURM job chain:
#                          setup → run array → parse
#
# Usage (from the repo root on CHPC):
#   ./sim/chpc/run_sim.sh scenario=n200.baseline              # required; defaults from config.sh
#   ./sim/chpc/run_sim.sh scenario=n200.baseline run_start=10 run_end=40
#   ./sim/chpc/run_sim.sh scenario=n200.baseline npredictors=25,50,206 run_end=100
#   ./sim/chpc/run_sim.sh scenario=n200.baseline skip_setup=1  # skip setup job (SO must exist)
#
# Notes:
#   * scenario=<name> is REQUIRED. The name must match ^[A-Za-z0-9][A-Za-z0-9.-]*$
#     (no underscores -- they delimit fields in p<P>_<scenario>_run<N> filenames).
#   * A sim/scenarios/<name>.R file must exist; this is the authority for all DGP
#     parameters. Create one with: source sim/R/load.r; load_sim(); new_scenario("<name>")
#   * Overrides are key=value with NO spaces around '=' (run_end=40, not run_end = 40).
#   * Keys are case-insensitive; they map to the UPPER_CASE vars in config.sh.
#   * Add  dryrun=1     to print the sbatch commands without submitting.
#   * Add  overwrite=1  to re-submit run indices that already have .Rdata output.
#   * Add  skip_setup=1 to skip the setup job and submit the array directly.
#                        USE ONLY when src/bhCRR.so is already compiled and the
#                        renv library is known to be current.
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---- 1. Parse key=value overrides into the environment (before config.sh) ----
DRYRUN=0
OVERWRITE=0
SKIP_SETUP=0
for arg in "$@"; do
  case "$arg" in
    *=*)
      key="${arg%%=*}"; val="${arg#*=}"
      key="$(printf '%s' "$key" | tr '[:lower:]' '[:upper:]' | tr -d '[:space:]')"
      if [[ "$key" == "DRYRUN" ]];      then DRYRUN="$val";      continue; fi
      if [[ "$key" == "OVERWRITE" ]];   then OVERWRITE="$val";   continue; fi
      if [[ "$key" == "SKIP_SETUP" ]];  then SKIP_SETUP="$val";  continue; fi
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

# Scenario file must exist -- it is the authority for all DGP parameters.
scenario_file="$REPO_ROOT/sim/scenarios/${SCENARIO}.R"
[[ -f "$scenario_file" ]] || \
  die "Scenario file not found: $scenario_file
Create it with:
  source sim/R/load.r; load_sim(); new_scenario(\"${SCENARIO}\", dir = \"sim/scenarios\")"

# Optional cross-cluster routing: each CHPC cluster has its own Slurm controller.
# Empty CHPC_CLUSTER = submit to the login cluster's scheduler.
CLUSTER_FLAG=()
[[ -n "${CHPC_CLUSTER:-}" ]] && CLUSTER_FLAG=(--clusters="$CHPC_CLUSTER")

# Per-scenario scratch dir. batch_sweep() creates p<P>/ subdirs beneath this.
#   $SCRATCH_BASE/<scenario>/p25/p25_<scenario>_run1.Rdata
SCRATCH_RUN_DIR="$SCRATCH_BASE/${SCENARIO}"
export SCRATCH_RUN_DIR

mkdir -p "$REPO_ROOT/sim/chpc/logs"
if [[ "$DRYRUN" != "0" ]]; then
  echo "  [dryrun] directories NOT created"
else
  mkdir -p "$SCRATCH_RUN_DIR" "$REPO_ROOT/sim/results/$SCENARIO/parsed"
fi

# ---- 3b. Collision guard -----------------------------------------------------
# Refuse to submit if any requested run indices already have .Rdata output under
# $SCRATCH_RUN_DIR/p<P>/, unless overwrite=1 was passed.
# This check runs before any job is queued, so it fails fast without wasting
# queue time. batch_sweep() also checks on the compute node; keeping both is
# deliberate (the bash check is cheap and doesn't need the renv library).
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

# Resolve the USE_RENV display string for the banner
use_renv_display="${USE_RENV:-auto}"

cat <<EOF
------------------------------------------------------------
bhCRR CHPC submission
  scenario       : $SCENARIO
  scenario file  : $scenario_file
  repo_root      : $REPO_ROOT
  scratch (runs) : $SCRATCH_RUN_DIR
  results        : $REPO_ROOT/sim/results/$SCENARIO/parsed
  R module       : $R_MODULE
  USE_RENV       : $use_renv_display
  BHCRR_LOAD     : ${BHCRR_LOAD:-cluster}
  account/part.  : $CHPC_ACCOUNT / $CHPC_PARTITION  (cluster: ${CHPC_CLUSTER:-<login>})
  npredictors    : $NPREDICTORS
  runs           : $RUN_START .. $RUN_END  ($total_runs runs)
  runs/task      : $RUNS_PER_TASK  ->  array 0-$array_max ($ntasks tasks)
  skip_setup     : $SKIP_SETUP
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

# ---- 5. Submit setup job (unless skip_setup=1) --------------------------------
if [[ "$SKIP_SETUP" == "1" ]]; then
  echo "skip_setup=1: skipping setup job."
  echo "  NOTE: src/bhCRR.so must already exist and the renv library must be current."
  echo "  If in doubt, omit skip_setup=1 so the setup job runs first."
  array_dep=""
else
  setup_jid=$(submit sbatch --parsable \
    --account="$CHPC_ACCOUNT" --partition="$CHPC_PARTITION" "${CLUSTER_FLAG[@]}" \
    --time="$SB_SETUP_TIME" --mem="$SB_SETUP_MEM" --cpus-per-task="$SB_SETUP_CPUS" \
    --job-name="bhcrr_setup_${SCENARIO}" \
    --output="$REPO_ROOT/sim/chpc/logs/setup_%j.out" \
    --export=ALL \
    "$HERE/setup.slurm")
  setup_jid="${setup_jid%%;*}"
  echo "Submitted setup job  : job $setup_jid"
  array_dep="--dependency=afterok:${setup_jid}"
fi

# ---- 6. Submit the run array -------------------------------------------------
run_jid=$(submit sbatch --parsable \
  --account="$CHPC_ACCOUNT" --partition="$CHPC_PARTITION" "${CLUSTER_FLAG[@]}" \
  --time="$SB_TIME" --mem="$SB_MEM" --cpus-per-task="$SB_CPUS" \
  --array="0-${array_max}" \
  ${array_dep:+"$array_dep"} \
  --job-name="bhcrr_run_${SCENARIO}" \
  --output="$REPO_ROOT/sim/chpc/logs/run_%A_%a.out" \
  --export=ALL \
  "$HERE/run_array.slurm")
# With --clusters, `sbatch --parsable` returns "jobid;cluster"; keep just the id
# so the dependency below is a clean numeric job ID.
run_jid="${run_jid%%;*}"
echo "Submitted run array  : job $run_jid  (array 0-$array_max)"

# ---- 7. Submit the parse job, gated on the whole array succeeding -------------
parse_jid=$(submit sbatch --parsable \
  --account="$CHPC_ACCOUNT" --partition="$CHPC_PARTITION" "${CLUSTER_FLAG[@]}" \
  --time="$SB_PARSE_TIME" --mem="$SB_MEM" --cpus-per-task=1 \
  --dependency="afterok:${run_jid}" \
  --job-name="bhcrr_parse_${SCENARIO}" \
  --output="$REPO_ROOT/sim/chpc/logs/parse_%j.out" \
  --export=ALL \
  "$HERE/parse.slurm")
parse_jid="${parse_jid%%;*}"
echo "Submitted parse job  : job $parse_jid  (runs after array $run_jid succeeds)"

echo
squeue_cmd="squeue -u \$USER"
[[ -n "${CHPC_CLUSTER:-}" ]] && squeue_cmd="squeue -M $CHPC_CLUSTER -u \$USER"
echo "Track with:  $squeue_cmd"
echo "Results will appear in: $REPO_ROOT/sim/results/$SCENARIO/parsed/"
