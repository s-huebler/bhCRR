#!/usr/bin/env bash
# =============================================================================
# sim/timing/submit_timing.sh  --  submit the ClinTALL timing array
#
# Usage (from the repo root on CHPC):
#   ./sim/timing/submit_timing.sh tag=clintall-v1
#   ./sim/timing/submit_timing.sh tag=clintall-v1 reps=3 init=LASSO_bic
#   ./sim/timing/submit_timing.sh tag=clintall-v1 time=02:00:00 mem=16G
#   ./sim/timing/submit_timing.sh tag=clintall-v1 skip_setup=1
#   ./sim/timing/submit_timing.sh tag=clintall-v1 array=50,51,52  # retry specific tasks
#   ./sim/timing/submit_timing.sh dryrun=1 tag=clintall-v1
#
# Required:
#   tag          short label written into every result row (no spaces or slashes)
#
# Optional overrides (key=value, no whitespace around =):
#   reps         replicates per p value          (default 3)
#   init         INIT_METHOD: LASSO_bic or zero  (default LASSO_bic)
#   time         walltime per task               (default 01:00:00)
#   mem          memory per task                 (default 8G)
#   array        explicit 0-based task list/range, e.g. 0-29 or 50,51,52
#                default: 0..(G*REPS-1) for the first submission
#   throttle     max simultaneous tasks          (default 30)
#   constraint   SLURM constraint flag           (default empty)
#   skip_setup   1 = skip setup.slurm            (default 0)
#   dryrun       1 = print commands, no submit   (default 0)
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_SIM_CHPC="$(cd "$HERE/../chpc" && pwd)"

# ---- 1. Parse key=value overrides -------------------------------------------
TAG=""
REPS=""
INIT_OVERRIDE=""
SS0_OVERRIDE=""
SS1_OVERRIDE=""
TIME_OVERRIDE=""
MEM_OVERRIDE=""
ARRAY_OVERRIDE=""
THROTTLE_OVERRIDE=""
CONSTRAINT_OVERRIDE=""
SKIP_SETUP=0
DRYRUN=0

_ALLOWED_KEYS="TAG REPS INIT SS0 SS1 TIME MEM ARRAY THROTTLE CONSTRAINT SKIP_SETUP DRYRUN"

for arg in "$@"; do
  case "$arg" in
    *=*)
      key="${arg%%=*}"; val="${arg#*=}"
      key="$(printf '%s' "$key" | tr '[:lower:]' '[:upper:]' | tr -d '[:space:]')"
      if ! echo " $_ALLOWED_KEYS " | grep -qw "$key"; then
        echo "submit_timing.sh: unknown key '$key'" >&2
        echo "  Allowed keys: $(echo "$_ALLOWED_KEYS" | tr ' ' ', ')" >&2
        exit 1
      fi
      case "$key" in
        TAG)        TAG="$val"            ;;
        REPS)       REPS="$val"           ;;
        INIT)       INIT_OVERRIDE="$val"  ;;
        SS0)        SS0_OVERRIDE="$val"   ;;
        SS1)        SS1_OVERRIDE="$val"   ;;
        TIME)       TIME_OVERRIDE="$val"  ;;
        MEM)        MEM_OVERRIDE="$val"   ;;
        ARRAY)      ARRAY_OVERRIDE="$val" ;;
        THROTTLE)   THROTTLE_OVERRIDE="$val" ;;
        CONSTRAINT) CONSTRAINT_OVERRIDE="$val" ;;
        SKIP_SETUP) SKIP_SETUP="$val"     ;;
        DRYRUN)     DRYRUN="$val"         ;;
      esac
      ;;
    *)
      echo "submit_timing.sh: unrecognized argument '$arg'" >&2
      echo "  Every argument must be key=value with no whitespace." >&2
      exit 1
      ;;
  esac
done
unset _ALLOWED_KEYS

# ---- 2. Source shared config (provides REPO_ROOT, cluster, etc.) -----------
# shellcheck source=/dev/null
source "$REPO_SIM_CHPC/config.sh"

die() { echo "submit_timing.sh: $*" >&2; exit 1; }

# ---- 3. Validate required argument -----------------------------------------
[[ -n "$TAG" ]] || die "tag is required -- pass tag=<label> on the command line."
[[ "$TAG" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || \
  die "TAG '$TAG' is invalid. Use letters, digits, dots, hyphens, underscores."

# ---- 4. Resolve derived settings -------------------------------------------
REPS="${REPS:-3}"
INIT_METHOD="${INIT_OVERRIDE:-LASSO_bic}"
[[ "$INIT_METHOD" =~ ^(LASSO_bic|zero)$ ]] || \
  die "init must be 'LASSO_bic' or 'zero', got '$INIT_METHOD'."
SS0="${SS0_OVERRIDE:-0.03}"
SS1="${SS1_OVERRIDE:-0.3}"
TASK_TIME="${TIME_OVERRIDE:-01:00:00}"
TASK_MEM="${MEM_OVERRIDE:-8G}"
THROTTLE="${THROTTLE_OVERRIDE:-30}"
CONSTRAINT="${CONSTRAINT_OVERRIDE:-}"

check_allocation || exit 1
[[ -d "$REPO_ROOT" ]] || die "REPO_ROOT does not exist: $REPO_ROOT"

# Grid and input paths
GRID_FILE="$REPO_ROOT/sim/timing/grid.tsv"
[[ -f "$GRID_FILE" ]] || \
  die "grid.tsv not found: $GRID_FILE\nRun prep_clintall_timing.r first."

G=$(tail -n +2 "$GRID_FILE" | wc -l | tr -d ' ')
[[ "$G" -gt 0 ]] || die "grid.tsv is empty."

TIMING_INPUT="$REPO_ROOT/sim/timing/input/clintall_timing_input.rds"
[[ -f "$TIMING_INPUT" ]] || \
  die "Input file not found: $TIMING_INPUT\nRun prep_clintall_timing.r first."

TIMING_OUT_DIR="${SCRATCH_BASE}/timing/${TAG}"
COLLECT_OUT="${TIMING_OUT_DIR}/clintall_timing_${TAG}.csv"

# Array range
if [[ -n "$ARRAY_OVERRIDE" ]]; then
  ARRAY_SPEC="$ARRAY_OVERRIDE"
else
  NTASKS=$(( G * REPS ))
  ARRAY_SPEC="0-$(( NTASKS - 1 ))%${THROTTLE}"
fi

# Optional --constraint flag
CONSTRAINT_FLAG=()
[[ -n "$CONSTRAINT" ]] && CONSTRAINT_FLAG=(--constraint="$CONSTRAINT")

# Cluster routing
CLUSTER_FLAG=()
[[ -n "${CHPC_CLUSTER:-}" ]] && CLUSTER_FLAG=(--clusters="$CHPC_CLUSTER")

mkdir -p "$REPO_ROOT/sim/chpc/logs" || true

# ---- 5. Banner --------------------------------------------------------------
cat <<EOF
------------------------------------------------------------
bhCRR ClinTALL timing submission
  tag            : $TAG
  repo_root      : $REPO_ROOT
  timing input   : $TIMING_INPUT
  output dir     : $TIMING_OUT_DIR
  grid           : $GRID_FILE  (G=$G)
  R module       : $R_MODULE
  BHCRR_LOAD     : ${BHCRR_LOAD:-cluster}
  account/part.  : $CHPC_ACCOUNT / $CHPC_PARTITION  (cluster: ${CHPC_CLUSTER:-<login>})
  init_method    : $INIT_METHOD
  ss             : ($SS0, $SS1)
  reps           : $REPS  ->  G*REPS=$(( G * REPS )) tasks
  array_spec     : $ARRAY_SPEC
  time / mem     : $TASK_TIME / $TASK_MEM
  skip_setup     : $SKIP_SETUP
------------------------------------------------------------
EOF

submit() {
  if [[ "$DRYRUN" != "0" ]]; then
    echo "[dryrun] $*" >&2
    echo "DRYRUN_JOBID"
  else
    "$@"
  fi
}

if [[ "$DRYRUN" == "0" ]]; then
  mkdir -p "$TIMING_OUT_DIR"
fi

# ---- 6. Setup job (unless skip_setup=1) ------------------------------------
if [[ "$SKIP_SETUP" == "1" ]]; then
  echo "skip_setup=1: skipping setup job."
  array_dep=""
else
  setup_jid=$(submit sbatch --parsable \
    --account="$CHPC_ACCOUNT" --partition="$CHPC_PARTITION" "${CLUSTER_FLAG[@]}" \
    --time="$SB_SETUP_TIME" --mem="$SB_SETUP_MEM" --cpus-per-task="$SB_SETUP_CPUS" \
    --job-name="bhcrr_setup_timing" \
    --output="$REPO_ROOT/sim/chpc/logs/setup_%j.out" \
    --export=ALL \
    "$REPO_SIM_CHPC/setup.slurm")
  setup_jid="${setup_jid%%;*}"
  echo "Submitted setup job  : job $setup_jid"
  array_dep="--dependency=afterok:${setup_jid}"
fi

# ---- 7. Timing array --------------------------------------------------------
run_jid=$(submit sbatch --parsable \
  --account="$CHPC_ACCOUNT" --partition="$CHPC_PARTITION" "${CLUSTER_FLAG[@]}" \
  --time="$TASK_TIME" --mem="$TASK_MEM" --cpus-per-task=1 \
  --array="$ARRAY_SPEC" \
  ${array_dep:+"$array_dep"} \
  "${CONSTRAINT_FLAG[@]}" \
  --job-name="bhcrr_timing_${TAG}" \
  --output="$REPO_ROOT/sim/chpc/logs/timing_%A_%a.out" \
  --export=ALL,TIMING_TAG="$TAG",TIMING_INPUT="$TIMING_INPUT",TIMING_OUT_DIR="$TIMING_OUT_DIR",REPS="$REPS",INIT_METHOD="$INIT_METHOD",SS0="$SS0",SS1="$SS1" \
  "$HERE/timing.slurm")
run_jid="${run_jid%%;*}"
echo "Submitted timing array: job $run_jid  (array $ARRAY_SPEC)"

# ---- 8. Collect job (afterany so it runs even if some tasks fail) -----------
collect_jid=$(submit sbatch --parsable \
  --account="$CHPC_ACCOUNT" --partition="$CHPC_PARTITION" "${CLUSTER_FLAG[@]}" \
  --time="00:10:00" --mem=4G --cpus-per-task=1 \
  --dependency="afterany:${run_jid}" \
  --job-name="bhcrr_collect_${TAG}" \
  --output="$REPO_ROOT/sim/chpc/logs/collect_%j.out" \
  --export=ALL,TIMING_TAG="$TAG",TIMING_OUT_DIR="$TIMING_OUT_DIR",GRID_FILE="$GRID_FILE",REPS="$REPS",ARRAY_JOB_ID="$run_jid" \
  "$HERE/collect.slurm")
collect_jid="${collect_jid%%;*}"
echo "Submitted collect job : job $collect_jid  (runs after array $run_jid)"

echo
squeue_cmd="squeue -u \$USER"
[[ -n "${CHPC_CLUSTER:-}" ]] && squeue_cmd="squeue -M $CHPC_CLUSTER -u \$USER"
echo "Track with:  $squeue_cmd"
echo
echo "Check timing-array jobs with sacct:"
echo "  sacct -M ${CHPC_CLUSTER:-<cluster>} -j ${run_jid} --format=JobID,JobName,State,Elapsed,Timelimit,MaxRSS,ExitCode -X"
echo
echo "Results will appear in: $TIMING_OUT_DIR/"
echo "Collected CSV at:        ${COLLECT_OUT}"
echo
echo "Pull results to local with:"
echo "  scp -r '<user>@notchpeak.chpc.utah.edu:${TIMING_OUT_DIR}/clintall_timing_${TAG}.csv' \\"
echo "      ~/Documents/bhCRR/sim/timing/results/"
