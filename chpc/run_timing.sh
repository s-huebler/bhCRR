#!/usr/bin/env bash
# =============================================================================
# chpc/run_timing.sh  --  submit bhCRR timing jobs, one sbatch per config.
#
# Usage (from the repo root on CHPC):
#   ./chpc/run_timing.sh                                 # all seven configs
#   ./chpc/run_timing.sh tags=n200-p300,n200-p5000       # subset
#   ./chpc/run_timing.sh init_mode=zero                  # override any config var
#   ./chpc/run_timing.sh time_n1309_pall=48:00:00 mem_n1309_pall=192G
#   ./chpc/run_timing.sh dryrun=1                        # print commands, no submit
#
# Resource overrides use underscores for the tag part (hyphens are not valid
# in shell variable names):  n200-p5000 -> time_n200_p5000 / mem_n200_p5000
#
# Notes:
#   * Each config is a separate sbatch call with its own --time and --mem.
#     This lets small jobs (p=300) and large jobs (p=all) sit in the queue
#     without either wasting resources or hitting the walltime limit.
#   * Add dryrun=1 to preview the sbatch commands without submitting.
#   * Overrides are key=value with NO spaces around '='.
#   * Keys are case-insensitive.
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---- 1. Parse key=value overrides into the environment ----------------------
DRYRUN=0
TAGS_FILTER=""   # empty = all

for arg in "$@"; do
  case "$arg" in
    *=*)
      key="${arg%%=*}"; val="${arg#*=}"
      key="$(printf '%s' "$key" | tr '[:lower:]' '[:upper:]' | tr -d '[:space:]')"
      case "$key" in
        DRYRUN)  DRYRUN="$val";       continue ;;
        TAGS)    TAGS_FILTER="$val";  continue ;;
      esac
      export "$key=$val"
      ;;
    *)
      echo "run_timing.sh: ignoring unrecognized argument '$arg' (expected key=value)" >&2
      ;;
  esac
done

# ---- 2. Fill in defaults from timing_config.sh ------------------------------
# shellcheck source=/dev/null
source "$HERE/timing_config.sh"

# ---- 3. Per-tag resource defaults (override with time_<tag>=... mem_<tag>=...)
# Tags use underscores in variable names because hyphens are not valid.
: "${TIME_N200_P300:=00:30:00}"    ; : "${MEM_N200_P300:=8G}"
: "${TIME_N1309_P300:=01:00:00}"   ; : "${MEM_N1309_P300:=8G}"
: "${TIME_N200_P1000:=01:00:00}"   ; : "${MEM_N200_P1000:=16G}"
: "${TIME_N200_P5000:=04:00:00}"   ; : "${MEM_N200_P5000:=32G}"
: "${TIME_N200_P10000:=08:00:00}"  ; : "${MEM_N200_P10000:=64G}"
: "${TIME_N200_PALL:=12:00:00}"    ; : "${MEM_N200_PALL:=96G}"
: "${TIME_N1309_PALL:=24:00:00}"   ; : "${MEM_N1309_PALL:=128G}"

# ---- 4. Sanity checks -------------------------------------------------------
die() { echo "run_timing.sh: $*" >&2; exit 1; }
check_allocation || exit 1
[[ -d "$REPO_ROOT" ]] || die "REPO_ROOT does not exist: $REPO_ROOT"

# Warn if the timing data hasn't been materialized yet.
if [[ "$DRYRUN" == "0" && ! -d "$TIMING_DATA_DIR" ]]; then
  echo "WARNING: TIMING_DATA_DIR does not exist: $TIMING_DATA_DIR" >&2
  echo "         Run prep_timing_data.R (with TIMING_DATA_DIR=$TIMING_DATA_DIR)" >&2
  echo "         before submitting, or the jobs will fail at data load." >&2
fi

CLUSTER_FLAG=()
[[ -n "${CHPC_CLUSTER:-}" ]] && CLUSTER_FLAG=(--clusters="$CHPC_CLUSTER")

mkdir -p "$REPO_ROOT/chpc/logs" "$TIMING_OUT_DIR"

# ---- 5. Build the tag list --------------------------------------------------
# Canonical ordering -- indices match timing.slurm's _TIMING_TAGS array.
declare -a ALL_TAGS=(
  "n200-p300"     # 0
  "n1309-p300"    # 1
  "n200-p1000"    # 2
  "n200-p5000"    # 3
  "n200-p10000"   # 4
  "n200-pall"     # 5
  "n1309-pall"    # 6
)

# Filter to the requested subset (tags= arg), or use all.
declare -a SUBMIT_TAGS=()
if [[ -n "$TAGS_FILTER" ]]; then
  IFS=',' read -ra _req <<< "$TAGS_FILTER"
  for t in "${_req[@]}"; do
    found=0
    for canon in "${ALL_TAGS[@]}"; do
      [[ "$t" == "$canon" ]] && { SUBMIT_TAGS+=("$t"); found=1; break; }
    done
    [[ "$found" -eq 0 ]] && die "Unknown tag '$t'. Valid tags: ${ALL_TAGS[*]}"
  done
else
  SUBMIT_TAGS=("${ALL_TAGS[@]}")
fi

# ---- 6. Banner --------------------------------------------------------------
# Helper: look up time/mem for a given tag.
tag_var_suffix() { printf '%s' "${1//-/_}" | tr '[:lower:]' '[:upper:]'; }
tag_time() { local v="TIME_$(tag_var_suffix "$1")"; echo "${!v}"; }
tag_mem()  { local v="MEM_$(tag_var_suffix "$1")";  echo "${!v}"; }

cat <<EOF
------------------------------------------------------------
bhCRR timing submission
  repo_root       : $REPO_ROOT
  timing data     : $TIMING_DATA_DIR
  results         : $TIMING_OUT_DIR
  R module        : $R_MODULE   (renv restore: $USE_RENV)
  account/part.   : $CHPC_ACCOUNT / $CHPC_PARTITION  (cluster: ${CHPC_CLUSTER:-<login>})
  init_mode       : $INIT_MODE
  ss              : ($SS0, $SS1)
  theta_a/b       : $THETA_A / $THETA_B
  maxit           : $MAXIT   epsilon: $EPSILON
  -- per-tag resources --
EOF

for tag in "${SUBMIT_TAGS[@]}"; do
  printf "  %-16s  time=%-10s  mem=%s\n" "$tag" "$(tag_time "$tag")" "$(tag_mem "$tag")"
done
echo "------------------------------------------------------------"

# ---- 7. Submit helper -------------------------------------------------------
submit() {
  if [[ "$DRYRUN" != "0" ]]; then
    echo "[dryrun] $*" >&2
    echo "DRYRUN_JOBID"
  else
    "$@"
  fi
}

# ---- 8. Submit one sbatch per tag -------------------------------------------
for tag in "${SUBMIT_TAGS[@]}"; do
  t_wall="$(tag_time "$tag")"
  t_mem="$(tag_mem "$tag")"

  jid=$(submit sbatch --parsable \
    --account="$CHPC_ACCOUNT" --partition="$CHPC_PARTITION" "${CLUSTER_FLAG[@]}" \
    --time="$t_wall" --mem="$t_mem" --cpus-per-task=1 \
    --job-name="bhcrr_timing_${tag}" \
    --output="$REPO_ROOT/chpc/logs/timing_${tag}_%j.out" \
    --export=ALL,TIMING_TAG="$tag" \
    "$HERE/timing.slurm")
  jid="${jid%%;*}"
  echo "Submitted $tag : job $jid  (time=$t_wall  mem=$t_mem)"
done

echo
squeue_cmd="squeue -u \$USER"
[[ -n "${CHPC_CLUSTER:-}" ]] && squeue_cmd="squeue -M $CHPC_CLUSTER -u \$USER"
echo "Track with:  $squeue_cmd"
echo "Results:     $TIMING_OUT_DIR/timing_<tag>.rds"
echo "Logs:        $REPO_ROOT/chpc/logs/timing_<tag>_<jobid>.out"
echo "Grep summary lines with:  grep '^TIMING |' chpc/logs/timing_*.out"
