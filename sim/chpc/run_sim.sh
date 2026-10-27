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

# Allowlist of keys that run_sim.sh accepts. Derived from config.sh's sweep
# and resource knobs; does not include DGP keys (those belong in the scenario
# file) or internal-only shell vars.
_ALLOWED_KEYS="SCENARIO NPREDICTORS RUN_START RUN_END RUNS_PER_TASK \
  SB_TIME SB_MEM SB_CPUS SB_PARSE_TIME \
  SB_SETUP_TIME SB_SETUP_MEM SB_SETUP_CPUS \
  SCRATCH_BASE R_MODULE USE_RENV BHCRR_LOAD \
  CHPC_ACCOUNT CHPC_PARTITION CHPC_CLUSTER"

_DGP_KEYS="NOBS ZERO_GAP_TARGET BETA1_ACTIVE BETA2_ACTIVE ACTIVE_BLOCK \
  BLOCK_PROPS BLOCK_RHO LATENT_TYPE LATENT_Q CENS_RATE U_MIN U_MAX"

for arg in "$@"; do
  case "$arg" in
    *=*)
      key="${arg%%=*}"; val="${arg#*=}"
      key="$(printf '%s' "$key" | tr '[:lower:]' '[:upper:]' | tr -d '[:space:]')"
      if [[ "$key" == "DRYRUN" ]];      then DRYRUN="$val";      continue; fi
      if [[ "$key" == "OVERWRITE" ]];   then OVERWRITE="$val";   continue; fi
      if [[ "$key" == "SKIP_SETUP" ]];  then SKIP_SETUP="$val";  continue; fi
      # Check DGP keys first -- specific, actionable error message
      if echo " $_DGP_KEYS " | grep -qw "$key"; then
        echo "run_sim.sh: '$key' is a DGP parameter -- it belongs in sim/scenarios/\${SCENARIO}.R" >&2
        echo "  DGP parameters are no longer accepted on the command line." >&2
        echo "  Create or edit the scenario file: sim/scenarios/\${SCENARIO:-<name>}.R" >&2
        echo "  (new_scenario() copies the template; see sim/README.md)" >&2
        exit 1
      fi
      # Validate against allowlist
      if ! echo " $_ALLOWED_KEYS " | grep -qw "$key"; then
        echo "run_sim.sh: unknown key '$key' -- not in the allowed list." >&2
        echo "  Allowed keys: $(echo $_ALLOWED_KEYS | tr ' ' ', ')" >&2
        exit 1
      fi
      export "$key=$val"
      # Reject malformed comma lists: a leading, trailing, or doubled comma
      # signals a typo, most commonly from writing spaces after commas so the
      # shell split the argument. "npredictors=50, 200, 1000" becomes three
      # words; the first word lands here as "npredictors=50," with a trailing
      # comma.
      if [[ "$key" == "NPREDICTORS" ]]; then
        if [[ "$val" =~ ^,|,,|,$ ]]; then
          echo "run_sim.sh: NPREDICTORS value '$val' has an empty element (leading, trailing, or doubled comma)." >&2
          echo "  This usually means spaces were written after commas, causing the shell to" >&2
          echo "  split the argument into multiple words. Only the first word reached here." >&2
          echo "  Wrong:  npredictors=50, 200, 1000   <- shell sees 3 words; first is 'npredictors=50,'" >&2
          echo "  Right:  npredictors=50,200,1000" >&2
          exit 1
        fi
      fi
      ;;
    *)
      echo "run_sim.sh: unrecognized argument '$arg'" >&2
      echo "  Every argument must be a key=value token with no whitespace anywhere in the token." >&2
      echo "  Wrong:  npredictors=50, 200, 1000   <- shell splits this into 3 words; '200,' and '1000' land here" >&2
      echo "  Right:  npredictors=50,200,1000" >&2
      exit 1
      ;;
  esac
done
unset _ALLOWED_KEYS _DGP_KEYS

# ---- 2. Fill in anything not overridden from the single config file ----------
# shellcheck source=/dev/null
source "$HERE/config.sh"

# ---- 2b. Guard: DGP env vars must NOT be set --------------------------------
# If a shell still has the pre-12b config.sh sourced, all 12 DGP vars are
# exported. run_array.r would catch them, but only after the queue wait.
# Detect and fail here before any sbatch.
_dgp_vars=(NOBS ZERO_GAP_TARGET BETA1_ACTIVE BETA2_ACTIVE ACTIVE_BLOCK
           BLOCK_PROPS BLOCK_RHO LATENT_TYPE LATENT_Q CENS_RATE U_MIN U_MAX)
_set_dgp=()
for _v in "${_dgp_vars[@]}"; do
  [[ -n "${!_v:-}" ]] && _set_dgp+=("$_v")
done
if [[ ${#_set_dgp[@]} -gt 0 ]]; then
  echo "run_sim.sh: the following DGP environment variable(s) are set:" >&2
  printf '  %s\n' "${_set_dgp[@]}" >&2
  echo "  DGP parameters are defined in sim/scenarios/${SCENARIO:-<name>}.R" >&2
  echo "  Unset them before submitting:" >&2
  printf '  unset %s\n' "${_set_dgp[@]}" >&2
  exit 1
fi
unset _dgp_vars _set_dgp _v

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

# ---- 3b. Manifest write/check via R -----------------------------------------
# Write or check scenario.env before any sbatch so a DGP mismatch fails before
# anything queues. Uses scenario.r (no bhCRR package needed -- scenario.r has
# no intra-sim or bhCRR dependencies). Permitted on login nodes (lightweight
# data analysis, no compilation, ~2 s). Hard-fails if R is unavailable.
if [[ "$DRYRUN" == "0" ]]; then
  command -v Rscript >/dev/null 2>&1 || \
    { module load "$R_MODULE" 2>/dev/null || true; }
  if ! command -v Rscript >/dev/null 2>&1; then
    die "Rscript not found after 'module load $R_MODULE'. Ensure the R module is available before submitting."
  fi
  echo "[run_sim.sh] Writing/checking scenario manifest..."
  _r_tmp=$(mktemp /tmp/bhcrr_manifest_XXXXX.R)
  cat > "$_r_tmp" <<'REOF'
args <- commandArgs(trailingOnly = TRUE)
manifest_dir  <- args[1]
scenario_path <- args[2]
source(file.path(Sys.getenv("REPO_ROOT"), "sim", "R", "scenario.r"))
env_path <- file.path(manifest_dir, "scenario.env")
if (file.exists(env_path)) {
  check_manifest(scenario_path, manifest_dir)
  message("  manifest OK: ", env_path)
} else {
  if (!dir.exists(manifest_dir)) dir.create(manifest_dir, recursive = TRUE)
  write_manifest(scenario_path, manifest_dir)
  message("  manifest written: ", env_path)
}
REOF
  if ! Rscript --vanilla "$_r_tmp" "$SCRATCH_RUN_DIR" "$scenario_file"; then
    rm -f "$_r_tmp"
    die "Manifest write/check failed. Fix the scenario mismatch before submitting."
  fi
  rm -f "$_r_tmp"
  echo "[run_sim.sh] Manifest OK."
else
  echo "  [dryrun] manifest write/check skipped"
fi

# ---- 3c. Collision guard -----------------------------------------------------
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
