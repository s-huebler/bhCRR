# bhCRR on CHPC (SLURM)

Run the `batch_run_complex` simulation as a parallel SLURM job array on the
University of Utah CHPC, then automatically collect the per-run results into
`Sim/chpc_results/<scenario>/` to pull back to local for figures.

## Files

| File | Role |
|------|------|
| `config.sh` | The **one** file you edit/source. Paths, R module, SLURM account, and default sim params. |
| `run_sim.sh` | Entry point. Parses `key=value` overrides, submits the array + a dependent parse job. |
| `run_array.slurm` | Array task: runs a chunk of run indices, writes one `.Rdata` per run to scratch. |
| `parse.slurm` | Runs once after the array succeeds; builds `betas_p<P>_<scenario>`. |
| `R/run_batch.R` | CHPC-portable version of `Sim/code/batch_run_complex.R` (env-driven, renv). |
| `R/parse_batch.R` | CHPC-portable version of `Sim/code/parse_batch_run.R`; writes CSV + `.Rdata`. |

## One-time setup

1. Clone/pull the repo on CHPC and note its path (e.g. `~/bhCRR`).
2. Allocation is prefilled in `config.sh` for the `qiaox` / lonepeak setup:
   - `CHPC_ACCOUNT` — defaults to `qiaox`
   - `CHPC_PARTITION` — defaults to `lonepeak` (shared, queues faster than a
     whole node). Override examples are in the config comments.
   - `CHPC_CLUSTER` — defaults to `lonepeak`; each CHPC cluster runs its own
     Slurm controller. Set to `""` for your login cluster's scheduler.
   - `R_MODULE` — confirm with `module spider R` (renv.lock pins **R 4.5.1**)
   - `SCRATCH_BASE` — defaults to `/scratch/general/vast/$USER/bhCRR`
3. The first array task runs `renv::restore()`; subsequent tasks are fast.

## Everyday use

```bash
ssh <you>@notchpeak.chpc.utah.edu
cd ~/bhCRR
source chpc/config.sh
./chpc/run_sim.sh npredictors=15 nobs=200 run_start=1 run_end=100 \
  scenario=scenarioA sb_time=1:00:00
```

`scenario=<name>` is **required** on every submission. All other sim parameters
have defaults in `config.sh` and can be overridden the same way (`key=value`,
**no spaces** around `=`, keys are case-insensitive).

```bash
# Extend an existing scenario with more runs (collision guard will catch overlaps)
./chpc/run_sim.sh scenario=scenarioA npredictors=15 nobs=200 run_start=101 run_end=200

# Add a second npredictors value to the same scenario
./chpc/run_sim.sh scenario=scenarioA npredictors=50 nobs=200 run_start=1 run_end=100

# Preview without submitting
./chpc/run_sim.sh scenario=scenarioA nobs=200 dryrun=1
```

`dryrun=1` prints the banner and resolved sbatch commands without submitting
anything. It also checks the scenario manifest if one already exists, so it is
the way to confirm knobs before a real submission. Off-cluster (Mac), point
`SCRATCH_BASE` at a local writable path:

```bash
SCRATCH_BASE=/tmp/bhcrr_dryrun ./chpc/run_sim.sh scenario=scenarioA dryrun=1
```

Named-vector keys (`block_props`, `block_rho`, `latent_q`) use `name=val,name=val`
format. The parser splits only on the **first** `=`, so extra `=` inside the
value are safe:

```bash
./chpc/run_sim.sh scenario=scenarioA block_rho=indep=0.0,weak=0.5,strong=0.8
```

## Scenarios

A **scenario** is a named, immutable set of data-generating parameters. Pick a
short descriptive slug for each distinct DGP you want to study; the slug becomes
the directory name on scratch and in results.

### Slug rules

The scenario name must match `^[A-Za-z0-9][A-Za-z0-9.-]*$`:
letters, digits, dots, and hyphens only, starting with a letter or digit.
**No underscores** — underscores delimit fields in the
`p<P>_<scenario>_run<N>.Rdata` filename format, so allowing them in the name
would make that pattern ambiguous when parsing results.

Good: `n200-baseline`, `n200.strongBlock`, `v2`, `n500highCens`  
Bad: `n200_baseline` (underscore), `-start` (leading hyphen)

### Scenario signature

On first submission, `run_sim.sh` writes a manifest (`scenario.env`) to the
scenario's scratch directory recording the **12 DGP keys** that define the
scenario:

| Key | Default in `config.sh` |
|-----|------------------------|
| `NOBS` | `200` |
| `ZERO_GAP_TARGET` | `0.1` |
| `BETA1_ACTIVE` | `0.40,-0.50,0.60,0.75,-0.80` |
| `BETA2_ACTIVE` | `neg_beta1` |
| `ACTIVE_BLOCK` | `strong,strong,weak,weak,indep` |
| `BLOCK_PROPS` | `indep=0.50,weak=0.25,strong=0.25` |
| `BLOCK_RHO` | `indep=0.00,weak=0.30,strong=0.60` |
| `LATENT_TYPE` | `gaussian` |
| `LATENT_Q` | `indep=0.5,weak=0.4,strong=0.5` |
| `CENS_RATE` | `0.05` |
| `U_MIN` | `100` |
| `U_MAX` | `100` |

Note that **`n` (NOBS) is part of the scenario** — it describes the
data-generating setup, not the computational scope. It no longer appears in
filenames; the scenario slug implicitly encodes it.

The following are **free to vary within a scenario** and are not recorded in the
manifest:

- `NPREDICTORS` — which p values to run (one subdir per p under the scenario dir)
- `RUN_START`, `RUN_END`, `RUNS_PER_TASK` — how many runs and how to chunk them
- All SLURM/account knobs (`sb_time`, `sb_mem`, `sb_cpus`, `CHPC_ACCOUNT`, etc.)

### Resubmitting to the same scenario

Every resubmission to an existing scenario re-reads `scenario.env` and compares
each of the 12 signature keys against the values you just passed. **Any mismatch
is a hard stop** before `sbatch` is called:

```
run_sim.sh: scenario mismatch -- NOBS:
  recorded : 200
  current  : 500
run_sim.sh: Scenario 'scenarioA' has recorded settings that differ from what
you passed (see above). Fix the mismatch or choose a new scenario name with
scenario=<new>.
```

To resolve: either correct the override so it matches the recorded values, or
start a new scenario (`scenario=scenarioB nobs=500 ...`). There is no override
flag — the mismatch check is a hard stop by design.

### Collision guard

Before submitting, `run_sim.sh` checks whether any of the requested run indices
already have `.Rdata` output under `$SCRATCH_BASE/<scenario>/p<P>/`. If any do,
the submission is refused:

```
run_sim.sh: collision for p=15 -- existing runs: 1 2 3
  Pass overwrite=1 to overwrite, or run_start=4 to append.
run_sim.sh: Refusing to submit due to run collisions (see above).
```

To append new runs without touching existing ones, use `run_start=<next free>`.
To recompute existing indices (e.g., after a code fix), pass `overwrite=1`:

```bash
./chpc/run_sim.sh scenario=scenarioA npredictors=15 nobs=200 \
  run_start=1 run_end=10 overwrite=1
```

## Directory layout

### Scratch (per scenario)

```
$SCRATCH_BASE/
  scenarioA/
    scenario.env           ← manifest: 12 DGP keys + creation timestamp + git SHA
    p15/
      p15_scenarioA_run1.Rdata
      p15_scenarioA_run2.Rdata
      ...
    p50/
      p50_scenarioA_run1.Rdata
      ...
```

The `scenario.env` manifest is also copied into `Sim/chpc_results/<scenario>/`
by the parse job, so the results folder is self-describing without the scratch
tree.

### Results (in the repo)

```
Sim/chpc_results/
  scenarioA/
    scenario.env
    betas_p15_scenarioA.csv
    betas_p15_scenarioA.Rdata        ← loads object `betas_p15`
    model_performance_p15_scenarioA.csv
    false_positives_p15_scenarioA.csv
    betas_p50_scenarioA.csv
    ...
    model_performance_scenarioA.csv  ← stacked across all p, only when >1 p group
```

All three per-group tables carry a `scenario` column as their first column.
`nobs` is sourced from `result$meta$nobs` (not the filename) and is verified to
be identical across every run before the parse completes.

### A note on old scratch files

`parse_batch.R` only matches `p<P>_<scenario>_run<N>.Rdata`. Files in the old
`n<N>_p<P>_run<N>.Rdata` format are silently skipped with a log message. If you
have old scratch directories under `$SCRATCH_BASE/runs_n<N>/`, they will not be
parsed and can be removed once no longer needed.

## What happens end-to-end

1. `run_sim.sh` validates the scenario name and manifest, checks for collisions,
   then splits `run_start..run_end` into `runs_per_task`-sized chunks and submits
   a **job array** (`--array=0-N`). Job name: `bhcrr_run_<scenario>`.
2. Each array task writes `p<P>_<scenario>_run<N>.Rdata` into
   `$SCRATCH_BASE/<scenario>/p<P>/`.
3. A **parse job** (`bhcrr_parse_<scenario>`) is submitted with
   `--dependency=afterok` on the array. When every task succeeds it reads all
   matching `.Rdata` files, runs the consistency checks, and writes results to
   `Sim/chpc_results/<scenario>/`.

Runs at different `npredictors` go into separate `p<P>` subdirs because the
column-to-block map is derived from `npredictors`: `X_17` at p=25 is not `X_17`
at p=221, so runs at different p cannot be pooled.

## Pulling results to local

From your Mac:

```bash
scp -r '<you>@notchpeak.chpc.utah.edu:~/bhCRR/Sim/chpc_results/scenarioA' \
    ~/Documents/bhCRR/Sim/chpc_results/
```

Or individual files:

```bash
scp '<you>@notchpeak.chpc.utah.edu:~/bhCRR/Sim/chpc_results/scenarioA/betas_p*.csv' \
    ~/Documents/bhCRR/Sim/chpc_results/scenarioA/
```

## Monitoring & logs

```bash
squeue -u $USER
# or, for lonepeak specifically:
squeue -M lonepeak -u $USER
```

Per-task stdout/stderr land in `chpc/logs/run_<arrayjob>_<task>.out` and
`chpc/logs/parse_<job>.out`.

## Notes

- `set.seed()` is intentionally never called — each run draws fresh data.
- Scratch is not backed up; `Sim/chpc_results/` is the durable output.
- Each array task recompiles the Rcpp sources via `sourceCpp`. If this becomes a
  bottleneck, precompiling `src/*.so` once is a future optimization.
