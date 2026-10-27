# bhCRR on CHPC (SLURM)

Run the bhCRR simulation as a parallel SLURM job array on the University of
Utah CHPC, then automatically collect the per-run results into
`sim/results/<scenario>/parsed/` to pull back to local for figures.

## Files

| File | Role |
|------|------|
| `config.sh` | The **one** file you source before submitting. Cluster identity, paths, sweep scope, and renv policy. DGP parameters live in the scenario file, not here. |
| `run_sim.sh` | Entry point. Parses `key=value` overrides; submits setup → array → parse. |
| `setup.slurm` | One-shot prepare job: restores the renv library and compiles `src/bhCRR.so`. Must run before the array (array depends on it via `afterok`). |
| `run_array.slurm` | Array task: runs a chunk of run indices; writes one `.Rdata` per run to scratch. |
| `parse.slurm` | Runs once after the array succeeds; stacks `.Rdata` files into CSV summaries. |
| `R/run_array.r` | Thin wrapper: reads env vars, calls `batch_sweep()`. All logic in `sim/local/batch.R`. |
| `R/parse.r` | Thin wrapper: reads env vars, calls `parse_results()`. All logic in `sim/R/parse.r`. |

## One-time setup

1. Clone/pull the repo on CHPC and note its path (e.g. `~/Documents/bhCRR`).
2. Allocation is prefilled in `config.sh` for the `qiaox` / lonepeak setup:
   - `CHPC_ACCOUNT` — defaults to `biostat-division-np`
   - `CHPC_PARTITION` — defaults to `biostat-division-np` (shared, queues faster
     than a whole node). Override examples are in the config comments.
   - `CHPC_CLUSTER` — defaults to `notchpeak`; each CHPC cluster runs its own
     Slurm controller. Set to `""` for your login cluster's scheduler.
   - `R_MODULE` — confirm with `module spider R` (renv.lock pins **R 4.5.1**)
   - `SCRATCH_BASE` — defaults to `/scratch/general/vast/$USER/bhCRR`
3. The **setup job** handles renv restore and compilation (see [renv & compilation](#renv--compilation) below). It runs once before the array.

## Everyday use

```bash
ssh <you>@notchpeak.chpc.utah.edu
cd ~/Documents/bhCRR
source sim/chpc/config.sh          # optional; run_sim.sh sources it too
module load "$R_MODULE"            # needed if you run R interactively on login

./sim/chpc/run_sim.sh scenario=n200.baseline npredictors=25,50 run_start=1 run_end=100
```

`scenario=<name>` is **required** on every submission. All other parameters
have defaults in `config.sh` and can be overridden with `key=value` (**no
whitespace anywhere in the token**, keys are case-insensitive).

> **Warning:** comma-separated values must contain no spaces.
> `npredictors=50,200,1000` is correct; `npredictors=50, 200, 1000` is wrong
> — the shell splits the latter into three words and `run_sim.sh` exits with
> an error on the unexpected bare words.

```bash
# Extend an existing scenario with more runs (collision guard catches overlaps)
./sim/chpc/run_sim.sh scenario=n200.baseline npredictors=25 run_start=101 run_end=200

# Preview without submitting
./sim/chpc/run_sim.sh scenario=n200.baseline dryrun=1

# Skip the setup job when bhCRR.so is already compiled and library is current
./sim/chpc/run_sim.sh scenario=n200.baseline skip_setup=1 run_start=1 run_end=10
```

`dryrun=1` prints the banner and resolved sbatch commands without submitting
anything. Off-cluster (Mac), point `SCRATCH_BASE` at a local writable path:

```bash
SCRATCH_BASE=/tmp/bhcrr_dryrun ./sim/chpc/run_sim.sh scenario=n200.baseline dryrun=1
```

## Scenarios

A **scenario** is a named, immutable set of data-generating parameters defined
in `sim/scenarios/<name>.R`. Pick a short descriptive slug for each distinct
DGP you want to study; the slug becomes the directory name on scratch and in
results.

### Creating a scenario

```r
source("sim/R/load.r")
load_sim()
new_scenario("n200.baseline", dir = "sim/scenarios")
# Opens sim/scenarios/n200.baseline.R — fill in every [PLACEHOLDER] before submitting.
```

The scenario file assigns exactly one object named `scenario` (a plain R list)
and must not call `library()`, `require()`, `source()`, or `set.seed()`.
`batch_sweep()` reads it via `scenario_from_file()` which enforces these rules.

### Slug rules

The scenario name must match `^[A-Za-z0-9][A-Za-z0-9.-]*$`:
letters, digits, dots, and hyphens only, starting with a letter or digit.
**No underscores** — underscores delimit fields in the
`p<P>_<scenario>_run<N>.Rdata` filename format, so allowing them in the name
would make that pattern ambiguous when parsing results.

Good: `n200-baseline`, `n200.strongBlock`, `v2`  
Bad: `n200_baseline` (underscore), `-start` (leading hyphen)

### Scenario identity and manifest

The scenario file is the **single source of truth** for all DGP parameters
(sample size, coefficients, block structure, censoring). `config.sh` carries
none of them.

On first run `batch_sweep()` writes a manifest (`scenario.env`) to the
scenario's scratch directory. Every subsequent submission re-reads the manifest
and hard-stops on any DGP mismatch before `sbatch` is called. To resolve: fix
the scenario file or choose a new scenario name.

### Collision guard

Before submitting, `run_sim.sh` checks whether any requested run indices already
have `.Rdata` output under `$SCRATCH_BASE/<scenario>/p<P>/`. If any do, the
submission is refused. To append, use `run_start=<next free>`; to recompute
existing indices, pass `overwrite=1`.

## renv & compilation

CHPC Arbiter caps login-node jobs at 4 cores / 8 GB. A full `renv::restore()`
of 166 source packages and a `compile_dll()` exceed that budget. **Both run in
the setup job** (`setup.slurm`), submitted before the array and gated via
`--dependency=afterok`.

### USE_RENV — three states

| Value | Behavior |
|-------|----------|
| `auto` (default) | Restore only if `renv::status()` reports the project is not synchronized with the lockfile. Right for most resubmissions. |
| `1` | Always run `renv::restore(prompt=FALSE)`. Use after a `renv.lock` change. |
| `0` | Skip restore entirely. Only safe when the library is known to be current. |

Set on the `run_sim.sh` command line: `./sim/chpc/run_sim.sh ... use_renv=1`

### renv library location

The renv library lives under `~/.cache/R/renv/library/` on the compute nodes
(outside the repo). It is visible from every CHPC resource, so one setup job
is sufficient for the whole cluster. An absent `renv/library/` entry inside
the repo checkout is normal and expected.

### Why PPM binaries are disabled

`RENV_CONFIG_PPM_ENABLED=FALSE` is set in `config.sh`. CHPC's R is built
against an external LAPACK and has no `libRlapack.so`. Posit Package Manager's
pre-built Linux binaries link against an R that does have it and fail to load
on CHPC with "libRlapack.so: cannot open shared object file". Disabling PPM
forces source installs that link against this cluster's R.

The renv cellar (`renv/cellar/fastcmprsk_1.24.10.tar.gz`) is vendored in the
repo and reaches CHPC on `git pull`, so renv uses it instead of fetching from
a remote.

### Before submitting: always source config.sh and load the R module

```bash
source sim/chpc/config.sh
module load "$R_MODULE"
```

Running R on the login node after sourcing config.sh is fine for short
interactive checks. Never run `renv::restore()` or `pkgbuild::compile_dll()`
on the login node — submit the setup job instead.

## Directory layout

### Scratch (per scenario)

```
$SCRATCH_BASE/
  n200.baseline/
    scenario.env           <- manifest: DGP fingerprint + creation timestamp + git SHA
    p25/
      p25_n200.baseline_run1.Rdata
      p25_n200.baseline_run2.Rdata
      ...
    p50/
      p50_n200.baseline_run1.Rdata
      ...
```

### Results (in the repo)

```
sim/results/
  n200.baseline/
    parsed/
      betas_p25_n200.baseline.csv
      model_performance_p25_n200.baseline.csv
      false_positives_p25_n200.baseline.csv
      by_block_p25_n200.baseline.csv
      betas_p50_n200.baseline.csv
      ...
      model_performance_n200.baseline.csv   <- stacked across all p (>1 p group)
```

## What happens end-to-end

1. `run_sim.sh` validates the scenario name, checks the scenario file exists,
   checks for collisions, then submits:
   - **setup job** (`bhcrr_setup_<scenario>`): restores renv and compiles
     `src/bhCRR.so`.
   - **run array** (`bhcrr_run_<scenario>`, `--dependency=afterok:setup`):
     splits `run_start..run_end` into `runs_per_task`-sized chunks. Each task
     writes `p<P>_<scenario>_run<N>.Rdata` into `$SCRATCH_BASE/<scenario>/p<P>/`.
   - **parse job** (`bhcrr_parse_<scenario>`, `--dependency=afterok:array`):
     reads all matching `.Rdata` files, runs consistency checks, and writes
     CSVs to `sim/results/<scenario>/parsed/`.

2. `skip_setup=1` skips the setup job and submits the array directly. Use only
   when `src/bhCRR.so` already exists and the renv library is current.

## Pulling results to local

```bash
scp -r '<you>@notchpeak.chpc.utah.edu:~/Documents/bhCRR/sim/results/n200.baseline/parsed' \
    ~/Documents/bhCRR/sim/results/n200.baseline/
```

## Monitoring & logs

```bash
squeue -u $USER
# for a specific cluster:
squeue -M notchpeak -u $USER
```

Per-task stdout/stderr land in `sim/chpc/logs/`:
- `setup_<job>.out`
- `run_<arrayjob>_<task>.out`
- `parse_<job>.out`

## Notes

- `set.seed()` is never called at the global level — each run derives a
  deterministic seed from `(scenario, p, run_id)` inside `run_one()`.
- Scratch is not backed up; `sim/results/` is the durable output.
- Array tasks use `load_sim("cluster")` which calls `pkgload::load_all(compile=FALSE)`
  — it never recompiles. The setup job compiles once; tasks inherit the `.so`.
