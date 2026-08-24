# bhCRR on CHPC (SLURM)

Run the `batch_run_complex` simulation as a parallel SLURM job array on the
University of Utah CHPC, then automatically collect the per-parameter-set results
(`betas_p<P>`) into `Sim/chpc_results/` to pull back to local for figures.

## Files

| File | Role |
|------|------|
| `config.sh` | The **one** file you edit/source. Paths, R module, SLURM account, and default sim params. |
| `run_sim.sh` | Entry point. Parses `key=value` overrides, submits the array + a dependent parse job. |
| `run_array.slurm` | Array task: runs a chunk of run indices, writes one `.Rdata` per run to scratch. |
| `parse.slurm` | Runs once after the array succeeds; builds `betas_p<P>`. |
| `R/run_batch.R` | CHPC-portable version of `Sim/code/batch_run_complex.R` (env-driven, renv). |
| `R/parse_batch.R` | CHPC-portable version of `Sim/code/parse_batch_run.R`; writes CSV + `.Rdata`. |

## One-time setup

1. Clone/pull the repo on CHPC and note its path (e.g. `~/bhCRR`).
2. Allocation is prefilled in `chpc/config.sh` for the `qiaox` / lonepeak setup
   (carried over from the ODSiData CHPC config):
   - `CHPC_ACCOUNT` — defaults to `qiaox`
   - `CHPC_PARTITION` — defaults to `lonepeak-shared` (shared, queues faster than
     a whole node). Override examples are in the config comments.
   - `CHPC_CLUSTER` — defaults to `lonepeak`; each CHPC cluster runs its own Slurm
     controller, so this routes `sbatch --clusters=…`. Set to `""` for your login
     cluster's scheduler.
   The one thing to confirm once is:
   - `R_MODULE` — check with `module spider R` (renv.lock pins **R 4.5.1**)
   - `SCRATCH_BASE` — defaults to `/scratch/general/vast/$USER/bhCRR`; adjust if needed
3. Make sure `renv` can restore: the first array task runs `renv::restore()`.

You can override any of these per-submission without editing the file, e.g.
`CHPC_PARTITION=notchpeak-shared-short CHPC_CLUSTER=notchpeak ./chpc/run_sim.sh`.

## Everyday use

```bash
ssh <you>@notchpeak.chpc.utah.edu
cd ~/bhCRR
source chpc/config.sh          # optional — run_sim.sh sources it too
./chpc/run_sim.sh              # defaults: nobs=200, npredictors=221, runs 1..10
```

Override any default on the command line (`key=value`, **no spaces** around `=`):

```bash
./chpc/run_sim.sh run_start=10 run_end=40
./chpc/run_sim.sh nobs=140
./chpc/run_sim.sh npredictors=25,50,221 run_end=100 runs_per_task=5
./chpc/run_sim.sh nobs=140 dryrun=1     # print sbatch commands, submit nothing
```

`dryrun=1` is the way to confirm all scenario knobs (coefficients, block
structure, censoring) before a real submission — the banner prints every resolved
value. Off-cluster (Mac), point `SCRATCH_BASE` at a local writable path to avoid
the `/scratch` mkdir failing:

```bash
SCRATCH_BASE=/tmp/bhcrr_dryrun ./chpc/run_sim.sh dryrun=1
```

Named-vector keys (`block_props`, `block_rho`, `latent_q`) use `name=val,name=val`
format. The `run_sim.sh` parser splits only on the **first** `=`, so the extra `=`
inside the value is safe:

```bash
./chpc/run_sim.sh block_rho=indep=0.0,weak=0.5,strong=0.8
```

Keys map to the `UPPER_CASE` variables in `config.sh` (case-insensitive):
`nobs, npredictors, run_start, run_end, runs_per_task, zero_gap_target,
beta1_active, beta2_active, active_block, block_props, block_rho, latent_type,
latent_q, cens_rate, u_min, u_max`, plus SLURM knobs `sb_time, sb_mem, sb_cpus`.
The allocation vars `CHPC_ACCOUNT, CHPC_PARTITION, CHPC_CLUSTER` are also honored
if exported (e.g. `CHPC_PARTITION=notchpeak-shared-short ./chpc/run_sim.sh`).

## What happens

1. `run_sim.sh` splits `run_start..run_end` into `runs_per_task`-sized chunks and
   submits a **job array** (`--array=0-N`). Each task writes
   `n<nobs>_p<p>_run<run>.Rdata` to `$SCRATCH_BASE/runs_n<nobs>/`.
2. A **parse job** is submitted with `--dependency=afterok` on the array. When
   every task succeeds it groups all `.Rdata` by `npredictors` and writes, for
   each group, to `Sim/chpc_results/`:
   - `betas_p<P>_n<NOBS>.csv`
   - `betas_p<P>_n<NOBS>.Rdata` (loads an object named `betas_p<P>`)

Runs with the same `nobs` accumulate in the same scratch dir, so
`run_start=10 run_end=40` **adds** to an earlier `1..10` batch and the parse picks
up everything present.

## Pulling results to local

From your Mac:

```bash
scp '<you>@notchpeak.chpc.utah.edu:~/bhCRR/Sim/chpc_results/betas_p*_n*.csv' \
    ~/Documents/bhCRR/Sim/chpc_results/
```

## Monitoring & logs

```bash
squeue -u $USER
```

Per-task stdout/stderr land in `chpc/logs/run_<arrayjob>_<task>.out` and
`chpc/logs/parse_<job>.out`.

## Notes

- `set.seed()` is intentionally never called — each run draws fresh data.
- Scratch is not backed up; treat `Sim/chpc_results/` as the durable output.
- Each array task recompiles the Rcpp sources via `sourceCpp`. If that becomes a
  bottleneck, precompiling `src/*.so` once is a future optimization.
