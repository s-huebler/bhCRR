---
editor_options: 
  markdown: 
    wrap: 72
---

# sim/ — bhCRR Simulation Study

Simulation infrastructure for benchmarking bhCRR against penalised
Fine-Gray comparators (LASSO, aLASSO, SCAD, MCP) under block-correlated
covariates and competing risks.

------------------------------------------------------------------------

## Layout

```         
sim/
├── R/                   shared core: load, DGP, fitting, scoring, parsing
│   ├── load.r           bring bhCRR into scope (single source-of-truth)
│   ├── dgp.r            data-generating process: covariates + outcomes
│   ├── scenario.r       scenario constructor, file reader, manifest I/O
│   ├── fit.r            fit_bhcrr() wrapper (calls bhcrr_autotune)
│   ├── comparators.r    LASSO / aLASSO / SCAD / MCP fits
│   ├── score.r          all evaluation metrics for one run
│   ├── run_one.r        single unit of work called by both batch drivers
│   └── parse.r          many .Rdata run files → summary CSVs
│
├── local/               interactive drivers (function definitions, no SLURM)
│   ├── single.R         defines run_single(): one run, printed summary
│   └── batch.R          defines batch_sweep(): sweep over (p, run_id)
│
├── scenarios/           scenario definitions (plain R lists)
│   ├── _template.R      master template; copy with new_scenario(), never run directly
│   └── <name>.R         real scenarios created from the template
│
├── results/             run outputs (git-ignored except structure)
│   ├── <scenario>/
│   │   ├── scenario.env  manifest written by batch_sweep() / run_sim.sh
│   │   ├── runs/         one p<P>_<scenario>_run<N>.Rdata per run_one() call
│   │   └── parsed/       summary CSVs + scenario.env written by parse_results()
│   └── _legacy/         pre-restructure Local_Testing outputs (git-ignored)
│
├── chpc/                CHPC / SLURM infrastructure
│   ├── config.sh        the one file to source before submitting
│   ├── setup.slurm      one-shot prepare job (renv + compile_dll)
│   ├── run_sim.sh       orchestrator: three-job chain (setup → array → parse)
│   ├── run_array.slurm  array job: calls sim/chpc/R/run_array.r
│   ├── parse.slurm      post-array parse job
│   ├── R/               thin CHPC wrappers (source load.r, call sim/R/)
│   └── logs/            SLURM stdout/stderr (git-ignored)
│
└── _legacy/
    └── code/            archived pre-restructure scripts (read-only reference)
```

------------------------------------------------------------------------

## Scenarios

A **scenario** is a named set of DGP parameters. All parameters that
affect the data-generating process — `nobs`, coefficients, block
structure, censoring — are part of the scenario. Parameters that only
control how you sweep (number of predictors, run range) are not part of
the scenario.

**Changing any DGP parameter requires a new scenario name** so that old
and new results are never mixed in `sim/results/`.

### Creating a new scenario

``` r
source("sim/R/load.r")
load_sim()
new_scenario("myScenario.01", dir = "sim/scenarios")
# Opens sim/scenarios/myScenario.01.R with the name field pre-filled.
# Fill in every value marked [PLACEHOLDER] before running.
```

`new_scenario()` copies `sim/scenarios/_template.R` and fills in the
`name` field. `_template.R` carries inline documentation for every key
and uses deliberately small placeholder values that run end-to-end as a
smoke test.

**Name rule:** scenario names must match `^[A-Za-z0-9][A-Za-z0-9.-]*$`.
No underscores — they are field delimiters in the
`p<P>_<scenario>_run<N>` filename format. Use `.` or `-` as word
separators.

Each scenario file assigns exactly one object named `scenario` (a plain
list) and must contain no calls to `library()`, `require()`, `source()`,
or `set.seed()`. `scenario_from_file()` in `sim/R/scenario.r` enforces
this before any validation runs.

------------------------------------------------------------------------

## Local workflow

Both drivers in `sim/local/` are **functions**, not source-and-assign
scripts. Sourcing the file defines the function; you then call it. (The
old `SCENARIO_FILE <- ...; source("sim/local/single.R")` pattern was
retired because a long debugging session could feed a stale variable
into a new run.)

### Single interactive run

``` r
source("sim/local/single.R")

result <- run_single(
  scenario_file = "sim/scenarios/myScenario.01.R",
  p             = 25,
  run_id        = 1
)
str(result$scores)
```

`run_single()` prints a compact summary (bhCRR tuning, score tables,
timing, warning counts) and returns the full `run_one()` result
**invisibly** — assign it, or the returned object is lost.
`scenario_file` is required with no default. Warnings from `fastcmprsk`
and bhCRR are collected into `$warnings` rather than flooding the
console. `quiet = TRUE` suppresses all printed output;
`show_by_block = FALSE` omits the per-block score table from the summary
but keeps it in `$scores$by_block`.

**Quick look with a cheap control** (avoids the 10-fold × 2-rep
production cost — slow at larger p):

``` r
result <- run_single(
  "sim/scenarios/myScenario.01.R",
  p           = 25,
  run_id      = 1,
  fit_control = bhcrr_cv_control(nfolds = 3, ncv = 1),
  s0_seq      = c(0.04, 0.08),
  s1_seq      = c(0.4,  0.8),
  max_widen   = 0L
)
```

To truncate the fold fits further, add `fit_args = list(maxit = 20L)` to
the control. The default control is production settings;
`nfolds = 3, ncv = 1` trades tuning quality for speed during
development.

### Local batch sweep

``` r
source("sim/local/batch.R")

batch_sweep(
  scenario_file = "sim/scenarios/myScenario.01.R",
  npredictors   = c(25, 50),
  run_start     = 1,
  run_end       = 10
)
# Writes sim/results/myScenario.01/runs/p25/p25_myScenario.01_run1.Rdata etc.
```

`batch_sweep()` writes one `.Rdata` file per `(p, run_id)` pair and logs
progress to the console. Key behaviours:

-   **Collision guard / overwrite**: by default, existing files are
    skipped. Pass `overwrite = TRUE` to replace them.
-   **Failure records**: a run that errors saves a stub with
    `$failed = TRUE`, a reproducible seed, and the error message.
    `parse_results()` reports and skips failure stubs.
-   **`dryrun = TRUE`**: prints the file paths that *would* be written
    without running any fits.
-   `fit_control`, `s0_seq`, `s1_seq`, `max_widen`, and `seed` pass
    through to every `run_one()` call.

### Parsing results

``` r
source("sim/R/load.r")
load_sim()

parse_results(
  runs_root = "sim/results/myScenario.01/runs",
  scenario  = "myScenario.01"
)
# Writes CSVs to sim/results/myScenario.01/parsed/
```

`parse_results()` locates `scenario.env` automatically in both layouts:
local (`<scenario>/runs/` with the manifest at `<scenario>/`), and CHPC
(runs and manifest together directly in `$SCRATCH_RUN_DIR`). This
auto-search exists because the single-layout assumption was a defect in
an earlier implementation.

------------------------------------------------------------------------

## CHPC workflow

``` bash
ssh <you>@notchpeak.chpc.utah.edu
cd ~/bhCRR
source sim/chpc/config.sh          # optional; run_sim.sh sources it too

# Submit a 10-run array for myScenario.01, p = 25 and 50:
./sim/chpc/run_sim.sh \
    scenario=myScenario.01 \
    npredictors=25,50 \
    run_start=1 run_end=10
```

The **scenario file is the only source of DGP parameters** — every
coefficient, block structure, and censoring rate lives in
`sim/scenarios/<slug>.R`. `run_sim.sh` exits with an error on a DGP key
(e.g. `nobs=200`) and on any unknown key.

`run_sim.sh` accepts sweep and resource keys as `key=value` with no
spaces around the equals sign:

| Key                    | What it controls                           |
|------------------------|--------------------------------------------|
| `scenario`             | Scenario slug (required)                   |
| `npredictors`          | Comma-separated p values                   |
| `run_start`, `run_end` | Run index range                            |
| `runs_per_task`        | Runs bundled per array task                |
| `SB_*`                 | SLURM resource overrides (time, mem, cpus) |
| `dryrun=1`             | Print `sbatch` commands without submitting |

See `sim/chpc/README.md` for the full list.

------------------------------------------------------------------------

## Results layout

```         
sim/results/<scenario>/
├── scenario.env              scenario manifest written by batch_sweep()
├── runs/
│   ├── p25/
│   │   ├── p25_<scenario>_run1.Rdata
│   │   └── p25_<scenario>_run2.Rdata
│   └── p50/
│       └── ...
└── parsed/
    ├── scenario.env              copy written by parse_results()
    ├── betas_p25_<scenario>.csv
    ├── model_performance_p25_<scenario>.csv
    ├── false_positives_p25_<scenario>.csv
    └── model_performance_<scenario>.csv    # stacked across all p values
```

`scenario.env` exists in two places: written at `<scenario>/` by
`batch_sweep()` / `run_sim.sh` on the first run, and copied into
`parsed/` by `parse_results()` so the parsed output is self-contained.

The `p<P>/` level is load-bearing: the column-to-block map depends on
`p`, so betas, by_block, and false_positives are never stacked across
different p values.

Each `p<P>_<scenario>_run<N>.Rdata` contains a single object `result`
with fields: `$run_id`, `$p`, `$scenario`, `$data`, `$fit`,
`$comparators`, `$scores`, `$meta`. See `sim/R/run_one.r` for the exact
contract.

------------------------------------------------------------------------

## Testing convention

`sim/` is excluded from the package build by `^sim$` in `.Rbuildignore`.
`R CMD check` runs tests against the built package, where `sim/` does
not exist — so **sim tests must NOT live in `tests/testthat/`**.

Simulation tests live in `sim/tests/test_<module>.R`, each runnable
standalone:

``` r
testthat::test_file("sim/tests/test_dgp.R")
```

`sim/tests/run_all.R` runs them in dependency order. They are never
wired into `devtools::test()`.

**Test-only dependency:** `withr` is used for environment variable
isolation in the tests (`withr::with_envvar`). It is already in
`renv.lock` (pulled in by `devtools`), but it is not declared in the
package `DESCRIPTION` since it is only used here, not in the package
itself.

------------------------------------------------------------------------

## Migration status

Steps in implementation order; mark done as each lands.

| Step | Module | Status | Notes |
|----|----|----|----|
| 1 | contracts | **done** (0541e8e → 4c2f973) | headers only; corrected before migration |
| 2 | `sim/R/load.r` | **done** (71b345e) |  |
| 3 | `sim/R/scenario.r` + `sim/scenarios/_template.R` | **done** (5d9bce0) |  |
| 4 | `sim/R/dgp.r` | **done** (f808398) | block_signature extended with latent/lbal/uvar (PART C) |
| 5 | `sim/R/fit.r` | **done** (bcb8b15) | eval_time fixed at scenario midpoint; fg_oracle separated |
| 6 | `sim/R/comparators.r` — fits | **done** (this commit) | BIC selection; aLASSO uses unpenalized FG initial |
| 7 | revisit comparator tuning | parked | see comparators.r FUTURE WORK; legacy CV at `sim/_legacy/code/other_comparators_cv_section.R` |
| 8 | `sim/R/score.r` | **done** (821bb48) | per-block stratification; fg_oracle uniform shape |
| 9 | `sim/R/run_one.r` | **done** (20d9e6c) | deterministic seeding; fg_oracle separate field |
| 10 | `sim/local/single.R` + `sim/local/batch.R` | **done** | run_single(); batch_sweep(); failure records; collision guard |
| 11 | `sim/R/parse.r` | **done** (this commit) | pure stacker; fg_oracle as model row; by_block CSV added |
| 12 | `sim/chpc/R/run_array.r` + `sim/chpc/R/parse.r` wrappers + shell wiring | **done** | R wrappers thin; setup.slurm one-shot prepare; DGP moved to scenario file; timing harness parked in `sim/_legacy/chpc-timing/` |
| 13 | retire `sim/_legacy/` | not started | after end-to-end run confirmed |

Legacy files absorbed by each step:

| Legacy file | Absorbed at step | Replaced by |
|----|----|----|
| `_legacy/code/single_run_complex.R` | 9 + 10 | `sim/R/run_one.r` + `sim/local/single.R` |
| `_legacy/code/single_run.R` | 9 + 10 | same |
| `_legacy/code/batch_run_complex.R` | 9 + 10 | `sim/R/run_one.r` + `sim/local/batch.R` |
| `_legacy/code/batch_run.R` | 9 + 10 | same |
| `_legacy/code/batch_run_wolbers.R` | 8 | `sim/R/score.r` (metrics consolidated) |
| `_legacy/code/batch_run_cv_wolbers.R` | 8 | same |
| `_legacy/code/batch_run_comparison.R` | 6 + 7 | `sim/R/comparators.r` |
| `_legacy/code/single_run_comparison.R` | 6 + 7 | same |
| `_legacy/code/other_comparators_cv_section.R` | 6 + 7 | same |
| `_legacy/code/parse_batch_run.R` | 11 | `sim/R/parse.r` |
| `sim/chpc/R/parse_batch.R` | 12 | thin wrapper calling `sim/R/parse.r` |
| `sim/chpc/R/run_batch.R` | 12 | renamed to `sim/chpc/R/run_array.r` |
