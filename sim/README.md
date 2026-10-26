# sim/ — bhCRR Simulation Study

Simulation infrastructure for benchmarking bhCRR against penalised Fine-Gray
comparators (LASSO, aLASSO, SCAD, MCP) under block-correlated covariates and
competing risks.

---

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
├── local/               interactive drivers (no SLURM)
│   ├── single.R         one run; leaves `result` in the global environment
│   └── batch.R          sweep over (p, run_id); writes .Rdata per run
│
├── scenarios/           scenario definitions (plain R lists)
│   ├── _template.R      master template; copy with new_scenario(), never run directly
│   └── <name>.R         real scenarios created from the template
│
├── results/             run outputs (git-ignored except structure)
│   ├── <scenario>/
│   │   ├── runs/        one result_p<p>_run<id>.Rdata per run_one() call
│   │   └── parsed/      summary CSVs written by parse_results()
│   └── _legacy/         pre-restructure Local_Testing outputs (git-ignored)
│
├── chpc/                CHPC / SLURM infrastructure
│   ├── config.sh        the one file to source before submitting
│   ├── run_sim.sh       orchestrator: builds the job array
│   ├── run_array.slurm  array job: calls sim/chpc/R/run_array.r
│   ├── parse.slurm      post-array parse job
│   ├── R/               thin CHPC wrappers (source load.r, call sim/R/)
│   └── logs/            SLURM stdout/stderr (git-ignored)
│
└── _legacy/
    └── code/            archived pre-restructure scripts (read-only reference)
```

---

## Scenarios

A **scenario** is a named set of DGP parameters. All parameters that affect the
data-generating process — `nobs`, coefficients, block structure, censoring — are
part of the scenario. Parameters that only control how you sweep (number of
predictors, run range) are not part of the scenario.

**Changing any DGP parameter requires a new scenario name** so that old and new
results are never mixed in `sim/results/`.

### Creating a new scenario

```r
source("sim/R/load.r")
load_sim()
new_scenario("myScenario.01", dir = "sim/scenarios")
# Opens sim/scenarios/myScenario.01.R with the name field pre-filled.
# Fill in every value marked [PLACEHOLDER] before running.
```

`new_scenario()` copies `sim/scenarios/_template.R` and fills in the `name`
field. `_template.R` carries inline documentation for every key and uses
deliberately small placeholder values that run end-to-end as a smoke test.

**Name rule:** scenario names must match `^[A-Za-z0-9][A-Za-z0-9.-]*$`. No
underscores — they are field delimiters in the `p<P>_<scenario>_run<N>` filename
format. Use `.` or `-` as word separators.

Each scenario file assigns exactly one object named `scenario` (a plain list)
and must contain no calls to `library()`, `require()`, `source()`, or
`set.seed()`. `scenario_from_file()` in `sim/R/scenario.r` enforces this before
any validation runs.

---

## Local workflow

### Single interactive run

```r
# From R (project root or any working directory):
SCENARIO_FILE <- "sim/scenarios/myScenario.01.R"
P             <- 25
RUN_ID        <- 1
source("sim/local/single.R")
# `result` is now in the global environment
str(result$scores)
```

### Local batch sweep

```r
SCENARIO_FILE <- "sim/scenarios/myScenario.01.R"
NPREDICTORS   <- c(25, 50)
RUN_START     <- 1
RUN_END       <- 10
RESULTS_DIR   <- "sim/results/"
source("sim/local/batch.R")
# Writes sim/results/myScenario.01/runs/result_p25_run01.Rdata etc.
```

---

## CHPC workflow

```bash
ssh <you>@notchpeak.chpc.utah.edu
cd ~/bhCRR
source sim/chpc/config.sh          # optional; run_sim.sh sources it too

# Submit a 10-run array for myScenario.01, p = 25 and 50:
./sim/chpc/run_sim.sh \
    scenario=myScenario.01 \
    npredictors=25,50 \
    run_start=1 run_end=10
```

All DGP parameters that differ from `config.sh` defaults can be passed as
`key=value` overrides on the `run_sim.sh` command line. See
`sim/chpc/README.md` for full options.

---

## Results layout

```
sim/results/<scenario>/
├── runs/
│   ├── result_p25_run001.Rdata
│   ├── result_p25_run002.Rdata
│   └── ...
└── parsed/
    ├── betas_p25_scenarioA.csv
    ├── model_performance_p25_scenarioA.csv
    ├── false_positives_p25_scenarioA.csv
    └── model_performance_scenarioA.csv    # stacked across all p values
```

Each `result_p<p>_run<id>.Rdata` contains a single object `result` with fields:
`$run_id`, `$p`, `$scenario`, `$data`, `$fit`, `$comparators`, `$scores`,
`$meta`. See `sim/R/run_one.r` for the exact contract.

---

## Testing convention

`sim/` is excluded from the package build by `^sim$` in `.Rbuildignore`.
`R CMD check` runs tests against the built package, where `sim/` does not
exist — so **sim tests must NOT live in `tests/testthat/`**.

Simulation tests live in `sim/tests/test_<module>.R`, each runnable standalone:

```r
testthat::test_file("sim/tests/test_dgp.R")
```

`sim/tests/run_all.R` runs them in dependency order. They are never wired into
`devtools::test()`.

**Test-only dependency:** `withr` is used for environment variable isolation in
the tests (`withr::with_envvar`). It is already in `renv.lock` (pulled in by
`devtools`), but it is not declared in the package `DESCRIPTION` since it is only
used here, not in the package itself.

---

## Migration status

Steps in implementation order; mark done as each lands.

| Step | Module | Status | Notes |
|------|--------|--------|-------|
| 1 | contracts | **done** (0541e8e → 4c2f973) | headers only; corrected before migration |
| 2 | `sim/R/load.r` | **done** (71b345e) | |
| 3 | `sim/R/scenario.r` + `sim/scenarios/_template.R` | **done** (5d9bce0) | |
| 4 | `sim/R/dgp.r` | **done** (f808398) | block_signature extended with latent/lbal/uvar (PART C) |
| 5 | `sim/R/fit.r` | **done** (bcb8b15) | eval_time fixed at scenario midpoint; fg_oracle separated |
| 6 | `sim/R/comparators.r` — fits | **done** (this commit) | BIC selection; aLASSO uses unpenalized FG initial |
| 7 | revisit comparator tuning | parked | see comparators.r FUTURE WORK; legacy CV at `sim/_legacy/code/other_comparators_cv_section.R` |
| 8 | `sim/R/score.r` | **done** (821bb48) | per-block stratification; fg_oracle uniform shape |
| 9 | `sim/R/run_one.r` | **done** (this commit) | deterministic seeding; fg_oracle separate field |
| 10 | `sim/local/single.R` + `sim/local/batch.R` | not started | |
| 11 | `sim/R/parse.r` | not started | |
| 12 | `sim/chpc/R/run_array.r` + `sim/chpc/R/parse_batch.R` wrappers | not started | thin wrappers only |
| 13 | retire `sim/_legacy/` | not started | after end-to-end run confirmed |

Legacy files absorbed by each step:

| Legacy file | Absorbed at step | Replaced by |
|---|---|---|
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
