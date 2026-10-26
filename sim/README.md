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
│   └── scenarioA.R      baseline scenario (200 obs, block-correlated Gaussian)
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
│   ├── run_array.slurm  array job: calls sim/chpc/R/run_batch.R
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

Each scenario lives in `sim/scenarios/<name>.R` and returns a plain list named
`scenario`. The list is validated by `scenario_from_file()` in `sim/R/scenario.r`.

**Changing any DGP parameter requires a new scenario name** so that old and new
results are never mixed in `sim/results/`.

---

## Local workflow

### Single interactive run

```r
# From R (project root or any working directory):
SCENARIO_FILE <- "sim/scenarios/scenarioA.R"
P             <- 25
RUN_ID        <- 1
source("sim/local/single.R")
# `result` is now in the global environment
str(result$scores)
```

### Local batch sweep

```r
SCENARIO_FILE <- "sim/scenarios/scenarioA.R"
NPREDICTORS   <- c(25, 50)
RUN_START     <- 1
RUN_END       <- 10
RESULTS_DIR   <- "sim/results/"
source("sim/local/batch.R")
# Writes sim/results/scenarioA/runs/result_p25_run01.Rdata etc.
```

---

## CHPC workflow

```bash
ssh <you>@notchpeak.chpc.utah.edu
cd ~/bhCRR
source sim/chpc/config.sh          # optional; run_sim.sh sources it too

# Submit a 10-run array for scenarioA, p = 25 and 50:
./sim/chpc/run_sim.sh \
    scenario=scenarioA \
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

## Migration status

| Legacy file | Status | Replaced by |
|---|---|---|
| `_legacy/code/single_run_complex.R` | archived | `sim/local/single.R` + `sim/R/run_one.r` |
| `_legacy/code/single_run.R` | archived | same |
| `_legacy/code/batch_run_complex.R` | archived | `sim/local/batch.R` + `sim/R/run_one.r` |
| `_legacy/code/batch_run.R` | archived | same |
| `_legacy/code/batch_run_wolbers.R` | archived | `sim/R/score.r` (metrics consolidated) |
| `_legacy/code/batch_run_cv_wolbers.R` | archived | same |
| `_legacy/code/batch_run_comparison.R` | archived | `sim/R/comparators.r` |
| `_legacy/code/single_run_comparison.R` | archived | same |
| `_legacy/code/other_comparators_cv_section.R` | archived | `sim/R/comparators.r` |
| `_legacy/code/parse_batch_run.R` | archived | `sim/R/parse.r` |
| `chpc/R/parse_batch.R` | thin wrapper | calls `sim/R/parse.r` |
| `chpc/R/run_batch.R` | thin wrapper | calls `sim/R/run_one.r` |

All `sim/R/` files are currently **header/contract stubs** — implementations are
not yet written. The stubs define the public API and field contracts so that
`run_one.r` → `parse.r` dependencies can be understood before any code is
written.
