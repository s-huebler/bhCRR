---
editor_options:
  markdown:
    wrap: 72
---

# sim/timing/ — ClinTALL Fit-Time Scaling Study

How long does `fit_ssl_psdh` take at n = 1309 as the number of gene
predictors p grows?  This sub-study answers that with a SLURM array that
runs one timed fit per (p, rep) cell and collects the results into a
single CSV.

---

## Layout

```
sim/timing/
├── prep_clintall_timing.r     LOCAL only — materialise the input RDS + grid
├── fit_one_timing.r           one SLURM array task = one fit at (p, rep)
├── timing.slurm               SLURM job template (no resource flags)
├── submit_timing.sh           orchestrator: setup → array → collect
├── collect_timing.r           stack res_*.rds + fill missing from sacct
├── collect.slurm              SLURM wrapper for collect_timing.r
├── grid.tsv                   committed — p_index, p for every grid point
├── input_manifest.txt         committed — md5, dims, event counts, seed
├── input/                     gitignored — clintall_timing_input.rds
└── results/                   gitignored — per-run res_p*_rep*.rds + CSVs
```

---

## Workflow

### 1. Local: materialise the input

```r
# from the repo root
Rscript sim/timing/prep_clintall_timing.r
```

This reads `inst/ClinTall/dat_design_gexp.rds` and `dat_cr.rds`,
standardizes, draws the fixed column permutation (seed 20261030), writes
`sim/timing/input/clintall_timing_input.rds` (gitignored),
`sim/timing/grid.tsv` (committed), and
`sim/timing/input_manifest.txt` (committed).

Commit `grid.tsv` and `input_manifest.txt`, then push before submitting
to CHPC.

### 2. CHPC: pull and submit

```bash
ssh <you>@notchpeak.chpc.utah.edu
cd ~/Documents/bhCRR
git pull

# Copy the input file (gitignored, not in the repo)
scp ~/path/to/local/bhCRR/sim/timing/input/clintall_timing_input.rds \
    <you>@notchpeak.chpc.utah.edu:~/Documents/bhCRR/sim/timing/input/

source sim/chpc/config.sh   # optional; submit_timing.sh sources it
./sim/timing/submit_timing.sh tag=clintall-v1 reps=3 time=02:00:00 mem=16G
```

`submit_timing.sh` submits:
1. **setup job** (`bhcrr_setup_timing`): restores renv and compiles
   `src/bhCRR.so`.
2. **timing array** (`bhcrr_timing_<tag>`, `--dependency=afterok:setup`):
   one task per (p, rep) pair.  Indices map as:
   `i = (rep-1)*G + (p_index-1)`, so tasks `0..G-1` are rep 1.
3. **collect job** (`bhcrr_collect_<tag>`, `--dependency=afterany:array`):
   stacks all `res_*.rds` into `clintall_timing_<tag>.csv`.

Skip the setup job after the first run with `skip_setup=1`.

### 3. Monitor

```bash
squeue -M notchpeak -u $USER
```

Grep the greppable summary lines from task logs:
```bash
grep '^TIMING |' sim/chpc/logs/timing_*.out
```

### 4. Collect manually (if the collect job failed)

```bash
# On CHPC: write sacct output and run collect
sacct -M notchpeak -j <ARRAY_JOB_ID> \
  --format="JobID,JobName,State,Elapsed,Timelimit,MaxRSS,ExitCode" \
  --parsable2 > /scratch/general/vast/$USER/bhCRR/timing/clintall-v1/sacct_<id>.txt

TIMING_TAG=clintall-v1 \
TIMING_OUT_DIR=/scratch/general/vast/$USER/bhCRR/timing/clintall-v1 \
GRID_FILE=sim/timing/grid.tsv \
REPS=3 \
SACCT_FILE=/scratch/general/vast/$USER/bhCRR/timing/clintall-v1/sacct_<id>.txt \
  Rscript sim/timing/collect_timing.r
```

### 5. Pull results to local

```bash
scp '<you>@notchpeak.chpc.utah.edu:/scratch/general/vast/$USER/bhCRR/timing/clintall-v1/clintall_timing_clintall-v1.csv' \
    ~/Documents/bhCRR/sim/timing/results/
```

---

## Retry tiers

When tasks fail, resubmit with an explicit `array=` list of 0-based task
indices.  Use `skip_setup=1` to skip the compile step.

| Failure | Recommended resource bump | Example |
|---------|--------------------------|---------|
| `TIMEOUT` | Double the walltime; try the next tier | `time=06:00:00` → `time=24:00:00` |
| `OUT_OF_MEM` | Double the memory | `mem=16G` → `mem=32G` |
| Other error | Investigate the log; retry at same resources | — |

**Retry example** (tasks 3 and 7 timed out at 1 h):

```bash
./sim/timing/submit_timing.sh \
    tag=clintall-v1 \
    array=3,7 \
    time=06:00:00 \
    skip_setup=1
```

The collect script keeps only the latest successful row per (p, rep), so
retry outputs safely overlay the original failed rows.

---

## P-grid and array geometry

`grid.tsv` contains one row per p value, built by
`make_timing_grid(pmax, start = 500, step = 2000)`:

```
p ∈ {500, 2500, 4500, …, 22500, 24500, 24619}   (G = 14)
```

The current grid uses `step = 2000`.  To regenerate `grid.tsv` after
changing the step **without** touching `input/` or `input_manifest.txt`:

```bash
Rscript -e "
  p <- sort(unique(c(seq(500L, 24619L, by = 2000L), 24619L)))
  d <- data.frame(p_index = seq_along(p), p = p)
  write.table(d, 'sim/timing/grid.tsv', sep='\t', quote=FALSE, row.names=FALSE)
  cat('G=', length(p), '\n')
"
```

Do **not** source `prep_clintall_timing.r` for this — sourcing runs the
full script and overwrites `input/clintall_timing_input.rds`.

Array task `i` (0-based) maps to:

```
p   = grid$p[ (i %% G) + 1 ]
rep = i %/% G + 1
```

With G = 14 and reps = 3 the full array is `0..41` (42 tasks).

---

## What is NOT done here

- No cross-validation, no autotuning (`bhcrr_autotune`).
- No comparator models (LASSO, aLASSO, SCAD, MCP).
- The ClinTALL `.rds` files are **never committed**; they are gitignored
  under `inst/ClinTall/`.
