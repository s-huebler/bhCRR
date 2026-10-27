# chpc-timing — Parked timing harness

Parked 2026-10-27.

## Why parked

These files make up a separate EM-timing study (not the main simulation
pipeline). They depend on archived R/ files that were never ported:

- `R/dedupe_warnings.r` — folded into `R/utils.r` during the sim/ migration
- `R/generate_foldid.r` — archived; functionality superseded by `bhcrr_make_folds`
- `R/threshold.R` — archived; never used in the rewritten package

`run_timing.R` and `prep_timing_data.R` source those paths directly and
cannot run as-is.

## What reviving it would need

1. Replace the three `source()` calls in `run_timing.R` with `load_sim("cluster")`
   (same treatment the run and parse wrappers got in step 12).
2. Wire `timing_config.sh` and `timing.slurm` to use the same `USE_RENV` knob
   and `setup.slurm` pre-job that step 12b adds to the main pipeline.
   See `sim/chpc/README.md`'s renv section for the three-state `USE_RENV` design
   rather than duplicating it here.
3. Confirm that `prep_timing_data.R`'s ClinTALL data paths still exist under
   `inst/ClinTall/` after any subsequent repo reorganization.

## Files

| File | Role |
|------|------|
| `run_timing.R` | Fits ONE bhCRR model from a pre-materialized timing config `.rds`; writes a result `.rds`. No autotuning, no comparators. |
| `prep_timing_data.R` | Materializes nested gene subsets of ClinTALL data for timing experiments; writes one `.rds` per config. |
| `run_timing.sh` | Submission script for the timing array. |
| `timing_config.sh` | Config for the timing study (separate from `config.sh`; carries its own `USE_RENV:=0`). |
| `timing.slurm` | SLURM job template for timing runs; includes its own `renv::restore()` block. |
| `prep_timing.slurm` | sbatch wrapper that runs `prep_timing_data.R` (resource flags on the sbatch line). Never tracked before 2026-10-27; existed only in the CHPC checkout. |
| `timing_measurements_2026-10-27.txt` | The 36 `TIMING \|` lines (fit time, BIC-init time, peak memory) grepped from the 12 gitignored `timing_*.out` logs: the study's only empirical scaling data. |
