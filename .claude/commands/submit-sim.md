---
description: Submit a bhCRR simulation job array on CHPC
argument-hint: <n runs>, <only the settings that differ from defaults>, <scenario name>
---

Submit a bhCRR simulation on CHPC.

REQUEST: $ARGUMENTS

Work in ~/bhCRR. Follow these rules exactly.

## How to read the request

The request names ONLY the settings that differ from the defaults in
`chpc/config.sh`, plus the run count and the scenario name. Everything not
mentioned is a default and must NOT appear as an override on the command line —
do not restate a default, and never invent a value that was not given.
"All defaults" or "the default scenario" means the only overrides are the run
range and the scenario name.

Common phrasings map to `key=value` overrides like this:

| said | override |
|------|----------|
| p, predictors, npredictors | `npredictors=` |
| n, sample size, nobs | `nobs=` |
| beta1 / beta2 (active coefficients) | `beta1_active=` / `beta2_active=` |
| block correlations | `block_rho=` (`name=val,name=val`) |
| block proportions / labels | `block_props=` / `active_block=` |
| censoring rate, observation window | `cens_rate=` / `u_min=` / `u_max=` |
| latent settings | `latent_type=` / `latent_q=` |
| walltime | `sb_time=` |
| runs per array task | `runs_per_task=` |

`key=value` takes no spaces around the `=` and needs no quotes. The parser splits
on the FIRST `=`, so named-vector values like
`block_rho=indep=0.0,weak=0.5,strong=0.8` are safe as written.

If any part of the request is ambiguous, or you find yourself guessing at a
value, STOP and ask. Do not submit a job you had to interpret.

## Preflight

All of this is cheap shell — no R. See the Arbiter rule at the bottom.

1. `git pull`, and report the resulting HEAD short SHA.
2. Confirm `chpc/run_sim.sh` mentions `SCENARIO`. If it does not, the
   scenario-aware pipeline was never pushed to this checkout — stop and say so.
3. Confirm the renv library exists and is populated:
   `ls renv/library/R-4.5/x86_64-pc-linux-gnu | head` — and that `fastcmprsk`
   and `tidyverse` are in it. `USE_RENV` defaults to 0, so a missing or thin
   library means every array task fails. Do NOT try to fix it. Stop and say that
   `renv::restore()` needs to run inside an `salloc` session first.
4. Validate the scenario slug against `[A-Za-z0-9][A-Za-z0-9.-]*` — no
   underscores. If it doesn't match, stop and ask; never silently rewrite it.
5. Pick the run range. If the request gives an explicit start, use it. Otherwise
   list `$SCRATCH_BASE/<scenario>/p<P>/` for each requested `npredictors`, take
   the highest existing `run<N>`, and start at N+1 so this batch appends instead
   of colliding. Empty or missing directory means start at 1. State the range you
   chose in your report.
   An existing scenario directory is normal and expected — scenarios accumulate
   runs across submissions. It is not a reason to stop.

## Submit

Run `./chpc/run_sim.sh` once, for real, with the overrides you built. No dry run.

Stop and report instead of retrying if:

- the collision guard fires on existing run indices — report which, and what
  `run_start` would append cleanly
- `run_sim.sh` hard-fails on a `scenario.env` mismatch — show me the diff it
  printed. NEVER resolve this by picking a different scenario name, passing
  `overwrite=1`, or deleting `scenario.env`. A mismatch means this scenario name
  already stands for different settings, and only I can decide which is wrong.
- sbatch itself errors

## Report back

- The exact command you ran, and a one-line list of what it overrode
  (or "defaults only" when it overrode nothing but the run range and scenario)
- The run range used, and why, if it wasn't given explicitly
- Array job ID + task count, and the parse job ID
- Scratch path, and the `Sim/chpc_results/<scenario>/` files the parse job will
  write
- The `squeue` command to track it, with `-M $CHPC_CLUSTER` if that is set

Run `squeue` ONCE to confirm the jobs are queued. Do not poll, wait, or watch.

## Arbiter rule (non-negotiable)

You are on a login node governed by Arbiter: roughly 15 core-minutes at or under
4 GB before penalties escalate. Everything you do here is `git`, `ls`/`cat`, and
`sbatch`. Do NOT run `Rscript`, `R`, `renv::restore()`, `Rcpp::sourceCpp`,
`chpc/R/run_batch.R`, or `chpc/R/parse_batch.R` on this node — not even "just to
check". Compiling and model fitting belong in the array; a re-parse belongs in an
`salloc` job.
