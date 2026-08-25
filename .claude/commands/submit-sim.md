---
description: Submit a bhCRR simulation job array on CHPC
argument-hint: <n runs>, <only the settings that differ from defaults>, <scenario name>
---

Submit a bhCRR simulation on CHPC.

REQUEST: $ARGUMENTS

## Where you are

Work in the CURRENT directory. Do not `cd` to a hardcoded path — this repo is
not necessarily at `~/bhCRR`. Confirm you are in the right place with
`git rev-parse --show-toplevel` and report that path. If it is not a git repo,
or `chpc/run_sim.sh` is not present under it, stop and say so.

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

Each Bash call is a FRESH shell — nothing persists between them. Any command that
needs `$SCRATCH_BASE`, `$R_MODULE`, `$CHPC_CLUSTER` or another config value must
`source chpc/config.sh` in that same command.

1. `git pull`, and report the resulting HEAD short SHA.
2. Confirm `chpc/run_sim.sh` mentions `SCENARIO`. If it does not, the
   scenario-aware pipeline was never pushed to this checkout — stop and say so.
3. Check the renv library. `USE_RENV` defaults to 0, so the array tasks use
   whatever library `.Rprofile` activates; if that library is missing, every task
   fails.

   Load the same R the jobs will use BEFORE any R call — a bare `Rscript` on a
   login node is a different R than the array gets, searching a different library:

       source chpc/config.sh
       module load "$R_MODULE"

   Then ask R where its libraries are rather than assuming a path. On this
   machine renv (1.1+) keeps the project library OUTSIDE the repo, under
   `~/.cache/R/renv/library/<project>-<hash>/linux-rocky-8.10/R-4.5/...`, so an
   absent `renv/library/` directory inside the repo is NORMAL and means nothing.
   Never conclude anything from `ls renv/library`. Ask R:

       Rscript -e 'cat(.libPaths(), sep="\n")' \
               -e 'cat("fastcmprsk:", system.file(package="fastcmprsk"), "\n")' \
               -e 'cat("tidyverse:", system.file(package="tidyverse"), "\n")'

   Interpret the result:
   - Both packages print a real path -> good, proceed.
   - A package prints empty -> stop, and report it together with the full
     `.libPaths()` output so I can see which library R actually searched.
   - The command errors, or you cannot tell -> WARN loudly, say exactly what you
     saw, and PROCEED with the submission anyway. An inconclusive check is not a
     reason to block; a genuinely broken library shows up in the first task's log
     within seconds, which is cheaper than a false stop.

   Never try to fix the library yourself — `renv::restore()` is a compile job and
   belongs in an `salloc` session, not here.

4. Validate the scenario slug against `[A-Za-z0-9][A-Za-z0-9.-]*` — no
   underscores. If it doesn't match, stop and ask; never silently rewrite it.
5. Pick the run range. If the request gives an explicit start, use it. Otherwise
   list `$SCRATCH_BASE/<scenario>/p<P>/` for each requested `npredictors`, take
   the highest existing `run<N>`, and start at N+1 so this batch appends instead
   of colliding. Empty or missing directory means start at 1. State the range you
   chose in your report.
   An existing scenario directory is normal and expected — scenarios accumulate
   runs across submissions. It is not a reason to stop.
6. Sanity-check the walltime against the workload. `SB_TIME` is per ARRAY TASK,
   and each task runs `RUNS_PER_TASK` runs sequentially. Read both from
   `chpc/config.sh` (or the overrides) and, if `SB_TIME` divided by
   `RUNS_PER_TASK` leaves under ~5 minutes per run, say so plainly in your report
   before submitting. Do not change either value on your own — just flag it.

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
- The resolved repo root, and the run range used
- Array job ID + task count, and the parse job ID
- Scratch path, and the `Sim/chpc_results/<scenario>/` files the parse job will
  write
- Any walltime concern from preflight step 6
- The `squeue` command to track it, with `-M $CHPC_CLUSTER` if that is set

Run `squeue` ONCE to confirm the jobs are queued. Do not poll, wait, or watch.

## Arbiter rule

You are on a login node governed by Arbiter: roughly 15 core-minutes at or under
4 GB before penalties escalate. Almost everything you do here is `git`,
`ls`/`cat`, and `sbatch`.

The line is compute time, not the R binary. A sub-second `Rscript -e` that only
prints a path or a version is fine. What is forbidden here, always, is anything
that COMPILES or FITS: `renv::restore()`, `Rcpp::sourceCpp`,
`chpc/R/run_batch.R`, `chpc/R/parse_batch.R`, or any model fit. Those belong in
the array, or in an `salloc` session for a re-parse.
