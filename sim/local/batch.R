# sim/local/batch.R -- local sweep driver: (scenario, npredictors, run range)
#
# OWNS
#   batch_sweep(scenario_file, npredictors, run_start, run_end,
#               results_root, load_mode, seed, fit_control,
#               s0_seq, s1_seq, max_widen, overwrite, dryrun)
#                     sweep over every combination of p in npredictors and
#                     run_id in run_start:run_end; call run_one() for each;
#                     save the result to the canonical path; handle per-run
#                     failures without stopping the sweep
#
# USAGE
#   source("sim/local/batch.R")
#   batch_sweep(scenario_file = "sim/scenarios/myScenario.R",
#               npredictors   = c(25, 50),
#               run_start     = 1,
#               run_end       = 10)
#
# This file does NOT auto-run on source().
#
# ON-DISK LAYOUT
#   sim/results/<scenario>/runs/p<P>/p<P>_<scenario>_run<N>.Rdata
#   sim/results/<scenario>/parsed/
#
#   The p<P>/ level is load-bearing: the column-to-block map is a function
#   of p, so runs at different p must never be pooled. runs/ is gitignored
#   (sim/results/*/runs/**) so raw fits stay out of the repo; parsed/ stays
#   committed so summary tables are versioned.
#
#   The filename format p<P>_<scenario>_run<N>.Rdata is exactly what
#   sim/chpc/R/parse_batch.R matches with:
#     ^p(\d+)_([A-Za-z0-9][A-Za-z0-9.-]*)_run(\d+)\.Rdata$
#   Local and CHPC results are therefore parseable by the same script.
#
# INPUTS
#   scenario_file  character(1)   path to a scenarios/*.R file
#   npredictors    integer(1+)    vector of p values to sweep
#   run_start      integer(1)     first run index, inclusive; default 1
#   run_end        integer(1)     last run index, inclusive; default 10
#   results_root   character(1) or NULL
#                                 root of the results tree;
#                                 NULL -> file.path(repo_root(), "sim", "results")
#   load_mode      character(1)   passed to load_bhcrr(); "dev" or "cluster";
#                                 default: Sys.getenv("BHCRR_LOAD", "dev")
#   seed           logical(1) or integer(1)
#                                 passed to run_one(); TRUE derives a
#                                 deterministic seed per (scenario, p, run_id)
#   fit_control    bhcrr_cv_control or NULL   passed to run_one()
#   s0_seq         numeric or NULL            passed to run_one()
#   s1_seq         numeric or NULL            passed to run_one()
#   max_widen      integer(1)                 passed to run_one(); default 1L
#   overwrite      logical(1)     FALSE (default): refuse the sweep if any
#                                 output file already exists; see collision guard
#   dryrun         logical(1)     TRUE: validate, check manifest, print what
#                                 would be written, but fit nothing and write
#                                 nothing; FALSE (default): run the sweep
#
# OUTPUTS
#   One file per run:
#     <results_root>/<scenario>/runs/p<P>/p<P>_<scenario>_run<N>.Rdata
#   Each contains a single object `result` (output of run_one()), saved with
#   compress = "xz" to match sim/chpc/R/run_batch.R.
#   Failed runs produce a failure record at the same path (see FAILURE RECORD).
#
#   Manifest on first run:
#     <results_root>/<scenario>/scenario.env
#   On subsequent runs with the same scenario, check_manifest() validates it.
#
# FAILURE RECORD
#   If run_one() errors, batch_sweep() saves the following at the canonical
#   path so the failure is visible on disk rather than appearing as a silent gap:
#     list(failed        = TRUE,
#          error         = conditionMessage(e),
#          meta          = list(run_id, p, scenario_name, git_sha, seed,
#                               timestamp))
#   The seed recorded is the seed run_one() would have used; the failure is
#   therefore reproducible via run_one(seed = <seed>).
#   parse.r (step 11) must detect $failed == TRUE and skip such records while
#   reporting the count; it must NOT treat them as corrupt files.
#
# PARALLELISM
#   This driver runs one call at a time. Parallelising the local sweep is a
#   separate question and is not addressed here. run_one() sets a seed derived
#   from (scenario, p, run_id), so a future parallel sweep would still be
#   reproducible per run without coordination.
#
# ERROR POLICY
#   batch_sweep() is the only sim/ code that catches errors. It wraps each
#   run_one() call in tryCatch, saves a failure record, and continues the
#   sweep. run_one() itself does not catch; it propagates so this driver owns
#   all failure handling.
#
# DEPENDS ON
#   sim/R/load.r      repo_root(), load_bhcrr()
#   sim/R/scenario.r  scenario_from_file(), write_manifest(), check_manifest()
#   sim/R/run_one.r   run_one()  (transitively all of sim/R/)
#
# REPLACES
#   sim/_legacy/code/batch_run_complex.R   (local batch use case)
#   sim/_legacy/code/batch_run.R           (idem, simpler variant)
#   sim/_legacy/code/batch_run_wolbers.R   (idem, Wolbers scoring variant)

# ---------------------------------------------------------------------------
# batch_sweep
# ---------------------------------------------------------------------------

batch_sweep <- function(scenario_file,
                        npredictors,
                        run_start    = 1L,
                        run_end      = 10L,
                        results_root = NULL,
                        load_mode    = Sys.getenv("BHCRR_LOAD", "dev"),
                        seed         = TRUE,
                        fit_control  = NULL,
                        s0_seq       = NULL,
                        s1_seq       = NULL,
                        max_widen    = 1L,
                        overwrite    = FALSE,
                        dryrun       = FALSE) {

  t_sweep_start <- proc.time()[["elapsed"]]

  # ---- Load bhCRR and sim modules (idempotent) ----
  # Call load_bhcrr, not load_sim -- load_sim sources run_one.r, and run_one.r
  # calls load_bhcrr internally; calling load_sim here would re-source modules
  # that are already in scope, which is harmless but wasteful. Drivers that
  # call batch_sweep after a fresh source("sim/R/load.r") + load_sim() call
  # will already have everything attached; load_bhcrr is idempotent.
  source_dir <- tryCatch(
    dirname(normalizePath(sys.frame(0)$ofile, mustWork = FALSE)),
    error = function(e) NULL
  )
  if (!is.null(source_dir) && file.exists(file.path(source_dir, "..", "R", "load.r"))) {
    source(file.path(source_dir, "..", "R", "load.r"))
    load_sim(load_mode, force = FALSE)
  } else {
    # Fallback: walk up from working directory
    .wd <- normalizePath(getwd(), mustWork = TRUE)
    .p <- .wd
    repeat {
      candidate <- file.path(.p, "sim", "R", "load.r")
      if (file.exists(candidate)) { source(candidate); load_sim(load_mode, force = FALSE); break }
      candidate2 <- file.path(.p, "R", "load.r")
      if (file.exists(candidate2)) { source(candidate2); load_sim(load_mode, force = FALSE); break }
      .pp <- dirname(.p)
      if (identical(.pp, .p)) stop("Cannot locate sim/R/load.r; set working directory inside the bhCRR repo.")
      .p <- .pp
    }
  }

  # ---- C1: Resolve and validate scenario ----
  spec <- scenario_from_file(scenario_file)
  scenario_name <- spec$name

  # ---- Resolve results root ----
  if (is.null(results_root))
    results_root <- file.path(repo_root(), "sim", "results")
  scenario_dir <- file.path(results_root, scenario_name)
  runs_dir     <- file.path(scenario_dir, "runs")
  manifest_dir <- scenario_dir   # scenario.env lives directly in scenario_dir

  # ---- C2: Manifest: write on first run, check on subsequent ----
  manifest_path <- file.path(manifest_dir, "scenario.env")
  if (!dryrun) {
    if (file.exists(manifest_path)) {
      check_manifest(spec, manifest_dir)   # hard-stops on any signature mismatch
    } else {
      dir.create(manifest_dir, showWarnings = FALSE, recursive = TRUE)
      write_manifest(spec, manifest_dir)
      message(sprintf("[batch] manifest written to %s", manifest_path))
    }
  }

  # ---- C3: Collision guard ----
  # Collect all planned output paths and check for pre-existing files.
  run_ids  <- seq.int(as.integer(run_start), as.integer(run_end))
  all_tags <- expand.grid(p = npredictors, run_id = run_ids,
                          stringsAsFactors = FALSE)

  planned <- mapply(function(p, run_id) {
    tag <- sprintf("p%d_%s_run%d", p, scenario_name, run_id)
    file.path(runs_dir, sprintf("p%d", p), paste0(tag, ".Rdata"))
  }, all_tags$p, all_tags$run_id, SIMPLIFY = FALSE)
  names(planned) <- sprintf("p%d_run%d", all_tags$p, all_tags$run_id)

  if (!overwrite) {
    existing <- Filter(file.exists, planned)
    if (length(existing) > 0L) {
      # Group by p for a readable message
      existing_tags <- names(existing)
      for (pval in sort(unique(all_tags$p))) {
        collide <- existing_tags[grepl(sprintf("^p%d_run", pval), existing_tags)]
        if (length(collide) > 0L) {
          run_nums <- as.integer(sub("^p\\d+_run", "", collide))
          message(sprintf(
            "[batch] COLLISION: p=%d already has run(s): %s",
            pval, paste(sort(run_nums), collapse = ", ")))
        }
      }
      stop(sprintf(
        "Output files already exist for %d run(s). ",
        length(existing)),
        "Use overwrite = TRUE to overwrite, ",
        "or run_start = <next free index> to append.")
    }
  }

  # ---- dryrun: print plan and exit ----
  if (dryrun) {
    message(sprintf("[batch] DRY RUN -- scenario: %s  p: %s  runs: %d..%d",
                    scenario_name, paste(npredictors, collapse = ","),
                    run_start, run_end))
    message(sprintf("[batch] results_root: %s", results_root))
    message("[batch] Files that would be written:")
    for (path in planned)
      message(sprintf("  %s", path))
    message("[batch] No fitting performed.")
    return(invisible(list(dryrun = TRUE, planned = unname(planned))))
  }

  # ---- C4 / C5: Sweep ----
  n_attempted <- 0L
  n_succeeded <- 0L
  n_failed    <- 0L

  # Derive git SHA once for failure records
  git_sha <- tryCatch({
    trimws(system2("git", c("-C", repo_root(), "rev-parse", "--short", "HEAD"),
                   stdout = TRUE, stderr = FALSE))
  }, error = function(e) "unknown", warning = function(e) "unknown")
  if (!nzchar(git_sha)) git_sha <- "unknown"

  for (p in npredictors) {
    out_dir_p <- file.path(runs_dir, sprintf("p%d", p))
    dir.create(out_dir_p, showWarnings = FALSE, recursive = TRUE)

    for (run_id in run_ids) {
      n_attempted <- n_attempted + 1L

      tag      <- sprintf("p%d_%s_run%d", p, scenario_name, run_id)
      out_file <- file.path(out_dir_p, paste0(tag, ".Rdata"))

      message(sprintf("[%s] %s starting...",
                      format(Sys.time(), "%H:%M:%S"), tag))

      # Derive the seed we will pass so the failure record matches
      run_seed <- if (isTRUE(seed)) {
        # Mirror run_one.r's .derive_seed logic exactly
        s_str <- paste0(scenario_name, "|p=", as.integer(p),
                        "|run=", as.integer(run_id))
        bytes <- utf8ToInt(s_str)
        m <- 2147483647.0
        h <- 5381.0
        for (b in bytes) h <- (h * 33.0 + b) %% m
        as.integer(h)
      } else if (isFALSE(seed)) {
        NA_integer_
      } else {
        as.integer(seed)
      }

      result <- tryCatch(
        run_one(spec, p = p, run_id = run_id,
                load_mode   = load_mode,
                seed        = seed,
                fit_control = fit_control,
                s0_seq      = s0_seq,
                s1_seq      = s1_seq,
                max_widen   = max_widen),
        error = function(e) {
          message(sprintf("  !! %s FAILED: %s", tag, conditionMessage(e)))
          list(failed = TRUE,
               error  = conditionMessage(e),
               meta   = list(run_id       = as.integer(run_id),
                             p            = as.integer(p),
                             scenario_name = scenario_name,
                             git_sha      = git_sha,
                             seed         = run_seed,
                             timestamp    = Sys.time()))
        }
      )

      save(result, file = out_file, compress = "xz")
      message(sprintf("  -> saved %s", out_file))

      if (isTRUE(result$failed)) {
        n_failed    <- n_failed + 1L
      } else {
        n_succeeded <- n_succeeded + 1L
      }
      rm(result)   # drop immediately; do not accumulate in memory
    }
  }

  elapsed <- as.numeric(proc.time()[["elapsed"]] - t_sweep_start)

  # ---- C5: Summary ----
  message(sprintf(
    "\n[batch] Done. attempted=%d  succeeded=%d  failed=%d  elapsed=%.1fs",
    n_attempted, n_succeeded, n_failed, elapsed))
  message(sprintf("[batch] Output directory: %s", runs_dir))

  invisible(list(
    scenario_name = scenario_name,
    npredictors   = npredictors,
    run_start     = as.integer(run_start),
    run_end       = as.integer(run_end),
    n_attempted   = n_attempted,
    n_succeeded   = n_succeeded,
    n_failed      = n_failed,
    elapsed       = elapsed,
    runs_dir      = runs_dir
  ))
}
