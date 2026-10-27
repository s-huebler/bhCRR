# sim/local/single.R -- interactive driver: one run, result returned to caller
#
# OWNS
#   run_single(scenario_file, p, run_id, load_mode, seed,
#              fit_control, s0_seq, s1_seq, max_widen,
#              show_by_block, quiet)
#                     call run_one() for a single (scenario, p, run_id) triple,
#                     capture warnings, print a compact summary, and return the
#                     result invisibly
#
# WHY A FUNCTION (not a source-and-assign script)
#   The original scaffold used set-globals-then-source (SCENARIO_FILE, P, RUN_ID).
#   A function is testable and cannot silently inherit a stale variable from a
#   session that has already been used to debug several runs. batch.R uses the
#   same convention.
#
# USAGE
#   source("sim/local/single.R")
#   result <- run_single(scenario_file = "sim/scenarios/myScenario.R", p = 25, run_id = 1)
#
#   Quick look with cheap fit_control (avoids the production 10-fold x 2-CV cost):
#   ctrl <- bhcrr_cv_control(nfolds = 3, ncv = 1,
#             fit_args = list(maxit = 20L, epsilon = 1e-3))
#   result <- run_single("sim/scenarios/myScenario.R", p = 25, run_id = 1,
#                        fit_control = ctrl, s0_seq = c(0.04, 0.08),
#                        s1_seq = c(0.4, 0.8), max_widen = 0L)
#
# INPUTS
#   scenario_file  character(1) or NULL
#                                path to a scenarios/*.R file; NULL is a
#                                HARD ERROR -- there is no default scenario.
#                                Silently running a toy fixture would produce
#                                numbers that look real.
#   p              integer(1)   number of predictors; default 25
#   run_id         integer(1)   run index; default 1
#   load_mode      character(1) "dev" (default) or "cluster"
#   seed           logical(1) or integer(1)
#                                TRUE derives a deterministic seed from
#                                (scenario, p, run_id); FALSE uses no seed;
#                                any integer uses that seed directly
#   fit_control    bhcrr_cv_control or NULL
#                                NULL uses fit_bhcrr()'s production defaults
#                                (nfolds = 10, ncv = 2). For a quick look,
#                                pass a cheap control (see USAGE above).
#   s0_seq         numeric or NULL   override the s0 tuning grid
#   s1_seq         numeric or NULL   override the s1 tuning grid
#   max_widen      integer(1)        max grid-widening attempts; default 1L
#   show_by_block  logical(1)   TRUE (default): print the per-block score table
#   quiet          logical(1)   TRUE: suppress all printed output
#
# RETURNS (invisibly)
#   The full run_one() result list with one additional field:
#     $warnings  character vector of warning messages collected during the run
#                (empty character(0) if none)
#   The complete result is the return value; `result <- run_single(...)` gives
#   access to every field the printed summary omits.
#
# WARNING CAPTURE
#   run_one() does not capture or suppress warnings. batch_sweep() also does
#   not capture them (deliberately -- it processes many runs and adding a
#   per-run capture would obscure the sweep log). This driver wraps run_one()
#   in withCallingHandlers() so warnings are collected and reported in the
#   summary rather than flooding the console or being truncated by SLURM at 50.
#   The first CHPC array run printed "There were 50 or more warnings" in both
#   tasks and the content was invisible. If the warning summary here proves
#   useful we may move it into run_one() later.
#
# WRITES NOTHING
#   run_single() does not write any files and does not create any directories.
#   It never touches sim/results/.
#
# DEPENDS ON
#   sim/R/load.r    repo_root(), load_sim()
#   sim/R/run_one.r run_one()  (transitively all of sim/R/)
#
# REPLACES
#   sim/_legacy/code/single_run_complex.R
#   sim/_legacy/code/single_run.R

# ---------------------------------------------------------------------------
# run_single
# ---------------------------------------------------------------------------

run_single <- function(scenario_file = NULL,
                       p             = 25L,
                       run_id        = 1L,
                       load_mode     = Sys.getenv("BHCRR_LOAD", "dev"),
                       seed          = TRUE,
                       fit_control   = NULL,
                       s0_seq        = NULL,
                       s1_seq        = NULL,
                       max_widen     = 1L,
                       show_by_block = TRUE,
                       quiet         = FALSE) {

  # ---- B2: Bootstrap loader ----
  # The only place in sim/local/ where a source() into sim/R/ is used; something
  # must bootstrap the loader before any sim/R/ function is available.
  .wd <- normalizePath(getwd(), mustWork = TRUE)
  .p  <- .wd
  repeat {
    candidate <- file.path(.p, "sim", "R", "load.r")
    if (file.exists(candidate)) { source(candidate); break }
    .pp <- dirname(.p)
    if (identical(.pp, .p))
      stop("Cannot locate sim/R/load.r. Run from inside the bhCRR repo.")
    .p <- .pp
  }
  load_sim(load_mode, force = FALSE)

  # ---- B1: scenario_file required ----
  if (is.null(scenario_file)) {
    avail <- list.files(
      file.path(repo_root(), "sim", "scenarios"),
      pattern = "^[^_].*\\.R$", full.names = FALSE
    )
    stop(
      "scenario_file is required -- there is no default scenario.\n",
      "  Pass a path to a sim/scenarios/<slug>.R file.\n",
      "  To create a new scenario: new_scenario(\"<name>\", dir = \"sim/scenarios\")\n",
      if (length(avail))
        paste0("  Existing scenarios: ",
               paste(avail, collapse = ", "))
      else
        "  (no scenario files found in sim/scenarios/ yet)",
      call. = FALSE
    )
  }

  if (!file.exists(scenario_file))
    stop("Scenario file not found: ", scenario_file, call. = FALSE)

  # ---- B4: Warning capture ----
  warns <- character(0)
  t0 <- proc.time()[["elapsed"]]

  result <- withCallingHandlers(
    run_one(
      scenario    = scenario_file,
      p           = as.integer(p),
      run_id      = as.integer(run_id),
      load_mode   = load_mode,
      seed        = seed,
      fit_control = fit_control,
      s0_seq      = s0_seq,
      s1_seq      = s1_seq,
      max_widen   = as.integer(max_widen)
    ),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  result$warnings <- warns

  # ---- B5: Summary ----
  if (!quiet) {
    .sep <- function(w = 60) cat(strrep("-", w), "\n")

    .sep()
    cat(sprintf("run_single()  scenario: %s   p: %d   run_id: %d\n",
                result$scenario$name, result$p, result$run_id))
    cat(sprintf("  seed used: %s\n",
                if (is.na(result$meta$seed)) "none (seed = FALSE)"
                else as.character(result$meta$seed)))

    # -- bhCRR tuning --
    .sep()
    cat("bhCRR tuning:\n")
    cat(sprintf("  s0: %g   s1: %g   n_attempts: %s   all_pass: %s\n",
                result$fit$s0, result$fit$s1,
                result$fit$n_attempts,
                result$fit$all_pass))
    if (!is.null(result$fit$flags) && length(result$fit$flags)) {
      cat("  flags:\n")
      print(result$fit$flags)
    } else {
      cat("  flags: (none)\n")
    }

    # -- eval horizon --
    eval_time_val <- tryCatch(
      unique(result$scores$overall$eval_time),
      error = function(e) NA
    )
    if (length(eval_time_val) == 1L && !is.na(eval_time_val))
      cat(sprintf("  eval_time: %g\n", eval_time_val))

    # -- overall scores --
    .sep()
    cat("Scores (overall, all models):\n")
    ov <- result$scores$overall
    # Drop eval_time from display if constant (already printed above)
    if (!is.null(ov$eval_time) && length(unique(ov$eval_time)) == 1L)
      ov$eval_time <- NULL
    print(ov, row.names = FALSE, digits = 4)

    # -- by_block --
    if (show_by_block && !is.null(result$scores$by_block)) {
      .sep()
      cat("Scores (by block):\n")
      print(result$scores$by_block, row.names = FALSE, digits = 4)
    }

    # -- timing --
    .sep()
    setup_s <- tryCatch(as.numeric(result$comparators$setup_seconds),
                        error = function(e) NA_real_)
    cat(sprintf("Timing:\n"))
    cat(sprintf("  elapsed (run_one): %.1f s\n", result$meta$elapsed_seconds))
    if (!is.na(setup_s))
      cat(sprintf("  comparator setup: %.1f s  (shared across lasso/alasso/scad/mcp; count once)\n",
                  setup_s))

    # -- warning summary --
    .sep()
    if (length(warns) == 0L) {
      cat("Warnings: none\n")
    } else {
      tbl <- sort(table(warns), decreasing = TRUE)
      cat(sprintf("Warnings: %d total, %d distinct:\n",
                  length(warns), length(tbl)))
      for (i in seq_along(tbl))
        cat(sprintf("  x%d  %s\n", tbl[[i]], names(tbl)[i]))
    }

    .sep()
    cat("Complete result is the return value. Use `result <- run_single(...)` to access all fields.\n")
  }

  invisible(result)
}
