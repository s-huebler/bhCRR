# sim/R/parse.r -- many result files -> CSV and .Rdata summary tables
#
# OWNS
#   parse_results(runs_root, scenario_name, out_dir = NULL, manifest_dir = NULL)
#                       discover .Rdata run files under runs_root/p<P>/,
#                       validate them, assemble long tables from the
#                       precomputed metrics in $scores and the saved
#                       coefficient vectors, write CSVs and .Rdata files
#
# DESIGN: PURE STACKER
#   parse_results() COMPUTES NO METRICS. Every metric (wolbers_c, tp/fp/tn/fn,
#   precision, recall, etc.) arrives precomputed in each run's $scores$overall
#   and $scores$by_block. parse.r discovers, validates, assembles, and writes.
#   Assembling is allowed; recomputing is not. Adding a new metric requires only
#   a change to score.r -- this file never needs to change for new metrics.
#
# INPUTS
#   runs_root     character(1)  directory CONTAINING the p<P>/ subdirs:
#                               local: sim/results/<scenario>/runs
#                               CHPC:  $SCRATCH_RUN_DIR
#   scenario_name character(1)  scenario slug, e.g. "scenarioA"
#   out_dir       character(1) or NULL
#                               where to write the summary files;
#                               NULL -> file.path(repo_root(), "sim", "results",
#                                        scenario_name, "parsed")
#   manifest_dir  character(1) or NULL
#                               directory containing scenario.env for the
#                               cross-check; NULL auto-searches runs_root
#                               then dirname(runs_root) to handle both the
#                               local layout (<scenario>/runs/) and the CHPC
#                               layout ($SCRATCH_RUN_DIR flat)
#
# CALL SITES
#   Local:  source("sim/R/load.r"); load_sim()
#           parse_results("sim/results/scenarioA/runs", "scenarioA")
#   CHPC:   sim/chpc/R/parse.r (thin wrapper; step 12)
#
# RETURNS
#   invisible(list) with:
#     $written        character  paths to every file written
#     $n_runs         integer    number of successfully parsed runs
#     $n_failed       integer    number of failure records skipped
#     $n_skipped      integer    number of non-matching files skipped
#     $failed_summary data.frame p, run_id, error for each failure record
#
# OUTPUT FILES (all written to out_dir)
#
#   betas_p<P>_<scenario>.csv / .Rdata
#     Long, one row per predictor × model × run:
#       scenario, predictor, Block, BlockRho, True, Model, Estimate, Bias,
#       nobs, run_number
#     The .Rdata holds an object named betas_p<P> (legacy convention).
#     Block / BlockRho come from $meta$block_map, joined BY PREDICTOR NAME.
#     Estimate is always the raw-scale coefficient (coef_raw / $fit$coef),
#     never coef_std. Bias = Estimate - True.
#
#   model_performance_p<P>_<scenario>.csv
#     One row per model per run: all $scores$overall columns (model renamed
#     to Model for legacy compatibility), with scenario, run_number, nobs,
#     npredictors prepended, plus: s0, s1, n_attempts, all_pass (bhCRR-only),
#     converged (per model), iterations (bhCRR-only), setup_seconds
#     (comparators-only). eval_time travels as a column of $scores$overall
#     already, so every row records the horizon it was scored at.
#
#   by_block_p<P>_<scenario>.csv
#     Long per-block frame, stacked $scores$by_block across runs, with
#     scenario and run_number prepended. model/block/block_rho renamed to
#     Model/Block/BlockRho. NOT widened: the long shape is more analysis-
#     friendly than legacy's wide TP_strong/FP_strong/... columns.
#
#   false_positives_p<P>_<scenario>.csv
#     For each (Model, Block, BlockRho), sum fp and (fp + tn) across runs:
#       scenario, Model, Block, BlockRho, n_false_pos, denominator, fp_rate
#     fp_rate = sum(fp) / sum(fp + tn). Pooling counts then dividing is
#     correct because the block map is deterministic and the beta pattern is
#     fixed: every run of a (scenario, p) group has the same number of true
#     zeros per block. fp/(fp+tn) = 1 - specificity restricted to the block's
#     true zeros, which is the key diagnostic for selection leakage across
#     correlated neighbours.
#
#   model_performance_<scenario>.csv
#     Only when more than one p-group is present: model_performance stacked
#     across all p values. Per-block and betas files are never stacked across
#     p because the block map is a function of p -- a predictor name refers
#     to a different covariate at p=25 vs p=50.
#
#   scenario.env (copied from runs_root/../scenario.env if present)
#     Makes out_dir self-describing without the runs tree.
#
# FG_ORACLE NOTE
#   fg_oracle is a MODEL, not a column. Legacy's betas table carried an FG
#   column holding the oracle estimate with NA outside the active set. score.r
#   makes fg_oracle one of six model rows. This file drops the FG column and
#   lets fg_oracle appear as its own rows with Model == "fg_oracle", exactly
#   like every other method. Its per-predictor estimates are reachable by
#   filtering Model == "fg_oracle".
#
# HARD STOP 1 -- scenario isolation (B2)
#   If any filename's scenario token or any loaded run's $meta$scenario_name
#   differs from scenario_name, stop naming both values and the offending file.
#
# HARD STOP 2 -- block signature consistency (B3)
#   Within one p-group, every run's $meta$block_signature must be identical.
#   On mismatch, print both signatures in full (they are human-readable
#   exactly for this purpose) and the files they came from, then stop.
#   Signatures across p-groups are expected to differ and are not compared.
#
# HARD STOP 3 -- nobs consistency (B4)
#   $meta$design_settings$nobs (cross-checked against $scenario$nobs) must
#   be identical across every run in the scenario. Stop naming both values
#   and the files.
#
# FAILURE RECORDS (B6)
#   batch_sweep() writes failure records with $failed == TRUE at the canonical
#   path. These are NOT corrupt files. Detected by $failed, skipped, counted,
#   and reported in the summary (error message and count). They never appear
#   in any output table and never trigger a hard stop.
#
# GLOBAL ENVIRONMENT NOTE
#   parse_results() loads run files into isolated environments (new.env()) and
#   never assigns anything to the global environment. Legacy's pattern of
#   assigning each loaded result to a named variable in globalenv() is
#   deliberately removed here: it pollutes the workspace and is incompatible
#   with being called from a SLURM script.
#
# DEPENDS ON
#   sim/R/load.r      repo_root()
#   sim/R/scenario.r  check_manifest()
#   bhCRR:: (none; parsing is pure R)
#
# REPLACES
#   sim/_legacy/code/parse_batch_run.R   (entire file)
#   sim/chpc/R/parse_batch.R             (entire file; becomes a thin wrapper
#                                         at migration step 12)

# ---------------------------------------------------------------------------
# .parse_rx -- canonical filename regex (identical to parse_batch.R:147)
# ---------------------------------------------------------------------------

.parse_rx <- "^p(\\d+)_([A-Za-z0-9][A-Za-z0-9.-]*)_run(\\d+)\\.Rdata$"

# ---------------------------------------------------------------------------
# parse_results
# ---------------------------------------------------------------------------

parse_results <- function(runs_root, scenario_name, out_dir = NULL, manifest_dir = NULL) {

  # ---- Resolve out_dir ----
  if (is.null(out_dir))
    out_dir <- file.path(repo_root(), "sim", "results", scenario_name, "parsed")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  message(sprintf("parse_results  scenario=%s", scenario_name))
  message(sprintf("  runs_root: %s", runs_root))
  message(sprintf("  out_dir:   %s", out_dir))

  # ---- B1: Discover run files ----
  all_files <- list.files(runs_root, pattern = "\\.Rdata$",
                          full.names = TRUE, recursive = TRUE)
  if (length(all_files) == 0L)
    stop("No .Rdata files found under ", runs_root,
         " (searched recursively, including p<P>/ subdirs).")

  rx_matches <- regmatches(basename(all_files),
                           regexec(.parse_rx, basename(all_files)))
  matched <- vapply(rx_matches, length, integer(1)) > 0L

  n_skipped_nm <- sum(!matched)
  if (n_skipped_nm > 0L)
    message(sprintf(
      "Skipping %d file(s) that do not match p<P>_<scenario>_run<N>.Rdata:\n  %s",
      n_skipped_nm,
      paste(basename(all_files)[!matched], collapse = "\n  ")))
  if (!any(matched))
    stop("No files matching p<P>_<scenario>_run<N>.Rdata found under ", runs_root, ".")

  files      <- all_files[matched]
  rx_matches <- rx_matches[matched]

  file_p    <- as.integer(vapply(rx_matches, function(m) m[[2L]], character(1)))
  file_scen <- vapply(rx_matches, function(m) m[[3L]], character(1))
  file_run  <- as.integer(vapply(rx_matches, function(m) m[[4L]], character(1)))

  # B2: HARD STOP -- scenario token mismatch
  bad_sc <- file_scen != scenario_name
  if (any(bad_sc))
    stop(sprintf(
      "Scenario token mismatch in %d file(s) -- expected '%s':\n%s",
      sum(bad_sc), scenario_name,
      paste(sprintf("  %s  (token: '%s')",
                    basename(files)[bad_sc], file_scen[bad_sc]),
            collapse = "\n")))

  # Sort by p then run
  ord      <- order(file_p, file_run)
  files    <- files[ord]
  file_p   <- file_p[ord]
  file_run <- file_run[ord]

  # ---- B5: manifest check setup ----
  # Resolve manifest_dir. NULL triggers an automatic search because the two
  # layouts differ:
  #   Local: runs_root = <scenario>/runs/  -> manifest at <scenario>/
  #          (manifest_dir = dirname(runs_root))
  #   CHPC:  runs_root = $SCRATCH_RUN_DIR  -> manifest alongside runs
  #          (manifest_dir = runs_root itself)
  # Search order: runs_root first, then dirname(runs_root); use the first
  # that contains scenario.env. The CHPC wrapper passes manifest_dir
  # explicitly rather than relying on this search.
  if (is.null(manifest_dir)) {
    if (file.exists(file.path(runs_root, "scenario.env"))) {
      manifest_dir <- runs_root
    } else {
      manifest_dir <- dirname(runs_root)
    }
  }
  manifest_path <- file.path(manifest_dir, "scenario.env")
  has_manifest  <- file.exists(manifest_path)
  if (!has_manifest)
    warning("No scenario.env found in ", manifest_dir,
            "; skipping manifest cross-check.")
  manifest_checked <- FALSE

  # ---- Main loop ----
  # Accumulators
  beta_tables  <- list()   # "p<P>" -> data.frame
  perf_tables  <- list()   # "p<P>" -> data.frame
  block_tables <- list()   # "p<P>" -> data.frame
  block_sigs   <- list()   # "p<P>" -> list(sig, file)
  nobs_check   <- NULL
  nobs_file    <- NULL
  n_failed     <- 0L
  failed_info  <- list()
  n_runs       <- 0L
  model_names  <- c("bhcrr","fg_oracle","lasso","alasso","scad","mcp")

  for (i in seq_along(files)) {
    f     <- files[i]
    p_f   <- file_p[i]
    run_f <- file_run[i]
    key   <- sprintf("p%d", p_f)
    nm    <- tools::file_path_sans_ext(basename(f))

    # Load into isolated environment -- never globalenv (see GLOBAL ENVIRONMENT NOTE)
    e <- new.env(parent = emptyenv())
    load(f, envir = e)
    res <- e$result

    # B6: failure record
    if (isTRUE(res$failed)) {
      n_failed <- n_failed + 1L
      err_msg  <- as.character(
        if (!is.null(res$error)) res$error else "unknown error")
      failed_info[[n_failed]] <- list(
        p      = p_f,
        run_id = run_f,
        error  = err_msg)
      message(sprintf("  [skip] %s -- failure record: %s", nm, err_msg))
      next
    }

    # B2: meta scenario_name check (belt-and-suspenders; filename already checked)
    if (!is.null(res$meta$scenario_name) &&
        res$meta$scenario_name != scenario_name)
      stop(sprintf(
        "Scenario mismatch inside '%s': $meta$scenario_name == '%s', expected '%s'.",
        nm, res$meta$scenario_name, scenario_name))

    # B4: HARD STOP -- nobs consistency
    nobs_f <- tryCatch(
      as.integer(res$meta$design_settings$nobs),
      error = function(e) NA_integer_
    )
    if (is.null(nobs_f) || is.na(nobs_f))
      nobs_f <- tryCatch(as.integer(res$scenario$nobs), error = function(e) NA_integer_)
    if (is.null(nobs_check)) {
      nobs_check <- nobs_f
      nobs_file  <- nm
    } else if (!identical(nobs_check, nobs_f)) {
      stop(sprintf(
        "nobs mismatch: '%s' has nobs=%d; earlier runs had nobs=%d (first in '%s').",
        nm, nobs_f, nobs_check, nobs_file))
    }

    # B3: HARD STOP -- block signature consistency within p-group
    sig <- res$meta$block_signature
    if (!is.null(sig)) {
      if (is.null(block_sigs[[key]])) {
        block_sigs[[key]] <- list(sig = sig, file = nm)
      } else if (!identical(block_sigs[[key]]$sig, sig)) {
        stop(sprintf(paste0(
          "block_signature mismatch within group %s:\n",
          "  %s:\n    %s\n",
          "  %s:\n    %s\n",
          "Per-block summaries would pool incomparable strata. ",
          "These runs cannot be summarised together."),
          key,
          block_sigs[[key]]$file, block_sigs[[key]]$sig,
          nm, sig))
      }
    }

    # B5: manifest cross-check (once, on first valid run)
    if (has_manifest && !manifest_checked && !is.null(res$scenario)) {
      tryCatch(
        check_manifest(res$scenario, manifest_dir),
        error = function(e)
          stop("Manifest mismatch for run '", nm, "': ", conditionMessage(e))
      )
      manifest_checked <- TRUE
    }

    n_runs <- n_runs + 1L

    # ---- C1: Betas table ----
    bm <- res$meta$block_map   # data.frame: name, block, block_rho, ...
    if (is.null(bm))
      stop(sprintf("%s: $meta$block_map is NULL; cannot join Block/BlockRho.", nm))

    true_beta  <- as.numeric(res$data$beta1)
    pred_names <- names(res$data$beta1)
    if (is.null(pred_names))
      pred_names <- paste0("X_", seq_along(true_beta))

    # Join block info by predictor name
    bm_idx <- match(pred_names, bm$name)
    if (anyNA(bm_idx))
      stop(sprintf(
        "%s: predictor(s) absent from block_map: %s",
        nm, paste(pred_names[is.na(bm_idx)], collapse = ", ")))

    block_col     <- bm$block[bm_idx]
    block_rho_col <- bm$block_rho[bm_idx]

    # Coefficient vectors (all raw-scale, all named X_1..X_p)
    model_coefs <- list(
      bhcrr     = res$fit$coef,
      fg_oracle = res$fg_oracle$coef_raw,
      lasso     = res$comparators$lasso$coef_raw,
      alasso    = res$comparators$alasso$coef_raw,
      scad      = res$comparators$scad$coef_raw,
      mcp       = res$comparators$mcp$coef_raw
    )

    betas_run <- do.call(rbind, lapply(model_names, function(mod) {
      est_named <- model_coefs[[mod]]
      # Align to pred_names order (all vectors are named; use names)
      est <- if (!is.null(names(est_named))) {
        as.numeric(est_named[pred_names])
      } else {
        as.numeric(est_named)
      }
      data.frame(
        scenario   = scenario_name,
        predictor  = pred_names,
        Block      = block_col,
        BlockRho   = block_rho_col,
        True       = true_beta,
        Model      = mod,
        Estimate   = est,
        Bias       = est - true_beta,
        nobs       = nobs_f,
        run_number = run_f,
        stringsAsFactors = FALSE
      )
    }))

    beta_tables[[key]] <- rbind(beta_tables[[key]], betas_run)

    # ---- C2: model_performance (stacked $scores$overall) ----
    ov <- res$scores$overall
    # Rename model -> Model for legacy column-name contract
    names(ov)[names(ov) == "model"] <- "Model"

    # Prepend metadata
    ov <- cbind(
      data.frame(scenario    = scenario_name,
                 run_number  = run_f,
                 nobs        = nobs_f,
                 npredictors = p_f,
                 stringsAsFactors = FALSE),
      ov,
      stringsAsFactors = FALSE
    )

    # bhCRR-specific columns (NA for all non-bhCRR rows)
    ov$s0            <- NA_real_
    ov$s1            <- NA_real_
    ov$n_attempts    <- NA_integer_
    ov$all_pass      <- NA
    ov$converged     <- NA
    ov$iterations    <- NA_integer_
    ov$setup_seconds <- NA_real_

    bhcrr_idx <- which(ov$Model == "bhcrr")
    if (length(bhcrr_idx)) {
      ov$s0[bhcrr_idx]         <- as.numeric(res$fit$s0)
      ov$s1[bhcrr_idx]         <- as.numeric(res$fit$s1)
      ov$n_attempts[bhcrr_idx] <- as.integer(res$fit$n_attempts)
      ov$all_pass[bhcrr_idx]   <- isTRUE(res$fit$all_pass)
      ov$converged[bhcrr_idx]  <- isTRUE(res$fit$converged)
      ov$iterations[bhcrr_idx] <- tryCatch(
        as.integer(res$fit$fit$iterations), error = function(e) NA_integer_)
    }

    oracle_idx <- which(ov$Model == "fg_oracle")
    if (length(oracle_idx))
      ov$converged[oracle_idx] <- isTRUE(res$fg_oracle$converged)

    # Comparators: converged per model; setup_seconds shared
    for (comp_nm in c("lasso","alasso","scad","mcp")) {
      idx <- which(ov$Model == comp_nm)
      if (length(idx)) {
        ov$converged[idx]     <- isTRUE(res$comparators[[comp_nm]]$converged)
        ov$setup_seconds[idx] <- as.numeric(res$comparators$setup_seconds)
      }
    }

    rownames(ov) <- NULL
    perf_tables[[key]] <- rbind(perf_tables[[key]], ov)

    # ---- C3: by_block (stacked $scores$by_block) ----
    bb <- res$scores$by_block
    names(bb)[names(bb) == "model"]     <- "Model"
    names(bb)[names(bb) == "block"]     <- "Block"
    names(bb)[names(bb) == "block_rho"] <- "BlockRho"
    bb$scenario   <- scenario_name
    bb$run_number <- run_f
    bb <- bb[, c("scenario","run_number","Model","Block","BlockRho",
                 setdiff(names(bb),
                         c("scenario","run_number","Model","Block","BlockRho"))),
             drop = FALSE]
    rownames(bb) <- NULL
    block_tables[[key]] <- rbind(block_tables[[key]], bb)
  }

  # ---- Write per-group tables ----
  written <- character(0)

  # C4: derive false_positives from stacked by_block
  fp_tables <- lapply(names(block_tables), function(key) {
    bb <- block_tables[[key]]
    # Sum fp and (fp + tn) across runs; pooling counts is valid because
    # the block map is deterministic (same n_true_zeros per block per run)
    agg <- stats::aggregate(
      cbind(sum_fp = bb$fp, sum_denom = bb$fp + bb$tn) ~
        Model + Block + BlockRho,
      data = bb,
      FUN  = sum
    )
    agg$fp_rate      <- agg$sum_fp / agg$sum_denom
    agg$n_false_pos  <- agg$sum_fp
    agg$denominator  <- agg$sum_denom
    agg$sum_fp       <- NULL
    agg$sum_denom    <- NULL
    agg$scenario     <- scenario_name
    agg[, c("scenario","Model","Block","BlockRho",
            "n_false_pos","denominator","fp_rate")]
  })
  names(fp_tables) <- names(block_tables)

  for (key in names(beta_tables)) {
    p_label <- as.integer(sub("^p", "", key))

    # -- betas CSV + Rdata --
    base_b   <- sprintf("betas_%s_%s", key, scenario_name)
    csv_b    <- file.path(out_dir, paste0(base_b, ".csv"))
    rdata_b  <- file.path(out_dir, paste0(base_b, ".Rdata"))
    write.csv(beta_tables[[key]], csv_b, row.names = FALSE)
    # .Rdata object named betas_p<P> (legacy convention per sim/chpc/README.md)
    assign(paste0("betas_", key), beta_tables[[key]])
    save(list = paste0("betas_", key), file = rdata_b)
    written  <- c(written, csv_b, rdata_b)
    message(sprintf("  wrote %s  (%d rows)",
                    basename(csv_b), nrow(beta_tables[[key]])))

    # -- model_performance --
    if (!is.null(perf_tables[[key]])) {
      rownames(perf_tables[[key]]) <- NULL
      csv_mp <- file.path(out_dir,
                          sprintf("model_performance_%s_%s.csv", key, scenario_name))
      write.csv(perf_tables[[key]], csv_mp, row.names = FALSE)
      written <- c(written, csv_mp)
      message(sprintf("  wrote %s  (%d rows)",
                      basename(csv_mp), nrow(perf_tables[[key]])))
    }

    # -- by_block --
    if (!is.null(block_tables[[key]])) {
      rownames(block_tables[[key]]) <- NULL
      csv_bb <- file.path(out_dir,
                          sprintf("by_block_%s_%s.csv", key, scenario_name))
      write.csv(block_tables[[key]], csv_bb, row.names = FALSE)
      written <- c(written, csv_bb)
      message(sprintf("  wrote %s  (%d rows)",
                      basename(csv_bb), nrow(block_tables[[key]])))
    }

    # -- false_positives --
    if (!is.null(fp_tables[[key]]) && nrow(fp_tables[[key]]) > 0L) {
      rownames(fp_tables[[key]]) <- NULL
      csv_fp <- file.path(out_dir,
                          sprintf("false_positives_%s_%s.csv", key, scenario_name))
      write.csv(fp_tables[[key]], csv_fp, row.names = FALSE)
      written <- c(written, csv_fp)
      message(sprintf("  wrote %s  (%d rows)",
                      basename(csv_fp), nrow(fp_tables[[key]])))
    }
  }

  # C5: stacked model_performance across all p (only when >1 group)
  if (length(perf_tables) > 1L) {
    all_perf       <- do.call(rbind, perf_tables)
    rownames(all_perf) <- NULL
    csv_stack      <- file.path(out_dir,
                                sprintf("model_performance_%s.csv", scenario_name))
    write.csv(all_perf, csv_stack, row.names = FALSE)
    written <- c(written, csv_stack)
    message(sprintf("  wrote %s  (%d rows, all groups)",
                    basename(csv_stack), nrow(all_perf)))
  }

  # B7: copy scenario.env into out_dir
  if (has_manifest) {
    dst <- file.path(out_dir, "scenario.env")
    file.copy(manifest_path, dst, overwrite = TRUE)
    written <- c(written, dst)
    message(sprintf("  copied scenario.env -> %s", out_dir))
  }

  # ---- Summary ----
  if (n_failed > 0L) {
    failed_df <- do.call(rbind, lapply(failed_info, as.data.frame,
                                       stringsAsFactors = FALSE))
    # Distinct error messages with counts
    err_counts <- sort(table(failed_df$error), decreasing = TRUE)
    message(sprintf(
      "\n[parse] %d failure record(s) skipped. Run IDs per p:", n_failed))
    for (pv in sort(unique(failed_df$p))) {
      ids <- sort(failed_df$run_id[failed_df$p == pv])
      message(sprintf("  p=%d: runs %s", pv, paste(ids, collapse = ", ")))
    }
    message("  Distinct errors:")
    for (j in seq_along(err_counts))
      message(sprintf("    (x%d) %s", err_counts[[j]], names(err_counts)[j]))
  } else {
    failed_df <- data.frame(p = integer(0), run_id = integer(0),
                            error = character(0), stringsAsFactors = FALSE)
  }

  message(sprintf("\n[parse] Done. runs=%d  failed=%d  skipped_names=%d",
                  n_runs, n_failed, n_skipped_nm))
  message(sprintf("  Output: %s", out_dir))

  invisible(list(
    written        = written,
    n_runs         = n_runs,
    n_failed       = n_failed,
    n_skipped      = n_skipped_nm,
    failed_summary = failed_df
  ))
}
