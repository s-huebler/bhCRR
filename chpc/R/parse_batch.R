#################### CHPC Parse Batch Run ################
# CHPC/SLURM-portable adaptation of Sim/code/parse_batch_run.R.
#
# Reads every p<P>_<scenario>_run<N>.Rdata under $SCRATCH_RUN_DIR, which now
# points at one scenario directory containing p<P> subdirs. Files that do not
# match the ^p<P>_<scenario>_run<N>.Rdata$ pattern are skipped with a message
# (this includes old n<N>_p<P>_run<N> files); the parse errors out only if
# zero files match the pattern.
#
# All outputs go to Sim/chpc_results/<scenario>/. Per-group tables are produced
# for each npredictors value (grouping key stays "p<P>"):
#
#   betas_p<P>_<scenario>.csv / .Rdata
#     Long tidy table, one row per predictor per model per run:
#       scenario, predictor, Block, BlockRho, True, FG, Model, Estimate,
#       Bias, nobs, run_number
#     Block / BlockRho come from result$meta$block_map, joined by predictor
#     name. Runs without a block map have Block = NA.
#
#   model_performance_p<P>_<scenario>.csv
#     One row per model per run: scenario, run_id, nobs, npredictors,
#     run_number, tuning params (s0/s1 for SSL, lambda for comparators),
#     convergence, iteration count, and TP/FP/TN/FN overall and per
#     correlation stratum when a block map is present.
#
#   false_positives_p<P>_<scenario>.csv
#     Null-predictor false-positive rates by Block x Model -- the key
#     diagnostic for selection leakage across correlated neighbours.
#
# If more than one npredictors group is present, model_performance_<scenario>.csv
# (stacked across groups) is also written.
#
# Checks performed:
#   - Filename scenario token must equal SCENARIO (hard stop on mismatch).
#   - nobs is read from res$meta$nobs (not the filename) and must be identical
#     across every run in the scenario; hard stop if not, naming both values.
#   - result$meta$block_signature must be identical within each p group.
#   - result$meta$scenario_params is compared field-by-field (all.equal,
#     tolerance 1e-8) against the first run's values; hard stop on any
#     difference, naming the field, both values, and the run_id. Runs saved
#     before scenario_params existed have it NULL -- warned and skipped.
#
# $SCRATCH_RUN_DIR/scenario.env is copied into out_dir if present, making the
# results folder self-describing without the scratch tree.
#
# Config from the environment (set by config.sh / parse.slurm):
#   REPO_ROOT, SCRATCH_RUN_DIR, SCENARIO


#################### Config from environment ################
get_env <- function(name, default = NULL) {
  v <- Sys.getenv(name, unset = NA_character_)
  if (is.na(v) || v == "") {
    if (is.null(default)) stop("Required environment variable not set: ", name)
    return(default)
  }
  v
}

repo_root <- get_env("REPO_ROOT")
in_dir    <- get_env("SCRATCH_RUN_DIR")
scenario  <- get_env("SCENARIO")
out_dir   <- file.path(repo_root, "Sim", "chpc_results", scenario)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Fallback true betas, only used if a file predates meta$beta1 being saved.
beta1_active_fallback <- c(0.40, -0.50, 0.60, 0.75, -0.80)

message(sprintf("parse_batch.R  scenario=%s", scenario))
message(sprintf("  in=%s", in_dir))
message(sprintf("  out=%s", out_dir))


#################### Performance-table settings ################
# A predictor counts as "selected" when abs(Estimate) > selection_tol. Every
# method here returns exact zeros for unselected predictors, so 0 is the honest
# default; raise it only if you deliberately want a magnitude-based rule, and
# say so when you report the numbers.
selection_tol <- 0

# Concavity parameters are fixed in the batch script and are NOT written into
# the saved result objects, so they cannot be recovered from disk. These values
# are an ANNOTATION copied from chpc/R/run_batch.R -- if you change gamma
# there, change it here too, or set entries to NA_real_ rather than let the
# column assert something the run did not do.
# LASSO-ONLY: was c(LASSO=NA, aLASSO=NA, SCAD=3.7, MCP=2.7); keep in sync with run_batch.R
comparator_gamma <- c(LASSO = NA_real_)


#################### Small helpers ################
# base R gained %||% in 4.4.0; define it only if this session predates that.
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}

# Confusion counts for one fitted model against the true cause-1 coefficients.
# "Positive" means the method put a non-zero coefficient on the predictor;
# "true" means the data-generating coefficient was non-zero.
confusion_counts <- function(true_beta, estimate, tol = 0) {
  keep <- !is.na(estimate) & !is.na(true_beta)
  sel  <- keep & abs(estimate) > tol
  act  <- keep & true_beta != 0
  c(TP = sum(sel & act),
    FP = sum(sel & !act),
    TN = sum(!sel & !act & keep),
    FN = sum(!sel & act))
}

# One row of confusion counts: overall, then one set per correlation block.
# block_levels is empty for runs with no block map, which simply drops the
# per-block columns rather than inventing strata.
performance_counts <- function(sub, block_levels, tol = 0) {
  out <- as.list(confusion_counts(sub$True, sub$Estimate, tol))
  for (b in block_levels) {
    idx <- !is.na(sub$Block) & sub$Block == b
    cnt <- confusion_counts(sub$True[idx], sub$Estimate[idx], tol)
    names(cnt) <- paste0(names(cnt), "_", b)
    out <- c(out, as.list(cnt))
  }
  out
}

# rbind that tolerates a differing column set (e.g. a run with no block map
# stacked onto runs that have one). Missing cells become NA instead of an error.
rbind_fill <- function(a, b) {
  if (is.null(a)) return(b)
  for (cc in setdiff(names(a), names(b))) b[[cc]] <- NA
  for (cc in setdiff(names(b), names(a))) a[[cc]] <- NA
  rbind(a, b[, names(a), drop = FALSE])
}


#################### Locate & order files ################
# recursive = TRUE: run_batch.R writes into per-npredictors p<P> subdirectories.
# Grouping is driven by the p in each FILENAME. Point in_dir at the scenario
# root to parse every p, or at one p<P> subdir to scope it.
all_files <- list.files(in_dir, pattern = "\\.Rdata$", full.names = TRUE,
                        recursive = TRUE)
if (length(all_files) == 0) {
  stop("No .Rdata files found under ", in_dir,
       " (searched recursively, including any p<P> subdirectories).")
}

# Match only the new p<P>_<scenario>_run<N>.Rdata format. Files that do not
# match (including old n<N>_p<P>_run<N>.Rdata files) are skipped.
rx <- "^p(\\d+)_([A-Za-z0-9][A-Za-z0-9.-]*)_run(\\d+)\\.Rdata$"
rx_matches <- regmatches(basename(all_files), regexec(rx, basename(all_files)))
matched    <- vapply(rx_matches, length, integer(1)) > 0

if (any(!matched)) {
  skipped_names <- basename(all_files)[!matched]
  message(sprintf("Skipping %d file(s) that do not match p<P>_<scenario>_run<N>.Rdata:\n  %s",
                  length(skipped_names),
                  paste(skipped_names, collapse = "\n  ")))
}
if (!any(matched)) {
  stop("No files matching p<P>_<scenario>_run<N>.Rdata found under ", in_dir, ".")
}

files      <- all_files[matched]
rx_matches <- rx_matches[matched]

# Verify every file's scenario token against SCENARIO before loading anything.
file_scenarios <- vapply(rx_matches, function(m) m[[3L]], character(1))
bad_sc <- file_scenarios != scenario
if (any(bad_sc)) {
  stop(sprintf(
    "Scenario token mismatch in %d file(s) -- expected '%s':\n%s",
    sum(bad_sc), scenario,
    paste(sprintf("  %s  (token: '%s')", basename(files)[bad_sc],
                  file_scenarios[bad_sc]),
          collapse = "\n")))
}

meta_tbl <- data.frame(
  npredictors = as.integer(vapply(rx_matches, function(m) m[[2L]], character(1))),
  run         = as.integer(vapply(rx_matches, function(m) m[[4L]], character(1))),
  stringsAsFactors = FALSE
)
ord      <- order(meta_tbl$npredictors, meta_tbl$run)
files    <- files[ord]
meta_tbl <- meta_tbl[ord, , drop = FALSE]


#################### Cycle through files ################
beta_tables         <- list()          # per-npredictors long data frames
perf_tables         <- list()          # per-npredictors model x run performance rows
skipped             <- character(0)    # files with no fitted model (failed runs)
block_sigs          <- list()          # block_signature seen per npredictors group
no_block_map        <- character(0)    # runs with no saved block map
nobs_check          <- NULL            # nobs seen on first loaded run
scenario_params_ref <- NULL            # scenario_params from first run with the field

for (i in seq_along(files)) {
  f  <- files[i]
  nm <- tools::file_path_sans_ext(basename(f))   # e.g. "p25_n200base_run1"

  e <- new.env()
  load(f, envir = e)             # loads object `result`
  res <- e$result

  p_f   <- meta_tbl$npredictors[i]
  run_f <- meta_tbl$run[i]

  # nobs from meta (no longer encoded in the filename).
  nobs_f <- as.integer(res$meta$nobs)
  if (is.null(nobs_f) || is.na(nobs_f)) {
    warning(sprintf("%s has no meta$nobs.", nm))
    nobs_f <- NA_integer_
  } else if (is.null(nobs_check)) {
    nobs_check <- nobs_f
  } else if (!identical(nobs_check, nobs_f)) {
    stop(sprintf(
      "nobs mismatch across runs in scenario '%s': '%s' has nobs=%d but earlier runs have nobs=%d.",
      scenario, nm, nobs_f, nobs_check))
  }

  # scenario_params consistency (field-by-field; skip for pre-field runs).
  sp <- res$meta$scenario_params
  if (is.null(sp)) {
    warning(sprintf(
      "%s has no meta$scenario_params (saved before this field existed) -- skipping scenario-param check.",
      nm))
  } else if (is.null(scenario_params_ref)) {
    scenario_params_ref <- sp
  } else {
    for (fld in names(scenario_params_ref)) {
      if (!fld %in% names(sp)) {
        stop(sprintf(
          "scenario_params mismatch in run '%s': field '%s' is present in the reference run but missing here.",
          nm, fld))
      }
      cmp <- all.equal(scenario_params_ref[[fld]], sp[[fld]], tolerance = 1e-8)
      if (!isTRUE(cmp)) {
        stop(sprintf(
          "scenario_params mismatch in run '%s', field '%s':\n  reference: %s\n  this run:  %s",
          nm, fld,
          paste(scenario_params_ref[[fld]], collapse = ", "),
          paste(sp[[fld]], collapse = ", ")))
      }
    }
  }

  # Failed run (batch saved only meta with an error) -> skip the table.
  if (is.null(res$final_mod)) {
    skipped <- c(skipped, nm)
    warning(sprintf("%s has no final_mod (failed run) -- skipped.", nm))
    next
  }

  coefs     <- res$final_mod$coefficients
  true_beta <- if (!is.null(res$meta$beta1)) {
    res$meta$beta1
  } else {
    c(beta1_active_fallback, rep(0, p_f - length(beta1_active_fallback)))
  }

  # FG (low-dim Fine-Gray on the true active set) estimates. Shared across all
  # model rows for a given predictor/run.
  fg_estimates <- rep(NA_real_, p_f)
  active_idx <- which(true_beta != 0)
  fg_estimates[active_idx] <- res$fg_true$coef

  predictors <- coefs$Variable

  # SSL estimates (from result$final_mod).
  ssl_df <- data.frame(
    predictor = predictors,
    True      = true_beta,
    FG        = fg_estimates,
    Model     = "SSL",
    Estimate  = coefs$Estimate,
    stringsAsFactors = FALSE
  )

  # Comparator estimates from selected_models.
  # beta_raw is named by predictor and already on the original coefficient scale.
  comp_df <- NULL
  if (!is.null(res$selected_models)) {
    comp_df <- do.call(rbind, lapply(names(res$selected_models), function(mn) {
      beta_raw <- res$selected_models[[mn]]$beta_raw
      data.frame(
        predictor = predictors,
        True      = true_beta,
        FG        = fg_estimates,
        Model     = mn,
        Estimate  = as.numeric(beta_raw[predictors]),
        stringsAsFactors = FALSE
      )
    }))
  } else {
    warning(sprintf("%s has no selected_models -- only SSL estimates recorded.", nm))
  }

  df <- rbind(ssl_df, comp_df)
  df$Bias       <- df$Estimate - df$True
  df$nobs       <- nobs_f
  df$run_number <- run_f

  key <- paste0("p", p_f)

  # Correlation stratum, joined by predictor name. Runs written before the
  # block structure existed have no map, so their Block stays NA rather than
  # being guessed at.
  bm <- res$meta$block_map
  if (!is.null(bm)) {
    hit <- match(df$predictor, bm$name)
    if (anyNA(hit)) {
      warning(sprintf("%s: %d predictor(s) absent from block_map -- Block set NA.",
                      nm, sum(is.na(hit))))
    }
    df$Block    <- bm$block[hit]
    df$BlockRho <- bm$rho[hit]
  } else {
    no_block_map <- c(no_block_map, nm)
    df$Block    <- NA_character_
    df$BlockRho <- NA_real_
  }

  # The block map must be constant across every run pooled into a group.
  # It is deterministic given the scenario, so a mismatch means the block spec
  # changed mid-batch and these runs cannot be summarised together.
  sig <- res$meta$block_signature
  if (!is.null(sig)) {
    if (is.null(block_sigs[[key]])) {
      block_sigs[[key]] <- sig
    } else if (!identical(block_sigs[[key]], sig)) {
      stop(sprintf(paste0("%s has a different block assignment from earlier ",
                          "runs in group %s.\n  earlier: %s\n  this run: %s\n",
                          "Per-block summaries would pool incomparable strata. ",
                          "Separate these runs or re-run the batch."),
                   nm, key, block_sigs[[key]], sig))
    }
  }

  df$scenario <- scenario
  df <- df[, c("scenario", "predictor", "Block", "BlockRho", "True", "FG",
               "Model", "Estimate", "Bias", "nobs", "run_number")]

  beta_tables[[key]] <- rbind(beta_tables[[key]], df)

  # One performance row per model for this run: tuning parameters, convergence,
  # iterations, and selection confusion counts. Tuning parameters that do not
  # apply to a method stay NA (s0/s1 are SSL-only; lambda is comparator-only;
  # gamma is SCAD/MCP-only).
  block_levels <- if (!is.null(bm)) sort(unique(as.character(bm$block))) else character(0)

  for (mn in unique(df$Model)) {
    sub <- df[df$Model == mn, , drop = FALSE]

    if (identical(mn, "SSL")) {
      # Tuned (s0, s1) from bhcrr_autotune, with the fitted object's own ss as
      # a fallback for runs saved before result$best existed.
      s0_v   <- res$best$s0 %||% res$final_mod$ss[1]
      s1_v   <- res$best$s1 %||% res$final_mod$ss[2]
      lam_v  <- NA_real_
      gam_v  <- NA_real_
      conv_v <- if (is.null(res$final_mod$conv)) NA else isTRUE(res$final_mod$conv)
      # EM iterations completed. fit_ssl_psdh() increments this at the END of
      # each pass and breaks before incrementing on the converging pass, so it
      # is one less than the number of EM passes actually run.
      iter_v <- res$final_mod$iterations %||% NA_integer_
    } else {
      sm     <- res$selected_models[[mn]]
      s0_v   <- NA_real_
      s1_v   <- NA_real_
      lam_v  <- sm$lambda %||% NA_real_
      gam_v  <- unname(comparator_gamma[mn])
      conv_v <- if (is.null(sm$converged)) NA else as.logical(sm$converged)
      # fastCrrp() reports convergence per lambda but not an iteration count,
      # and select_by_bic() does not record one, so this is NA for the
      # comparators unless the batch script starts saving it.
      iter_v <- sm$iterations %||% NA_integer_
    }

    row <- data.frame(
      scenario    = scenario,
      run_id      = nm,
      nobs        = nobs_f,
      npredictors = p_f,
      run_number  = run_f,
      Model       = mn,
      s0          = as.numeric(s0_v),
      s1          = as.numeric(s1_v),
      lambda      = as.numeric(lam_v),
      gamma       = as.numeric(gam_v),
      converged   = as.logical(conv_v),
      iterations  = as.integer(iter_v),
      stringsAsFactors = FALSE
    )
    row <- cbind(row, as.data.frame(performance_counts(sub, block_levels,
                                                       selection_tol)))
    rownames(row) <- NULL

    perf_tables[[key]] <- rbind_fill(perf_tables[[key]], row)
  }
}

if (length(no_block_map) > 0)
  message(sprintf("%d run(s) carried no block map (Block is NA): %s",
                  length(no_block_map), paste(no_block_map, collapse = ", ")))
for (key in names(block_sigs)) {
  message(sprintf("Block assignment for %s: %s", key, block_sigs[[key]]))
}


#################### Per-block false positives ################
# The diagnostic the block structure exists for: among predictors whose true
# cause-1 coefficient is 0, how often does each method give a non-zero estimate,
# split by correlation stratum. A strong-block rate well above the independent
# one is selection leaking across correlated neighbours.
false_positive_by_block <- function(df) {
  fp <- df[df$True == 0 & !is.na(df$Block), ]
  if (nrow(fp) == 0) return(NULL)
  agg <- aggregate(list(n = fp$Estimate != 0),
                   by = list(Block = fp$Block, BlockRho = fp$BlockRho,
                             Model = fp$Model),
                   FUN = function(z) c(rate = mean(z), count = sum(z),
                                       total = length(z)))
  out <- do.call(data.frame, agg)
  names(out) <- c("Block", "BlockRho", "Model", "fp_rate", "n_false_pos",
                  "n_null_predictors")
  out[order(out$Model, out$BlockRho), ]
}

fp_tables <- lapply(beta_tables, false_positive_by_block)
fp_tables <- fp_tables[!vapply(fp_tables, is.null, logical(1))]
for (k in names(fp_tables)) {
  fp_tables[[k]] <- cbind(scenario = scenario, fp_tables[[k]],
                          stringsAsFactors = FALSE)
}


#################### Write per-group tables ################
written <- character(0)
for (key in names(beta_tables)) {
  tbl <- beta_tables[[key]]

  # ---- betas ----
  base       <- sprintf("betas_%s_%s", key, scenario)
  csv_path   <- file.path(out_dir, paste0(base, ".csv"))
  rdata_path <- file.path(out_dir, paste0(base, ".Rdata"))

  write.csv(tbl, csv_path, row.names = FALSE)

  # Save under the object name betas_p<p> so load() gives you the familiar name.
  assign(paste0("betas_", key), tbl)
  save(list = paste0("betas_", key), file = rdata_path)

  written <- c(written, csv_path, rdata_path)
  message(sprintf("  wrote %s  (%d rows, runs: %s)",
                  basename(csv_path), nrow(tbl),
                  paste(range(tbl$run_number), collapse = "-")))

  # ---- model performance ----
  if (!is.null(perf_tables[[key]])) {
    rownames(perf_tables[[key]]) <- NULL
    perf_path <- file.path(out_dir,
                           sprintf("model_performance_%s_%s.csv", key, scenario))
    write.csv(perf_tables[[key]], perf_path, row.names = FALSE)
    written <- c(written, perf_path)
    message(sprintf("  wrote %s  (%d rows)",
                    basename(perf_path), nrow(perf_tables[[key]])))
  }

  # ---- false positives ----
  if (!is.null(fp_tables[[key]])) {
    fp_path <- file.path(out_dir,
                         sprintf("false_positives_%s_%s.csv", key, scenario))
    write.csv(fp_tables[[key]], fp_path, row.names = FALSE)
    written <- c(written, fp_path)
    message(sprintf("  wrote %s  (%d rows)",
                    basename(fp_path), nrow(fp_tables[[key]])))
  }
}


#################### Stacked performance across groups ################
# If more than one npredictors group is present, write a combined performance
# table so figures that compare across p don't need to load and rbind manually.
if (length(perf_tables) > 1L) {
  all_perf <- NULL
  for (key in names(perf_tables)) {
    all_perf <- rbind_fill(all_perf, perf_tables[[key]])
  }
  if (!is.null(all_perf)) {
    rownames(all_perf) <- NULL
    stack_path <- file.path(out_dir,
                            sprintf("model_performance_%s.csv", scenario))
    write.csv(all_perf, stack_path, row.names = FALSE)
    written <- c(written, stack_path)
    message(sprintf("  wrote %s  (%d rows, all groups)",
                    basename(stack_path), nrow(all_perf)))
  }
}


#################### Copy scenario.env ################
scenario_env_src <- file.path(in_dir, "scenario.env")
if (file.exists(scenario_env_src)) {
  scenario_env_dst <- file.path(out_dir, "scenario.env")
  file.copy(scenario_env_src, scenario_env_dst, overwrite = TRUE)
  written <- c(written, scenario_env_dst)
  message(sprintf("  copied scenario.env -> %s", out_dir))
}


#################### Summary ################
message(sprintf("Built %d group table(s): %s",
                length(beta_tables),
                paste(paste0("betas_", names(beta_tables)), collapse = ", ")))
if (length(perf_tables) > 0)
  message(sprintf("Built %d performance table(s): %s",
                  length(perf_tables),
                  paste(paste0("perf_", names(perf_tables)), collapse = ", ")))
if (length(fp_tables) > 0)
  message(sprintf("Built %d false-positive table(s): %s",
                  length(fp_tables),
                  paste(paste0("fp_", names(fp_tables)), collapse = ", ")))
if (length(skipped) > 0)
  message(sprintf("Skipped %d failed run(s): %s",
                  length(skipped), paste(skipped, collapse = ", ")))
message(sprintf("Wrote %d file(s) to %s", length(written), out_dir))
message("Parse complete.")
