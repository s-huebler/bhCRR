#################### Parse Batch Run ################
# Cycles through Sim/Local_Testing/, pulls each fitted model out of its .Rdata
# file, and:
#   (1) assigns the model object to a variable named exactly like the file
#       (e.g. file n200_p25_run1.Rdata -> object `n200_p25_run1`), and
#   (2) builds one long tidy data frame per npredictors group, with 1 row per
#       predictor per model per run and columns:
#         predictor, Block, BlockRho, True, FG, Model, Estimate, Bias, nobs,
#         run_number
#   (3) builds a second, wide data frame with 1 row per model per run holding
#       that fit's tuning parameters, convergence status, iteration count and
#       selection confusion counts (TP/FP/TN/FN), both overall and per
#       correlation block when the run carries a block map. See the
#       "Per-run model performance" section at the bottom for the exact
#       column list.
#
# Block / BlockRho come from result$meta$block_map, which batch_run_complex.R
# saves per run. They are joined BY PREDICTOR NAME, never by position, and are
# added only when the run actually carries a block map, so older result files
# predating the block structure still parse (their Block is NA).
#
# Runs pooled into one group table must share result$meta$block_signature. The
# assignment is deterministic within a scenario, so a mismatch means the block
# specification itself changed partway through the batch and the per-block
# summaries would be pooling incomparable strata. That is reported, not averaged
# over.
#
# The Model column identifies which method each Estimate came from
#   (SSL, LASSO, aLASSO, SCAD, MCP), mirroring the mod_comparison_long table at
#   the bottom of single_run.R. True is the padded cause-1 coefficient, FG is the
#   low-dimensional Fine-Gray estimate fit on the true active set (NA elsewhere),
#   and Bias = Estimate - True.
#
# Each per-group table is assigned to a variable `betas_p25`, `betas_p50`, ...
# and all are collected in the named list `beta_tables`.
#
# Note: the SSL fit lives in result$final_mod (cause-1 subdistribution hazard),
# the BIC-selected comparators live in result$selected_models (each with a
# beta_raw vector in original predictor units), and result$meta$beta1 holds the
# true (padded) cause-1 coefficients.


repo_root <- "/Users/sophiehuebler/Documents/bhCRR"
in_dir    <- file.path(repo_root, "Sim/Local_Testing")

# Fallback true betas, only used if a file predates meta$beta1 being saved.
beta1_active_fallback <- c(0.40, -0.50, 0.60, 0.75, -0.80)


#################### Performance-table settings ################
# A predictor counts as "selected" when abs(Estimate) > selection_tol. Every
# method here returns exact zeros for unselected predictors, so 0 is the honest
# default; raise it only if you deliberately want a magnitude-based rule, and
# say so when you report the numbers.
selection_tol <- 0

# Concavity parameters are fixed in the batch script and are NOT written into
# the saved result objects, so they cannot be recovered from disk. These values
# are an ANNOTATION copied from batch_run_complex.R -- if you change gamma
# there, change it here too, or set the entries to NA_real_ rather than let the
# column assert something the run did not do.
comparator_gamma <- c(LASSO  = NA_real_,
                      aLASSO = NA_real_,
                      SCAD   = 3.7,
                      MCP    = 2.7)

# Where to write the per-run performance table. Set to NULL to skip writing.
perf_csv <- file.path(in_dir, "model_performance.csv")


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
files <- list.files(in_dir, pattern = "\\.Rdata$", full.names = TRUE)
if (length(files) == 0) stop("No .Rdata files found in ", in_dir)

# Parse n / p / run from each filename and sort naturally (p, then run).
meta_tbl <- regmatches(basename(files),
                       regexec("n(\\d+)_p(\\d+)_run(\\d+)", basename(files)))
meta_tbl <- do.call(rbind, lapply(meta_tbl, function(m) as.integer(m[2:4])))
colnames(meta_tbl) <- c("nobs", "npredictors", "run")
ord   <- order(meta_tbl[, "npredictors"], meta_tbl[, "run"])
files <- files[ord]
meta_tbl <- meta_tbl[ord, , drop = FALSE]


#################### Cycle through files ################
model_objects <- character(0)   # names of model objects created in .GlobalEnv
beta_tables   <- list()         # per-npredictors long data frames
perf_tables   <- list()         # per-npredictors model x run performance rows
skipped       <- character(0)   # files with no fitted model (failed runs)
block_sigs    <- list()         # block_signature seen per npredictors group
no_block_map  <- character(0)   # runs with no saved block map

for (i in seq_along(files)) {
  f  <- files[i]
  nm <- tools::file_path_sans_ext(basename(f))   # e.g. "n200_p25_run1"

  e <- new.env()
  load(f, envir = e)             # loads object `result`
  res <- e$result

  nobs_f <- meta_tbl[i, "nobs"]
  p_f    <- meta_tbl[i, "npredictors"]
  run_f  <- meta_tbl[i, "run"]

  # Failed run (batch_run saved only meta with an error) -> skip the table,
  # but still surface it.
  if (is.null(res$final_mod)) {
    skipped <- c(skipped, nm)
    warning(sprintf("%s has no final_mod (failed run) -- skipped.", nm))
    next
  }

  # (1) Expose the model object under a variable named exactly like the file.
  assign(nm, res$final_mod, envir = .GlobalEnv)
  model_objects <- c(model_objects, nm)

  # (2) Add this run's coefficients to its npredictors group.
  coefs     <- res$final_mod$coefficients
  true_beta <- if (!is.null(res$meta$beta1)) {
    res$meta$beta1
  } else {
    c(beta1_active_fallback, rep(0, p_f - length(beta1_active_fallback)))
  }

  # (3) Extract FG (low-dim Fine-Gray on the true active set) estimates. These
  # are shared across every model row for a given predictor/run.
  fg_estimates <- rep(NA_real_, p_f)
  active_idx <- which(true_beta != 0)
  fg_estimates[active_idx] <- res$fg_true$coef

  predictors <- coefs$Variable

  # (4) SSL estimates (from result$final_mod).
  ssl_df <- data.frame(
    predictor = predictors,
    True      = true_beta,
    FG        = fg_estimates,
    Model     = "SSL",
    Estimate  = coefs$Estimate,
    stringsAsFactors = FALSE
  )

  # (5) Comparator estimates (LASSO / aLASSO / SCAD / MCP) from selected_models.
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

  # (6) Stack SSL + comparators into one long table for this run.
  df <- rbind(ssl_df, comp_df)
  df$Bias       <- df$Estimate - df$True
  df$nobs       <- nobs_f
  df$run_number <- run_f

  key <- paste0("p", p_f)

  # (7) Correlation stratum, joined by predictor name. Runs written before the
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

  # (8) The block map must be constant across every run pooled into a group.
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

  df <- df[, c("predictor", "Block", "BlockRho", "True", "FG", "Model",
               "Estimate", "Bias", "nobs", "run_number")]

  beta_tables[[key]] <- rbind(beta_tables[[key]], df)

  # (9) One performance row per model for this run: tuning parameters,
  # convergence, iterations, and the selection confusion counts. Tuning
  # parameters that do not apply to a method stay NA (s0/s1 are SSL-only;
  # lambda is comparator-only; gamma is SCAD/MCP-only).
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


#################### Expose per-group tables ################
# Assign each group table to its own variable (betas_p25, betas_p50, ...).
for (key in names(beta_tables)) {
  assign(paste0("betas_", key), beta_tables[[key]], envir = .GlobalEnv)
}

message(sprintf("Loaded %d model objects: %s",
                length(model_objects),
                paste(model_objects, collapse = ", ")))
message(sprintf("Built %d group table(s): %s",
                length(beta_tables),
                paste(paste0("betas_", names(beta_tables)), collapse = ", ")))
if (length(skipped) > 0)
  message(sprintf("Skipped %d failed run(s): %s",
                  length(skipped), paste(skipped, collapse = ", ")))
if (length(no_block_map) > 0)
  message(sprintf("%d run(s) carried no block map (Block is NA): %s",
                  length(no_block_map), paste(no_block_map, collapse = ", ")))
for (key in names(block_sigs)) {
  message(sprintf("Block assignment for %s: %s", key, block_sigs[[key]]))
}


#################### Per-run model performance ################
# One row per model per run:
#   run_id, nobs, npredictors, run_number, Model
#   s0, s1        -- SSL spike/slab scales (NA for every other method)
#   lambda        -- BIC-selected shrinkage for LASSO/aLASSO/SCAD/MCP (NA for SSL)
#   gamma         -- concavity for SCAD/MCP (NA elsewhere); see comparator_gamma
#   converged     -- did the fit converge
#   iterations    -- EM iterations for SSL; NA for the comparators (not recorded)
#   TP/FP/TN/FN   -- selection confusion counts against the true cause-1 support
#   TP_<block>, FP_<block>, TN_<block>, FN_<block>
#                 -- the same counts within each correlation stratum, present
#                    only for runs that carried a block map
#
# Each group is exposed as perf_p25, perf_p50, ... and all groups are stacked
# into model_performance.
for (key in names(perf_tables)) {
  rownames(perf_tables[[key]]) <- NULL
  assign(paste0("perf_", key), perf_tables[[key]], envir = .GlobalEnv)
}

model_performance <- NULL
for (key in names(perf_tables)) {
  model_performance <- rbind_fill(model_performance, perf_tables[[key]])
}
if (!is.null(model_performance)) rownames(model_performance) <- NULL
assign("model_performance", model_performance, envir = .GlobalEnv)

if (!is.null(model_performance)) {
  message(sprintf("Built performance table: %d row(s) across %d group(s): %s",
                  nrow(model_performance), length(perf_tables),
                  paste(paste0("perf_", names(perf_tables)), collapse = ", ")))

  if (!is.null(perf_csv)) {
    utils::write.csv(model_performance, perf_csv, row.names = FALSE)
    message(sprintf("  -> wrote %s", perf_csv))
  }
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
for (key in names(fp_tables)) {
  assign(paste0("fp_", key), fp_tables[[key]], envir = .GlobalEnv)
}
if (length(fp_tables) > 0) {
  message(sprintf("Built %d false-positive table(s): %s",
                  length(fp_tables),
                  paste(paste0("fp_", names(fp_tables)), collapse = ", ")))
}
