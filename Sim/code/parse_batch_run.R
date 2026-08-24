#################### Parse Batch Run ################
# Cycles through Sim/Local_Testing/, pulls each fitted model out of its .Rdata
# file, and:
#   (1) assigns the model object to a variable named exactly like the file
#       (e.g. file n200_p25_run1.Rdata -> object `n200_p25_run1`), and
#   (2) builds one long tidy data frame per npredictors group, with 1 row per
#       predictor per model per run and columns:
#         predictor, Block, BlockRho, True, FG, Model, Estimate, Bias, nobs,
#         run_number
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
