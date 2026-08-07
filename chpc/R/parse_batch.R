#################### CHPC Parse Batch Run ################
# CHPC/SLURM-portable adaptation of Sim/code/parse_batch_run.R.
#
# Reads every n<nobs>_p<p>_run<run>.Rdata in $SCRATCH_RUN_DIR, groups by
# npredictors, and builds one long tidy table per group with columns:
#   predictor, True, FG, Model, Estimate, Bias, nobs, run_number
# mirroring the original parse. Each group table is written to
# Sim/chpc_results/ as both a CSV and an .Rdata (holding an object named
# betas_p<p>) so you can pull them to local and make figures.
#
#   True  = padded cause-1 coefficient (meta$beta1)
#   FG    = low-dim Fine-Gray estimate on the true active set (NA elsewhere)
#   Model = SSL / LASSO_BIC / LASSO_Wolbers / aLASSO_* / SCAD_* / MCP_*
#   Bias  = Estimate - True
#
# Config comes from the environment (set by config.sh / parse.slurm):
#   REPO_ROOT, SCRATCH_RUN_DIR, NOBS


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
out_dir   <- file.path(repo_root, "Sim", "chpc_results")
nobs_env  <- as.integer(get_env("NOBS", "NA"))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Fallback true betas, only used if a file predates meta$beta1 being saved.
beta1_active_fallback <- c(0.40, -0.50, 0.60, 0.75, -0.80)

message(sprintf("parse_batch.R  in=%s", in_dir))
message(sprintf("  out=%s", out_dir))


#################### Locate & order files ################
files <- list.files(in_dir, pattern = "\\.Rdata$", full.names = TRUE)
if (length(files) == 0) stop("No .Rdata files found in ", in_dir)

meta_tbl <- regmatches(basename(files),
                       regexec("n(\\d+)_p(\\d+)_run(\\d+)", basename(files)))
meta_tbl <- do.call(rbind, lapply(meta_tbl, function(m) as.integer(m[2:4])))
colnames(meta_tbl) <- c("nobs", "npredictors", "run")
ord   <- order(meta_tbl[, "npredictors"], meta_tbl[, "run"])
files <- files[ord]
meta_tbl <- meta_tbl[ord, , drop = FALSE]


#################### Cycle through files ################
beta_tables <- list()         # per-npredictors long data frames
skipped     <- character(0)   # files with no fitted model (failed runs)

for (i in seq_along(files)) {
  f  <- files[i]
  nm <- tools::file_path_sans_ext(basename(f))   # e.g. "n200_p25_run1"

  e <- new.env()
  load(f, envir = e)             # loads object `result`
  res <- e$result

  nobs_f <- meta_tbl[i, "nobs"]
  p_f    <- meta_tbl[i, "npredictors"]
  run_f  <- meta_tbl[i, "run"]

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

  fg_estimates <- rep(NA_real_, p_f)
  active_idx <- which(true_beta != 0)
  fg_estimates[active_idx] <- res$fg_true$coef

  predictors <- coefs$Variable

  ssl_df <- data.frame(
    predictor = predictors,
    True      = true_beta,
    FG        = fg_estimates,
    Model     = "SSL",
    Estimate  = coefs$Estimate,
    stringsAsFactors = FALSE
  )

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
  df <- df[, c("predictor", "True", "FG", "Model",
               "Estimate", "Bias", "nobs", "run_number")]

  key <- paste0("p", p_f)
  beta_tables[[key]] <- rbind(beta_tables[[key]], df)
}


#################### Write per-group tables ################
written <- character(0)
for (key in names(beta_tables)) {
  tbl <- beta_tables[[key]]

  # nobs tag for the filename: use the value seen in the data (should be uniform
  # within a per-nobs scratch dir), falling back to the NOBS env var.
  ntag <- if (length(unique(tbl$nobs)) == 1L) unique(tbl$nobs) else nobs_env
  base <- sprintf("betas_%s_n%s", key, ntag)

  csv_path  <- file.path(out_dir, paste0(base, ".csv"))
  rdata_path<- file.path(out_dir, paste0(base, ".Rdata"))

  write.csv(tbl, csv_path, row.names = FALSE)

  # Save under the object name betas_p<p> so load() gives you the familiar name.
  assign(paste0("betas_", key), tbl)
  save(list = paste0("betas_", key), file = rdata_path)

  written <- c(written, csv_path)
  message(sprintf("  wrote %s  (%d rows, runs: %s)",
                  basename(csv_path), nrow(tbl),
                  paste(range(tbl$run_number), collapse = "-")))
}

message(sprintf("Built %d group table(s): %s",
                length(beta_tables),
                paste(paste0("betas_", names(beta_tables)), collapse = ", ")))
if (length(skipped) > 0)
  message(sprintf("Skipped %d failed run(s): %s",
                  length(skipped), paste(skipped, collapse = ", ")))
message("Parse complete.")
