#################### Parse Batch Run ################
# Cycles through Sim/Local_Testing/, pulls each fitted model out of its .Rdata
# file, and:
#   (1) assigns the model object to a variable named exactly like the file
#       (e.g. file n200_p25_run1.Rdata -> object `n200_p25_run1`), and
#   (2) builds one long tidy data frame per npredictors group, with 1 row per
#       predictor per run and columns:
#         predictor, True_beta, Estimate, nobs, run_number
#
# Each per-group table is assigned to a variable `betas_p25`, `betas_p50`, ...
# and all are collected in the named list `beta_tables`.
#
# Note: the fitted model lives in result$final_mod (cause-1 subdistribution
# hazard), and result$meta$beta1 holds the true (padded) cause-1 coefficients.


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

  # (3) Extract FG_True estimates
  fg_estimates <- rep(NA, p_f)
  active_idx <- which(true_beta != 0)
  fg_estimates[active_idx] <- res$fg_true$coef

  df <- data.frame(
    predictor  = coefs$Variable,
    True_beta  = true_beta,
    Estimate   = coefs$Estimate,
    FG_True    = fg_estimates,
    nobs       = nobs_f,
    run_number = run_f,
    stringsAsFactors = FALSE
  )

  key <- paste0("p", p_f)
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
