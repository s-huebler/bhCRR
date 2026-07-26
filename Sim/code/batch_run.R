#################### Batch Runner ################
# Runs the single_run workflow repeatedly across a set of npredictors values,
# saving one result file per run to Sim/Local_Testing/.
#
# NOTE: set.seed() is intentionally NOT called anywhere. Each call to
# simulateTwoCauseFineGrayModel() draws fresh data, so every run produces a
# different dataset and therefore a different fitted model.
#
# Files are named like: n200_p25_run1.Rdata
# Each file contains a single object `result` (a list); load with:
#   load("Sim/Local_Testing/n200_p25_run1.Rdata"); str(result, max.level = 1)


#################### Set Up (load everything once) ################
library(tidyverse)
library(fastcmprsk, lib.loc = "/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library")
library(survival)

repo_root <- "/Users/sophiehuebler/Documents/bhCRR"

source(file.path(repo_root, "R/helpers.r"))
source(file.path(repo_root, "R/utils.r"))
source(file.path(repo_root, "R/dedupe_warnings.r"))

source(file.path(repo_root, "R/update_betas.r"))
source(file.path(repo_root, "R/expected_inclusion_probs.r"))
source(file.path(repo_root, "R/expected_penalty_weights.r"))
source(file.path(repo_root, "R/update_mixture_prob.r"))

Rcpp::sourceCpp(file.path(repo_root, "src/cv_fastcrrp.cpp"))
Rcpp::sourceCpp(file.path(repo_root, "src/RcppExports.cpp"))
source(file.path(repo_root, "R/cv_fastCrrp_cpp.r"))
source(file.path(repo_root, "R/fit_ssl_psdh.r"))

source(file.path(repo_root, "R/generate_foldid.r"))
source(file.path(repo_root, "R/predict_from_ssl_psdh.r"))
source(file.path(repo_root, "R/wolbers_c.r"))
source(file.path(repo_root, "R/cv_ssl_psdh.r"))
source(file.path(repo_root, "R/tune_ssl_psdh.r"))
source(file.path(repo_root, "R/threshold.R"))
source(file.path(repo_root, "R/tuning_diagnostics.r"))

# Start from a clean, non-deterministic RNG state. Remove any existing global
# seed so each run draws fresh data (belt-and-suspenders with the seed-leak
# fixes now in the CV/sim functions).
if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)


#################### Batch Configuration ################
# ---- edit these to control the sweep ----
nobs             <- 200                 # sample size (fixed across runs)
npredictors_grid <- c(100)      # set of npredictors to sweep over
run_start        <- 4           # e.g., set to 1 for first batch, 11 for second
run_end          <- 10           # e.g., set to 10 for first batch, 20 for second

# Active (nonzero) coefficients, padded with zeros to reach npredictors.
beta1_active <- c(0.40, -0.50, 0.60, 0.75, -0.80)
beta2_active <- c(0,     0.3,  0,    0,   -0.2)

zero_gap_target <- 0.1   # clinically relevant minimum treatment effect

out_dir <- file.path(repo_root, "Sim/Local_Testing")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


#################### One Run ################
# Generates fresh data, fits the initial model, autotunes, refits, and returns
# a bundled result list. No set.seed() here on purpose.
run_once <- function(nobs, npredictors, beta1_active, beta2_active,
                     zero_gap_target) {

  ## ----- Data generation -----
  beta1 <- c(beta1_active, rep(0, npredictors - length(beta1_active)))
  beta2 <- c(beta2_active, rep(0, npredictors - length(beta2_active)))

  sim <- simulateTwoCauseFineGrayModel(nobs = nobs,
                                       beta1 = beta1,
                                       beta2 = beta2,
                                       X = NULL,
                                       u.min = 100,
                                       u.max = 100,
                                       p = 0.5,
                                       returnX = TRUE)

  sim_data <- cbind(data.frame(ID = 1:nobs,
                               TTE = sim$ftime,
                               Status = as.integer(sim$fstatus)),
                    sim$X)

  names(sim_data) <- c("ID", "TTE", "Status", paste0("X_", 1:(ncol(sim_data) - 3)))

  x <- as.matrix(sim_data %>% select(starts_with("X")))
  y <- as.matrix(sim_data %>%
                   mutate(Status = as.numeric(Status)) %>%
                   select(TTE, Status))

  ## ----- FG model -----
  active_idx <- which(beta1 != 0)
  active_x <- x[, active_idx, drop = FALSE]

  fg_true <- fastcmprsk::fastCrr(fastcmprsk::Crisk(sim$ftime, as.integer(sim$fstatus)) ~ active_x,
                                 variance = FALSE)

  ## ----- SSL model -----
  ssl_psdh_time1 <- Sys.time()

  # Initial model fit
  mod <- fit_ssl_psdh(x, y,
                      ss = c(0.04, 0.6),
                      initial_sparsity = 0.05,
                      theta_a = 1,
                      theta_b = 1,
                      maxit = 50,
                      epsilon = 1e-04,
                      init = NULL,
                      init_lam_path = 10^seq(log10(0.1), log10(0.001), length = 10),
                      inner_maxit_start = 1000)

  # Tuning
  autotune <- bhcrr_autotune(mod,
                             beta_min   = zero_gap_target,
                             beta_floor = 0.01,
                             nfolds     = 10,
                             ncv        = 2,
                             foldid     = NULL,
                             max_widen  = 1)

  best_s0 <- autotune$best$s0
  best_s1 <- autotune$best$s1

  final_mod <- fit_ssl_psdh(x, y,
                            ss = c(best_s0, best_s1),
                            initial_sparsity = 0.05,
                            theta_a = 1,
                            theta_b = 1,
                            maxit = 100,
                            epsilon = 1e-04,
                            init = NULL,
                            init_lam_path = 10^seq(log10(0.1), log10(0.001), length = 10),
                            inner_maxit_start = 1000)

  ssl_psdh_time2 <- Sys.time()
  ssl_psdh_time  <- ssl_psdh_time2 - ssl_psdh_time1

  ## ----- Bundle everything needed for downstream diagnostics -----
  list(
    meta = list(nobs = nobs,
                npredictors = npredictors,
                beta1 = beta1,
                beta2 = beta2,
                zero_gap_target = zero_gap_target,
                timestamp = Sys.time(),
                ssl_psdh_time = ssl_psdh_time),
    sim_data  = sim_data,
    x         = x,
    y         = y,
    mod       = mod,
    autotune  = autotune,
    best      = list(s0 = best_s0, s1 = best_s1),
    final_mod = final_mod,
    fg_true   = fg_true    # <-- Add this line
  )
}


#################### Sweep ################
for (npredictors in npredictors_grid) {

  # Inner loop now iterates over your specific run sequence
  for (run in run_start:run_end) {

    tag <- sprintf("n%d_p%d_run%d", nobs, npredictors, run)
    message(sprintf("[%s] %s starting...", format(Sys.time(), "%H:%M:%S"), tag))

    result <- tryCatch(
      run_once(nobs, npredictors, beta1_active, beta2_active, zero_gap_target),
      error = function(e) {
        message(sprintf("  !! %s FAILED: %s", tag, conditionMessage(e)))
        list(meta = list(nobs = nobs, npredictors = npredictors, run = run,
                         timestamp = Sys.time(), error = conditionMessage(e)))
      }
    )

    out_file <- file.path(out_dir, paste0(tag, ".Rdata"))
    save(result, file = out_file, compress = "xz")
    message(sprintf("  -> saved %s", out_file))
  }
}


message("Batch complete.")
