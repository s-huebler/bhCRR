#################### run_timing.R ####################
# Loads ONE timing config's .rds (written by prep_timing_data.R), fits ONE
# bhCRR model, and writes a self-contained result .rds to TIMING_OUT_DIR.
# No autotuning, no comparator models, no data generation.
#
# Environment variables (all required except those with defaults):
#   TIMING_TAG           which config to run, e.g. "n200-p5000"
#   TIMING_DATA_DIR      directory containing timing_<tag>.rds
#   TIMING_OUT_DIR       directory for result files
#   REPO_ROOT            path to the bhCRR repo checkout
#   SS0, SS1             spike-and-slab scales (defaults 0.001, 0.4)
#   THETA_A, THETA_B     EM hyperparameters (defaults 1, 1)
#   INITIAL_SPARSITY     default 0.05
#   MAXIT                default 50
#   EPSILON              default 1e-04
#   INIT_MODE            how fit_ssl_psdh is initialized; default "bic"
#                          "bic"  fit a LASSO path with fastCrrp over
#                                 init_lam_path and pass the coefficients from
#                                 the lowest-BIC lambda as init
#                          "cv"   init=NULL, so fit_ssl_psdh runs its own
#                                 5-fold cv_fastCrrp_cpp search internally
#                          "zero" init=rep(0,p), skipping any initial fit
#   INNER_MAXIT_START    default 1000


#################### Config from environment ####################

get_env <- function(name, default = NULL) {
  v <- Sys.getenv(name, unset = NA_character_)
  if (is.na(v) || v == "") {
    if (is.null(default)) stop("Required environment variable not set: ", name)
    return(default)
  }
  v
}

timing_tag      <- get_env("TIMING_TAG")
timing_data_dir <- get_env("TIMING_DATA_DIR")
timing_out_dir  <- get_env("TIMING_OUT_DIR")
repo_root       <- get_env("REPO_ROOT")
ss0             <- as.numeric(get_env("SS0",              "0.001"))
ss1             <- as.numeric(get_env("SS1",              "0.4"))
theta_a         <- as.numeric(get_env("THETA_A",          "1"))
theta_b         <- as.numeric(get_env("THETA_B",          "1"))
initial_sparsity <- as.numeric(get_env("INITIAL_SPARSITY", "0.05"))
maxit           <- as.integer(get_env("MAXIT",            "50"))
epsilon         <- as.numeric(get_env("EPSILON",          "1e-04"))
init_mode       <- get_env("INIT_MODE",                   "bic")
inner_maxit_start <- as.integer(get_env("INNER_MAXIT_START", "1000"))

if (!init_mode %in% c("bic", "cv", "zero"))
  stop("INIT_MODE must be 'bic', 'cv' or 'zero', got: '", init_mode, "'")

dir.create(timing_out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(timing_out_dir, paste0("timing_", timing_tag, ".rds"))

message("run_timing.R")
message("  TIMING_TAG      = ", timing_tag)
message("  TIMING_DATA_DIR = ", timing_data_dir)
message("  TIMING_OUT_DIR  = ", timing_out_dir)
message("  REPO_ROOT       = ", repo_root)
message(sprintf("  ss=(%s,%s)  theta_a=%s  theta_b=%s  maxit=%d  epsilon=%s",
                ss0, ss1, theta_a, theta_b, maxit, epsilon))
message("  init_mode       = ", init_mode)
message("  inner_maxit_start = ", inner_maxit_start)


#################### Load R sources ####################

message("\nLoading R sources (includes Rcpp compilation)...")
t_src_start <- Sys.time()

suppressPackageStartupMessages({
  library(fastcmprsk)
  library(survival)
})

src <- function(...) source(file.path(repo_root, ...), local = FALSE)
cpp <- function(f)   Rcpp::sourceCpp(file.path(repo_root, f))

src("R/helpers.r")
src("R/utils.r")
src("R/dedupe_warnings.r")
src("R/update_betas.r")
src("R/expected_inclusion_probs.r")
src("R/expected_penalty_weights.r")
src("R/update_mixture_prob.r")
cpp("src/cv_fastcrrp.cpp")
cpp("src/RcppExports.cpp")
src("R/cv_fastCrrp_cpp.R")
src("R/fit_ssl_psdh.r")
src("R/generate_foldid.r")
src("R/predict_from_ssl_psdh.r")
src("R/wolbers_c.r")
src("R/threshold.R")

t_src_wall <- as.numeric(Sys.time() - t_src_start, units = "secs")
message(sprintf("  sources loaded + compiled in %.1f s", t_src_wall))


#################### Load timing data ####################

message("\nLoading timing data...")
t_data_start <- Sys.time()

# Helper: write a minimal failure result then rethrow. Used for pre-fit errors
# (data not found, format wrong) so the job is traceable even if it never runs.
.write_preflight_failure <- function(msg) {
  message("ERROR: ", msg)
  tryCatch(
    saveRDS(list(
      tag         = timing_tag,
      n           = NA_integer_, p = NA_integer_,
      ss          = c(ss0, ss1), theta_a = theta_a, theta_b = theta_b,
      maxit       = maxit, epsilon = epsilon, init_mode = init_mode,
      timings     = NULL, memory = NULL,
      converged   = NA, iterations = NA_integer_, n_nonzero = NA_integer_,
      coefficients = NULL, event_mix = NULL,
      error       = msg,
      hostname    = Sys.info()[["nodename"]],
      slurm_job_id = Sys.getenv("SLURM_JOB_ID",        unset = NA_character_),
      slurm_array_task_id = Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA_character_),
      timestamp   = Sys.time(),
      r_version   = R.version$version.string
    ), out_file),
    error = function(e)
      message("WARNING: could not write failure result to ", out_file, ": ", conditionMessage(e))
  )
  message("Failure result written to: ", out_file)
  stop(msg, call. = FALSE)
}

data_file <- file.path(timing_data_dir, paste0("timing_", timing_tag, ".rds"))
if (!file.exists(data_file))
  .write_preflight_failure(paste0(
    "Timing data file not found: ", data_file,
    "\nRun prep_timing_data.R first and verify TIMING_DATA_DIR."
  ))

data_obj <- tryCatch(
  readRDS(data_file),
  error = function(e) .write_preflight_failure(paste0(
    "readRDS() failed on ", data_file, ".\n",
    "The file may have been written with save() instead of saveRDS().\n",
    "Re-run prep_timing_data.R and retry.\n",
    "Original error: ", conditionMessage(e)
  ))
)

x <- data_obj$x
y <- data_obj$y

t_data_wall <- as.numeric(Sys.time() - t_data_start, units = "secs")

n <- nrow(x)
p <- ncol(x)
message(sprintf("  loaded: n=%d  p=%d  in %.1f s", n, p, t_data_wall))

stopifnot(
  is.matrix(x), is.numeric(x),
  is.matrix(y), ncol(y) == 2L,
  identical(rownames(x), rownames(y)),
  nrow(x) == data_obj$meta$n,
  ncol(x) == data_obj$meta$p
)

event_mix <- table(y[, "status_relapse"])
message("  event mix: ", paste(names(event_mix), as.integer(event_mix), sep = "=", collapse = "  "))

if (sum(y[, "status_relapse"] == 1L) == 0L)
  stop("No cause-1 events in this config -- cannot fit.")


#################### Initialization ####################

init_lam_path <- 10^seq(log10(0.1), log10(0.001), length.out = 10L)

init_arg      <- NULL
init_cv_wall  <- NULL     # set when init_mode == "cv"
init_bic_wall <- NULL     # set when init_mode == "bic"
init_bic_sel  <- NULL     # provenance for the BIC-selected initialization

if (init_mode == "zero") {

  init_arg <- rep(0, p)
  message("\nInitialization: zero vector (no initial model is fit).")

} else if (init_mode == "cv") {

  # Time a standalone call to cv_fastCrrp_cpp() with the same arguments
  # fit_ssl_psdh uses internally, to estimate initialization cost separately
  # from the EM loop.
  # NOTE: this is a SEPARATE call. With init = NULL, fit_ssl_psdh runs
  # cv_fastCrrp_cpp AGAIN internally, so this time approximates the init cost
  # but is NOT subtracted from t_fit_wall below.
  message("\nTiming standalone cv_fastCrrp_cpp() (init cost estimate)...")
  t_init_start <- Sys.time()
  init_est <- cv_fastCrrp_cpp(x, y[, 1L], y[, 2L], k = 5L,
                              penalty = "LASSO",
                              lambda_path = init_lam_path,
                              tuning = "wolbers",
                              eval_quantile = 0.5)
  init_cv_wall <- as.numeric(Sys.time() - t_init_start, units = "secs")
  message(sprintf("  estimated init cost: %.1f s  (separate call -- fit will repeat this internally)",
                  init_cv_wall))
  rm(init_est)
  invisible(gc())
  init_arg <- NULL          # fit_ssl_psdh performs its own CV search

} else if (init_mode == "bic") {

  # Fit ONE LASSO path with fastCrrp over init_lam_path and take the
  # coefficients at the lowest-BIC lambda as the starting values.
  #
  # Why this exists: the "cv" path splits into 5 folds, and every n200-* config
  # here carries only 19 cause-1 events, leaving under 4 events per held-out
  # fold. BIC selection needs no fold splitting, so it does not degrade with
  # few events, and it fits the initial model exactly once instead of once per
  # fold plus once more inside fit_ssl_psdh.
  #
  # x is already column-standardized by prep_timing_data.R, so these
  # coefficients are on the scale fit_ssl_psdh operates on -- no
  # back-transformation is needed.
  #
  # The matrix form of the formula (~ x) is deliberate: a formula with p named
  # terms would be unusable at p = 24618.
  message("\nFitting LASSO initialization path, selecting by BIC...")
  t_init_start <- Sys.time()

  init_fit <- fastcmprsk::fastCrrp(
    fastcmprsk::Crisk(y[, 1L], y[, 2L], cencode = 0, failcode = 1) ~ x,
    penalty = "LASSO",
    lambda  = init_lam_path
  )

  beta_path <- as.matrix(init_fit$coef)
  if (nrow(beta_path) != p && ncol(beta_path) == p) beta_path <- t(beta_path)
  if (nrow(beta_path) != p)
    stop("fastCrrp() returned a coefficient matrix of unexpected dimension: ",
         paste(dim(beta_path), collapse = " x "), " (expected ", p, " rows).")

  # BIC over the path; df = number of nonzero coefficients at each lambda.
  df_path <- colSums(beta_path != 0)
  loglik  <- as.numeric(init_fit$logLik)
  if (length(loglik) != ncol(beta_path))
    stop("fastCrrp() logLik has length ", length(loglik),
         " but the coefficient path has ", ncol(beta_path), " columns.")
  bic_path <- -2 * loglik + log(n) * df_path

  # Consider only finite BIC values, and where fastCrrp reports per-lambda
  # convergence, only converged solutions.
  valid <- is.finite(bic_path)
  conv  <- init_fit$converged
  if (!is.null(conv) && length(conv) == length(bic_path))
    valid <- valid & as.logical(conv)
  if (!any(valid))
    stop("No finite, converged LASSO solution on the initialization path; ",
         "cannot select an initial estimate by BIC.")

  bic_for_selection <- bic_path
  bic_for_selection[!valid] <- Inf
  best_idx <- which.min(bic_for_selection)

  init_arg <- as.numeric(beta_path[, best_idx])
  if (length(init_arg) != p || any(!is.finite(init_arg)))
    stop("BIC-selected initial estimate is not a finite length-", p, " vector.")

  # fastCrrp may reorder the supplied path, so record the lambda it actually used.
  lam_used <- init_fit$lambda
  if (is.null(lam_used)) lam_used <- init_fit$lambda.path

  init_bic_wall <- as.numeric(Sys.time() - t_init_start, units = "secs")

  init_bic_sel <- list(
    index       = as.integer(best_idx),
    lambda      = if (is.null(lam_used)) NA_real_ else as.numeric(lam_used[best_idx]),
    bic         = as.numeric(bic_path[best_idx]),
    df          = as.integer(df_path[best_idx]),
    n_nonzero   = as.integer(sum(init_arg != 0)),
    lambda_path = if (is.null(lam_used)) NA_real_ else as.numeric(lam_used),
    bic_path    = as.numeric(bic_path),
    df_path     = as.integer(df_path),
    converged   = if (is.null(conv)) NA else as.logical(conv[best_idx])
  )

  message(sprintf(
    "  BIC-selected: index=%d  lambda=%.5g  bic=%.2f  df=%d  nnz=%d  (%.1f s)",
    best_idx, init_bic_sel$lambda, init_bic_sel$bic,
    init_bic_sel$df, init_bic_sel$n_nonzero, init_bic_wall))

  rm(init_fit, beta_path)
  invisible(gc())
}


#################### Fit ####################

message("\nFitting fit_ssl_psdh()...")
message(sprintf("  n=%d  p=%d  ss=(%s,%s)  theta_b=%s  init=%s  maxit=%d",
                n, p, ss0, ss1, theta_b, init_mode, maxit))

invisible(gc(reset = TRUE))   # reset peak-memory counters immediately before fit
pt_fit_start  <- proc.time()
t_fit_start   <- Sys.time()

error_msg <- NULL
mod <- tryCatch(
  fit_ssl_psdh(x, y,
               ss               = c(ss0, ss1),
               initial_sparsity = initial_sparsity,
               theta_a          = theta_a,
               theta_b          = theta_b,
               maxit            = maxit,
               epsilon          = epsilon,
               init             = init_arg,
               init_lam_path    = init_lam_path,
               inner_maxit_start = inner_maxit_start),
  error = function(e) {
    message("  !! fit_ssl_psdh() FAILED: ", conditionMessage(e))
    e
  }
)

pt_fit   <- proc.time() - pt_fit_start
t_fit_wall <- as.numeric(Sys.time() - t_fit_start, units = "secs")
gc_after <- gc()

if (inherits(mod, "error")) {
  error_msg <- conditionMessage(mod)
  mod       <- NULL
}


#################### Extract results ####################

converged  <- if (is.null(mod)) NA else mod$conv
iterations <- if (is.null(mod)) NA_integer_ else mod$iterations
n_nonzero  <- if (is.null(mod)) NA_integer_ else sum(mod$coefficients$Estimate != 0)
coefs      <- if (is.null(mod)) NULL else mod$coefficients  # data.frame Variable/Estimate

# R 4.5 added a "limit (Mb)" column, making gc() return 7 columns rather than 6.
# Avoid hardcoding column indices: find "max used" by name and take the "(Mb)"
# column immediately after it, which is always the "max used (Mb)" value.
.gc_max_used_col <- which(colnames(gc_after) == "max used")
.gc_mb_col       <- .gc_max_used_col + 1L   # "(Mb)" immediately follows "max used"
memory <- list(
  ncells_max_used_mb = gc_after["Ncells", .gc_mb_col],
  vcells_max_used_mb = gc_after["Vcells", .gc_mb_col],
  vcells_max_used_gb = gc_after["Vcells", .gc_mb_col] / 1024
)

timings <- list(
  src_compile_wall_sec  = t_src_wall,
  data_load_wall_sec    = t_data_wall,
  init_cv_wall_sec      = init_cv_wall,   # NULL unless init_mode == "cv"
  init_bic_wall_sec     = init_bic_wall,  # NULL unless init_mode == "bic"
  fit_wall_sec          = t_fit_wall,
  fit_proc_user_sec     = as.numeric(pt_fit["user.self"]),
  fit_proc_sys_sec      = as.numeric(pt_fit["sys.self"]),
  fit_proc_elapsed_sec  = as.numeric(pt_fit["elapsed"])
)


#################### One-line summary (greppable from SLURM log) ####################

message(sprintf(
  "TIMING | tag=%-12s | n=%4d p=%6d | init=%4s | fit=%.1fs | iter=%s conv=%s | nnz=%s | err=%s",
  timing_tag, n, p, init_mode,
  t_fit_wall,
  if (is.na(iterations)) "NA" else as.character(iterations),
  if (is.na(converged))  "NA" else as.character(converged),
  if (is.na(n_nonzero))  "NA" else as.character(n_nonzero),
  if (is.null(error_msg)) "none" else error_msg
))

if (!is.null(init_cv_wall))
  message(sprintf("TIMING | tag=%-12s | init_cv_est=%.1fs  (separate call; not subtracted from fit time)",
                  timing_tag, init_cv_wall))

if (!is.null(init_bic_wall))
  message(sprintf(
    "TIMING | tag=%-12s | init_bic=%.1fs lambda=%.5g df=%d nnz=%d  (fit BEFORE the model; not included in fit time)",
    timing_tag, init_bic_wall, init_bic_sel$lambda,
    init_bic_sel$df, init_bic_sel$n_nonzero))

message(sprintf("TIMING | tag=%-12s | peak_mem_vcells=%.1f MB (%.2f GB)",
                timing_tag, memory$vcells_max_used_mb, memory$vcells_max_used_gb))


#################### Save result ####################

result <- list(
  tag                = timing_tag,
  n                  = n,
  p                  = p,
  ss                 = c(ss0, ss1),
  theta_a            = theta_a,
  theta_b            = theta_b,
  maxit              = maxit,
  epsilon            = epsilon,
  init_mode          = init_mode,
  init_bic_selection = init_bic_sel,   # NULL unless init_mode == "bic"
  timings            = timings,
  memory             = memory,
  converged          = converged,
  iterations         = iterations,
  n_nonzero          = n_nonzero,
  coefficients       = coefs,
  event_mix          = event_mix,
  error              = error_msg,
  hostname           = Sys.info()[["nodename"]],
  slurm_job_id       = Sys.getenv("SLURM_JOB_ID",         unset = NA_character_),
  slurm_array_task_id = Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA_character_),
  timestamp          = Sys.time(),
  r_version          = R.version$version.string
)

saveRDS(result, out_file)
message("Result written to: ", out_file)
