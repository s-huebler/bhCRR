# sim/timing/fit_one_timing.r
#
# One SLURM array task = one fit_ssl_psdh call at one (p, rep).
#
# Environment variables:
#   TIMING_INPUT   path to clintall_timing_input.rds  (required)
#   TIMING_OUT_DIR directory for result files          (required)
#   TIMING_TAG     short label for this batch          (required)
#   REPS           number of replicates per p value    (default 3)
#   INIT_METHOD    "LASSO_bic" or "zero"               (default "LASSO_bic")
#   SS0            spike scale                         (default 0.002)
#   SS1            slab scale                          (default 0.3)
#   SLURM_ARRAY_TASK_ID  0-based task index (0 = p_grid[1], rep 1)
#
# Array geometry:
#   G   = number of rows in grid.tsv
#   i   = SLURM_ARRAY_TASK_ID (0-based)
#   p   = grid$p[ (i %% G) + 1 ]       indices 0..G-1 are rep 1
#   rep = i %/% G + 1

# ---------------------------------------------------------------------------
# 0. Load bhCRR
# ---------------------------------------------------------------------------
source(file.path(Sys.getenv("REPO_ROOT", getwd()), "sim", "R", "load.r"))
load_bhcrr(Sys.getenv("BHCRR_LOAD", "dev"))
root <- repo_root()

# ---------------------------------------------------------------------------
# 1. Read environment variables
# ---------------------------------------------------------------------------
.req_env <- function(name) {
  v <- Sys.getenv(name, unset = "")
  if (!nzchar(v)) stop("Required environment variable not set: ", name)
  v
}
.opt_env <- function(name, default) {
  v <- Sys.getenv(name, unset = "")
  if (nzchar(v)) v else default
}

timing_input  <- .req_env("TIMING_INPUT")
timing_out_dir <- .req_env("TIMING_OUT_DIR")
timing_tag    <- .req_env("TIMING_TAG")
reps          <- as.integer(.opt_env("REPS",        "3"))
init_method   <- .opt_env("INIT_METHOD", "LASSO_bic")
SS0           <- as.numeric(.opt_env("SS0",          "0.002"))
SS1           <- as.numeric(.opt_env("SS1",          "0.3"))

if (!init_method %in% c("LASSO_bic", "zero"))
  stop("INIT_METHOD must be 'LASSO_bic' or 'zero', got: '", init_method, "'")

task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "0"))

# ---------------------------------------------------------------------------
# 2. Read the p-grid
# ---------------------------------------------------------------------------
grid_path <- file.path(root, "sim", "timing", "grid.tsv")
if (!file.exists(grid_path))
  stop("grid.tsv not found: ", grid_path,
       "\nRun prep_clintall_timing.r first.")
grid <- read.table(grid_path, header = TRUE, sep = "\t")
G    <- nrow(grid)

p   <- grid$p[(task_id %% G) + 1L]
rep <- task_id %/% G + 1L

cat(sprintf("[timing] tag=%s  task_id=%d  G=%d  p=%d  rep=%d\n",
            timing_tag, task_id, G, p, rep))
cat(sprintf("[timing] init_method=%s  ss=(%g,%g)  reps=%d\n",
            init_method, SS0, SS1, reps))

dir.create(timing_out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(timing_out_dir,
                      sprintf("res_p%d_rep%d.rds", p, rep))

# ---------------------------------------------------------------------------
# 3. Read input data (timed separately)
# ---------------------------------------------------------------------------
cat(sprintf("[timing] Reading input from: %s\n", timing_input))
t_read_start <- proc.time()[["elapsed"]]
inp  <- readRDS(timing_input)
t_read <- proc.time()[["elapsed"]] - t_read_start
cat(sprintf("[timing] Read in %.2f s\n", t_read))

x_full   <- inp$x
y        <- inp$y
n        <- nrow(x_full)
n_events1 <- sum(y[, "status_relapse"] == 1L)

cat(sprintf("[timing] n=%d  pmax=%d  n_events1=%d\n", n, ncol(x_full), n_events1))

# ---------------------------------------------------------------------------
# 4. Untimed warm-up at p = 50, then gc
# ---------------------------------------------------------------------------
cat("[timing] Warm-up fit at p=50...\n")
invisible(tryCatch(
  fit_ssl_psdh(x_full[, 1:50], y,
               ss               = c(SS0, SS1),
               initial_sparsity = 0.05,
               theta_b          = 50L,
               maxit            = 5L,
               epsilon          = 1e-4,
               init_method      = "zero"),
  error = function(e) NULL
))
invisible(gc())

# ---------------------------------------------------------------------------
# 5. Timed fit
# ---------------------------------------------------------------------------
x_p <- x_full[, seq_len(p)]

cat(sprintf("[timing] Fitting p=%d with init_method=%s...\n", p, init_method))

err_msg <- NA_character_
fit     <- NULL

t_total_start <- proc.time()[["elapsed"]]

fit <- tryCatch(
  fit_ssl_psdh(x_p, y,
               ss                = c(SS0, SS1),
               initial_sparsity  = 0.05,
               theta_b           = ncol(x_p),
               maxit             = 50L,
               epsilon           = 1e-4,
               inner_maxit_start = 1000L,
               init_method       = init_method),
  error = function(e) {
    cat(sprintf("[timing] ERROR: %s\n", conditionMessage(e)))
    err_msg <<- conditionMessage(e)
    NULL
  }
)

t_total <- proc.time()[["elapsed"]] - t_total_start

# ---------------------------------------------------------------------------
# 6. Extract timing and metadata
# ---------------------------------------------------------------------------
if (!is.null(fit)) {
  t_init      <- if (is.null(fit$init_meta)) 0 else
                   fit$init_meta$elapsed_sec %||% 0
  t_em        <- t_total - t_init
  iterations  <- fit$iterations
  conv        <- fit$conv
  nnz         <- sum(fit$coefficients$Estimate != 0)
  init_lambda <- fit$init_meta$lambda     %||% NA_real_
  init_df     <- fit$init_meta$df         %||% NA_integer_
} else {
  t_init      <- NA_real_
  t_em        <- NA_real_
  iterations  <- NA_integer_
  conv        <- NA
  nnz         <- NA_integer_
  init_lambda <- NA_real_
  init_df     <- NA_integer_
}

# Peak RSS (Linux only)
peak_rss_mb <- tryCatch({
  lines  <- readLines("/proc/self/status")
  hwm    <- grep("^VmHWM:", lines, value = TRUE)
  if (length(hwm) == 1L)
    as.numeric(gsub("[^0-9]", "", hwm)) / 1024
  else
    NA_real_
}, error = function(e) NA_real_)

# CPU model (Linux only)
cpu_model <- tryCatch({
  lines <- readLines("/proc/cpuinfo")
  mn    <- grep("^model name", lines, value = TRUE)
  if (length(mn) >= 1L) trimws(sub(".*:\\s*", "", mn[1L])) else NA_character_
}, error = function(e) NA_character_)

bhcrr_hash <- tryCatch(
  trimws(system2("git", c("-C", root, "rev-parse", "--short", "HEAD"),
                 stdout = TRUE, stderr = FALSE)[[1]]),
  error = function(e) NA_character_
)

fastcmprsk_ver <- tryCatch(
  as.character(packageVersion("fastcmprsk")),
  error = function(e) NA_character_
)

# ---------------------------------------------------------------------------
# 7. One grep-able TIMING line to stdout
# ---------------------------------------------------------------------------
cat(sprintf(
  "TIMING | tag=%-20s | p=%6d rep=%d | init=%-10s | t_total=%.2fs t_init=%.2fs t_em=%.2fs | iter=%s conv=%s nnz=%s | err=%s\n",
  timing_tag, p, rep, init_method,
  t_total,
  if (is.na(t_init)) 0 else t_init,
  if (is.na(t_em))   t_total else t_em,
  if (is.na(iterations)) "NA" else as.character(iterations),
  if (is.na(conv))   "NA" else as.character(conv),
  if (is.na(nnz))    "NA" else as.character(nnz),
  if (is.na(err_msg)) "none" else err_msg
))

# ---------------------------------------------------------------------------
# 8. Write result row
# ---------------------------------------------------------------------------
row <- list(
  tag               = timing_tag,
  task_id           = task_id,
  p                 = p,
  rep               = rep,
  n                 = n,
  n_events1         = n_events1,
  init_method       = init_method,
  s0                = SS0,
  s1                = SS1,
  t_read            = t_read,
  t_init            = t_init,
  t_em              = t_em,
  t_total           = t_total,
  iterations        = iterations,
  conv              = conv,
  nnz               = nnz,
  init_lambda       = init_lambda,
  init_df           = init_df,
  err               = err_msg,
  peak_rss_mb       = peak_rss_mb,
  host              = Sys.info()[["nodename"]],
  cpu_model         = cpu_model,
  slurm_job_id      = Sys.getenv("SLURM_JOB_ID",        unset = NA_character_),
  time_limit        = Sys.getenv("SLURM_JOB_TIMELIMIT",  unset = NA_character_),
  R_version         = R.version$version.string,
  fastcmprsk_version = fastcmprsk_ver,
  bhCRR_git_hash    = bhcrr_hash,
  timestamp         = Sys.time()
)

saveRDS(row, out_file)
cat(sprintf("[timing] Wrote result to: %s\n", out_file))
