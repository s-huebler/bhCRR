# sim/timing/find_s0.r
#
# Feasibility check: find the smallest s0 in a candidate grid that yields a
# non-null bhCRR fit (nnz >= 3, conv = TRUE) at p = 500, then verify at
# full p = 24619.  s1 is fixed at 0.3; init_method = "LASSO_bic"; all other
# fit_ssl_psdh defaults.
#
# Run from the repo root (locally):
#   Rscript sim/timing/find_s0.r
#
# Writes sim/timing/results/find_s0.csv.

source(file.path(Sys.getenv("REPO_ROOT", getwd()), "sim", "R", "load.r"))
load_bhcrr(Sys.getenv("BHCRR_LOAD", "dev"))
root <- repo_root()

# ---------------------------------------------------------------------------
# Load input
# ---------------------------------------------------------------------------
input_path <- file.path(root, "sim", "timing", "input", "clintall_timing_input.rds")
if (!file.exists(input_path))
  stop("Input not found: ", input_path, "\nRun prep_clintall_timing.r first.")

cat("[find_s0] Loading input...\n")
inp    <- readRDS(input_path)
x_full <- inp$x
y      <- inp$y
n      <- nrow(x_full)
pmax   <- ncol(x_full)
cat(sprintf("[find_s0] n=%d  pmax=%d\n", n, pmax))

S1 <- 0.3

# ---------------------------------------------------------------------------
# Helper: fit at one (s0, p) and return a summary row
# ---------------------------------------------------------------------------
fit_one <- function(s0, p) {
  x_p <- x_full[, seq_len(p)]
  cat(sprintf("[find_s0] Fitting s0=%.3g  p=%d ...", s0, p))
  t0  <- proc.time()[["elapsed"]]
  err <- NA_character_
  fit <- tryCatch(
    fit_ssl_psdh(x_p, y,
                 ss               = c(s0, S1),
                 initial_sparsity = 0.05,
                 theta_b          = ncol(x_p),
                 maxit            = 50L,
                 epsilon          = 1e-4,
                 inner_maxit_start = 1000L,
                 init_method      = "LASSO_bic"),
    error = function(e) { err <<- conditionMessage(e); NULL }
  )
  t_total <- proc.time()[["elapsed"]] - t0

  if (is.null(fit)) {
    cat(sprintf(" ERROR: %s\n", err))
    return(data.frame(s0 = s0, p = p, nnz = NA_integer_,
                      iterations = NA_integer_, conv = NA,
                      max_pip = NA_real_, t_total = t_total,
                      err = err, stringsAsFactors = FALSE))
  }

  nnz     <- sum(fit$coefficients$Estimate != 0)
  iters   <- fit$iterations
  conv    <- fit$conv
  max_pip <- max(fit$pips)
  cat(sprintf(" nnz=%d  iter=%d  conv=%s  max_pip=%.4f  t=%.1fs\n",
              nnz, iters, conv, max_pip, t_total))
  data.frame(s0 = s0, p = p, nnz = nnz, iterations = iters,
             conv = conv, max_pip = max_pip, t_total = t_total,
             err = NA_character_, stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Step 1: Sweep s0 at p = 500
# ---------------------------------------------------------------------------
s0_grid <- c(0.01, 0.02, 0.03, 0.05, 0.08)
cat("\n--- Step 1: sweep s0 at p = 500 ---\n")
rows_500 <- lapply(s0_grid, fit_one, p = 500L)
tbl_500  <- do.call(rbind, rows_500)

# ---------------------------------------------------------------------------
# Step 2: Choose smallest s0 with nnz >= 3 and conv = TRUE
# ---------------------------------------------------------------------------
candidates <- tbl_500[!is.na(tbl_500$nnz) &
                      tbl_500$nnz >= 3L    &
                      !is.na(tbl_500$conv) &
                      tbl_500$conv == TRUE, ]

if (nrow(candidates) == 0L) {
  cat("\n*** No s0 in c(", paste(s0_grid, collapse=", "), ") gave nnz>=3 & conv=TRUE.\n")
  cat("    Extending grid to c(0.12, 0.2) and retrying...\n\n")
  extra <- c(0.12, 0.2)
  extra_rows <- lapply(extra, fit_one, p = 500L)
  tbl_500    <- rbind(tbl_500, do.call(rbind, extra_rows))
  candidates <- tbl_500[!is.na(tbl_500$nnz) &
                        tbl_500$nnz >= 3L    &
                        !is.na(tbl_500$conv) &
                        tbl_500$conv == TRUE, ]
  if (nrow(candidates) == 0L)
    stop("No s0 up to 0.2 gives nnz >= 3 and conv = TRUE at p = 500. ",
         "Manual investigation required.")
}

chosen_s0 <- min(candidates$s0)
cat(sprintf("\nChosen s0 = %.3g  (smallest with nnz>=3 & conv=TRUE at p=500)\n", chosen_s0))

# ---------------------------------------------------------------------------
# Step 3: Fit chosen s0 at full p
# ---------------------------------------------------------------------------
cat(sprintf("\n--- Step 3: fit chosen s0=%.3g at p=%d ---\n", chosen_s0, pmax))
row_full <- fit_one(chosen_s0, pmax)

# If nnz at full p is suspiciously large, try one step smaller
if (!is.na(row_full$nnz) && row_full$nnz > 120L) {
  smaller_candidates <- candidates[candidates$s0 < chosen_s0, ]
  if (nrow(smaller_candidates) > 0L) {
    alt_s0 <- max(smaller_candidates$s0)
    cat(sprintf("\nnnz=%d > 120 at full p; trying one step smaller: s0=%.3g\n",
                row_full$nnz, alt_s0))
    row_alt <- fit_one(alt_s0, pmax)
    row_full <- rbind(row_full, row_alt)
  } else {
    cat(sprintf("\nnnz=%d > 120 but no smaller s0 available.\n", row_full$nnz))
  }
}

# ---------------------------------------------------------------------------
# Step 4: Combine and print full table
# ---------------------------------------------------------------------------
all_rows <- rbind(tbl_500, row_full)
# Mark the chosen s0
all_rows$chosen <- all_rows$s0 == chosen_s0

cat("\n=== Full results table ===\n")
print(all_rows[, c("s0","p","nnz","iterations","conv","max_pip","t_total","chosen")],
      row.names = FALSE, digits = 4)

# ---------------------------------------------------------------------------
# Step 5: Write CSV
# ---------------------------------------------------------------------------
out_dir <- file.path(root, "sim", "timing", "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_csv <- file.path(out_dir, "find_s0.csv")
write.csv(all_rows, out_csv, row.names = FALSE)
cat(sprintf("\n[find_s0] Written to: %s\n", out_csv))
cat(sprintf("[find_s0] Recommended: ss=c(%.3g, %.3g)\n", chosen_s0, S1))
