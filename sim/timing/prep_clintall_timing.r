# sim/timing/prep_clintall_timing.r
#
# LOCAL ONLY.  Never run on CHPC (ClinTALL .rds are gitignored).
#
# Reads the two ClinTALL data files, aligns by sample ID, standardizes, draws
# the fixed column permutation, builds the p-grid, and writes:
#   sim/timing/input/clintall_timing_input.rds   (gitignored)
#   sim/timing/grid.tsv                          (committed)
#   sim/timing/input_manifest.txt                (committed)
#
# Run from the repo root:
#   Rscript sim/timing/prep_clintall_timing.r
#   Rscript sim/timing/prep_clintall_timing.r grid_step=2000   # coarser grid
#
# To regenerate ONLY grid.tsv without touching input/ (e.g. to coarsen the grid):
#   Rscript -e "
#     p <- sort(unique(c(seq(500L, 24619L, by = 2000L), 24619L)))
#     d <- data.frame(p_index = seq_along(p), p = p)
#     write.table(d, 'sim/timing/grid.tsv', sep='\t', quote=FALSE, row.names=FALSE)
#     cat('G=', length(p), '\n')
#   "
# Do NOT source this file for that purpose — sourcing runs the full script and
# overwrites input/clintall_timing_input.rds.

# ---------------------------------------------------------------------------
# Exported helpers (usable via source() or Rscript -e)
# ---------------------------------------------------------------------------

#' Build the p-grid for the timing study.
#'
#' Returns unique(c(seq(start, pmax, by = step), pmax)), sorted.
#' Default start=500, step=2000 gives the production grid.
make_timing_grid <- function(pmax, start = 500L, step = 2000L) {
  sort(unique(c(seq(start, pmax, by = step), pmax)))
}

#' Write a p-grid vector to grid.tsv.
write_grid <- function(p_grid, path) {
  grid_df <- data.frame(p_index = seq_along(p_grid), p = p_grid)
  write.table(grid_df, path, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("[prep] Wrote grid.tsv (%d rows): p = %d, ..., %d\n",
              length(p_grid), p_grid[1], p_grid[length(p_grid)]))
  invisible(grid_df)
}

# Parse optional grid_step=<n> from command-line args
.args      <- commandArgs(trailingOnly = TRUE)
.step_arg  <- grep("^grid_step=", .args, value = TRUE)
GRID_STEP  <- if (length(.step_arg)) as.integer(sub("^grid_step=", "", .step_arg[1])) else 2000L

source(file.path(Sys.getenv("REPO_ROOT", getwd()), "sim", "R", "load.r"))
load_bhcrr(Sys.getenv("BHCRR_LOAD", "dev"))

root <- repo_root()

# ---------------------------------------------------------------------------
# 1. Read ClinTALL data
# ---------------------------------------------------------------------------
cat("[prep] Reading ClinTALL data files...\n")
gexp <- readRDS(file.path(root, "inst", "ClinTall", "dat_design_gexp.rds"))
cr   <- readRDS(file.path(root, "inst", "ClinTall", "dat_cr.rds"))

# Align by sample ID; stop on any mismatch
shared <- intersect(rownames(gexp), cr$sample_id)
missing_from_cr   <- setdiff(rownames(gexp), cr$sample_id)
missing_from_gexp <- setdiff(cr$sample_id,   rownames(gexp))
if (length(missing_from_cr) > 0L)
  stop("IDs in gexp rownames not found in dat_cr$sample_id: ",
       paste(head(missing_from_cr, 5), collapse = ", "),
       if (length(missing_from_cr) > 5) " ..." else "")
if (length(missing_from_gexp) > 0L)
  stop("IDs in dat_cr$sample_id not found in gexp rownames: ",
       paste(head(missing_from_gexp, 5), collapse = ", "),
       if (length(missing_from_gexp) > 5) " ..." else "")

n_full <- length(shared)
cat(sprintf("[prep] Aligned %d samples\n", n_full))

x_raw <- gexp[shared, ]
cr_aln <- cr[match(shared, cr$sample_id), ]

# Coerce to numeric matrix
x_raw <- as.matrix(x_raw)
storage.mode(x_raw) <- "double"

# ---------------------------------------------------------------------------
# 2. Build outcome matrix
# ---------------------------------------------------------------------------
y <- cbind(tte = cr_aln$tte, status_relapse = cr_aln$status_relapse)
rownames(y) <- shared

n_cause1    <- sum(y[, "status_relapse"] == 1L)
n_cause2    <- sum(y[, "status_relapse"] == 2L)
n_censored  <- sum(y[, "status_relapse"] == 0L)

cat(sprintf("[prep] Event mix: cause-1=%d  cause-2=%d  censored=%d\n",
            n_cause1, n_cause2, n_censored))

if (n_cause1   != 120L) stop("Expected 120 cause-1 events, got ",   n_cause1)
if (n_cause2   !=  97L) stop("Expected 97 cause-2 events, got ",    n_cause2)
if (n_censored != 1092L) stop("Expected 1092 censored, got ",       n_censored)
cat("[prep] Event-count assertions passed.\n")

# ---------------------------------------------------------------------------
# 3. Drop NA and zero-variance genes; standardize
# ---------------------------------------------------------------------------
cat("[prep] Dropping NA / zero-variance genes and standardizing...\n")

any_na  <- apply(x_raw, 2L, anyNA)
col_var <- apply(x_raw, 2L, var, na.rm = TRUE)
keep    <- !any_na & col_var > 0
n_dropped <- sum(!keep)
cat(sprintf("[prep] Dropped %d genes (%d NA, %d zero-variance)\n",
            n_dropped,
            sum(any_na),
            sum(!any_na & col_var <= 0)))

x <- scale(x_raw[, keep])   # center and scale on full n
pmax_val <- ncol(x)
cat(sprintf("[prep] After filtering: n=%d  pmax=%d\n", nrow(x), pmax_val))

# ---------------------------------------------------------------------------
# 4. Fixed column permutation (nested subsets)
# ---------------------------------------------------------------------------
SEED <- 20261030L
set.seed(SEED)
gene_order <- sample(pmax_val)
x <- x[, gene_order]   # reorder so first p columns = nested subset

# ---------------------------------------------------------------------------
# 5. Build p-grid and write grid.tsv
# ---------------------------------------------------------------------------
p_grid <- make_timing_grid(pmax_val, start = 500L, step = GRID_STEP)
G      <- length(p_grid)
cat(sprintf("[prep] Grid: step=%d  G=%d  p values from %d to %d\n",
            GRID_STEP, G, p_grid[1], p_grid[G]))

grid_path <- file.path(root, "sim", "timing", "grid.tsv")
write_grid(p_grid, grid_path)

# ---------------------------------------------------------------------------
# 6. Save input .rds
# ---------------------------------------------------------------------------
input_dir  <- file.path(root, "sim", "timing", "input")
dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)
input_path <- file.path(input_dir, "clintall_timing_input.rds")

saveRDS(
  list(
    x          = x,
    y          = y,
    gene_order = gene_order,
    seed       = SEED,
    n_dropped  = n_dropped,
    created    = Sys.time()
  ),
  input_path,
  compress = "xz"
)
cat("[prep] Saved input RDS to:", input_path, "\n")

# ---------------------------------------------------------------------------
# 7. Write input_manifest.txt (committed)
# ---------------------------------------------------------------------------
md5_val <- tools::md5sum(input_path)
git_hash <- tryCatch(
  trimws(system2("git", c("-C", root, "rev-parse", "--short", "HEAD"),
                 stdout = TRUE, stderr = FALSE)),
  error = function(e) "unknown"
)

manifest_lines <- c(
  paste0("md5=",          md5_val[[1]]),
  paste0("dim_x=",        nrow(x), "x", ncol(x)),
  paste0("n_cause1=",     n_cause1),
  paste0("n_cause2=",     n_cause2),
  paste0("n_censored=",   n_censored),
  paste0("n_dropped=",    n_dropped),
  paste0("pmax=",         pmax_val),
  paste0("G=",            G),
  paste0("seed=",         SEED),
  paste0("git_hash=",     git_hash),
  paste0("created=",      format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
)

manifest_path <- file.path(root, "sim", "timing", "input_manifest.txt")
writeLines(manifest_lines, manifest_path)
cat("[prep] Wrote manifest to:", manifest_path, "\n")
cat("\nManifest contents:\n")
cat(paste(manifest_lines, collapse = "\n"), "\n")
cat("\n[prep] Done.\n")
