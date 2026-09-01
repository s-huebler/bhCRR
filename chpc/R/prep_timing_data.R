#################### prep_timing_data.R ####################
# Materialize nested gene subsets of the ClinTALL data for timing /
# feasibility experiments. Writes one .rds per config plus a manifest CSV
# into TIMING_DATA_DIR. No model fitting here.
#
# LADDER: n is held at the full cohort (1309) and only p is swept.
# An earlier ladder also varied n at 200, but every n=200 config converged to
# the null model regardless of p -- 19 cause-1 events carried no signal for the
# BIC-selected LASSO initialization, so df was 0 at every lambda. Sweeping n
# below the full cohort produced no usable timing information and was dropped.
#
# Environment variables (all have defaults):
#   CLINTALL_DIR     directory containing dat_design.rds / dat_cr.rds
#                    default: inst/ClinTall
#   TIMING_DATA_DIR  output directory for the subsets
#                    default: inst/ClinTall/timing
#   TIMING_SEED      integer seed
#                    default: 20260901


#################### Config from environment ####################

get_env <- function(name, default = NULL) {
  v <- Sys.getenv(name, unset = NA_character_)
  if (is.na(v) || v == "") {
    if (is.null(default)) stop("Required environment variable not set: ", name)
    return(default)
  }
  v
}

clintall_dir    <- get_env("CLINTALL_DIR",    "inst/ClinTall")
timing_data_dir <- get_env("TIMING_DATA_DIR", "inst/ClinTall/timing")
seed            <- as.integer(get_env("TIMING_SEED", "20260901"))

dir.create(timing_data_dir, showWarnings = FALSE, recursive = TRUE)

message("prep_timing_data.R")
message("  CLINTALL_DIR    = ", clintall_dir)
message("  TIMING_DATA_DIR = ", timing_data_dir)
message("  TIMING_SEED     = ", seed)


#################### Load data ####################

read_rds_safe <- function(path, label) {
  if (!file.exists(path))
    stop(label, " not found at: ", path,
         "\nCheck CLINTALL_DIR ('", dirname(path), "') and re-run Clean_Data.qmd if needed.")
  tryCatch(
    readRDS(path),
    error = function(e) stop(
      "readRDS() failed on ", label, " (", path, ").\n",
      "The file may have been written with save() instead of saveRDS() -- ",
      "re-run the '# Organizing CompRisk' section of Clean_Data.qmd and retry.\n",
      "Original error: ", conditionMessage(e)
    )
  )
}

message("\nLoading data...")
dat_cr     <- read_rds_safe(file.path(clintall_dir, "dat_cr.rds"),     "dat_cr")
dat_design <- read_rds_safe(file.path(clintall_dir, "dat_design.rds"), "dat_design")
message("  dat_cr:     ", nrow(dat_cr),     " x ", ncol(dat_cr))
message("  dat_design: ", nrow(dat_design), " x ", ncol(dat_design))


#################### Align by patient ID ####################

ids_cr     <- dat_cr$sample_id
ids_design <- rownames(dat_design)

if (!setequal(ids_cr, ids_design)) {
  n_only_cr     <- length(setdiff(ids_cr,     ids_design))
  n_only_design <- length(setdiff(ids_design, ids_cr))
  n_overlap     <- length(intersect(ids_cr,   ids_design))
  stop(
    "Patient ID mismatch between dat_cr and dat_design.\n",
    "  dat_cr rows:        ", length(ids_cr),     "\n",
    "  dat_design rows:    ", length(ids_design), "\n",
    "  Overlapping IDs:    ", n_overlap,           "\n",
    "  In dat_cr only:     ", n_only_cr,           "\n",
    "  In dat_design only: ", n_only_design
  )
}

# Align dat_design rows to dat_cr order (do NOT rely on row order)
dat_design_aligned <- dat_design[dat_cr$sample_id, , drop = FALSE]
stopifnot(identical(rownames(dat_design_aligned), dat_cr$sample_id))

N_patients <- nrow(dat_cr)
stopifnot(N_patients == 1309L)
message("\nAligned ", N_patients, " patients.")


#################### Build gene-only predictor pool ####################

clinical_cols <- c("Age", "Sex", "Race", "Ethnicity", "Insurance",
                   "Tumor.Specimen.Type", "CNS.Status", "Risk.Group", "Treatment.Arm")

missing_clinical <- setdiff(clinical_cols, colnames(dat_design_aligned))
if (length(missing_clinical) > 0)
  stop("Expected clinical columns not found in dat_design: ",
       paste(missing_clinical, collapse = ", "))

gene_cols <- setdiff(colnames(dat_design_aligned), clinical_cols)

non_numeric <- gene_cols[
  !vapply(dat_design_aligned[, gene_cols, drop = FALSE], is.numeric, logical(1))
]
if (length(non_numeric) > 0)
  stop("Non-numeric gene column(s) found (first 5): ",
       paste(non_numeric[seq_len(min(5L, length(non_numeric)))], collapse = ", "))

X_all <- as.matrix(dat_design_aligned[, gene_cols, drop = FALSE])

# Drop columns with zero or non-finite variance (constant columns cause
# divide-by-zero when standardizing; fit_ssl_psdh does not standardize internally)
gene_vars <- apply(X_all, 2L, var, na.rm = TRUE)
bad_var   <- !is.finite(gene_vars) | gene_vars <= 0
n_dropped <- sum(bad_var)
if (n_dropped > 0) {
  dropped_names <- gene_cols[bad_var]
  message("Dropping ", n_dropped, " zero/non-finite-variance column(s): ",
          paste(dropped_names[seq_len(min(5L, n_dropped))], collapse = ", "),
          if (n_dropped > 5L) paste0(" ... and ", n_dropped - 5L, " more") else "")
  X_all <- X_all[, !bad_var, drop = FALSE]
}

P_MAX <- ncol(X_all)
message("Gene columns: ", length(gene_cols), " total, ",
        n_dropped, " dropped, ", P_MAX, " remaining (P_MAX = ", P_MAX, ").")


#################### Outcome matrix ####################

y_full <- matrix(
  c(dat_cr$tte, dat_cr$status_relapse),
  nrow  = N_patients,
  ncol  = 2L,
  dimnames = list(dat_cr$sample_id, c("tte", "status_relapse"))
)

event_mix_full <- table(dat_cr$status_relapse)
message("\nEvent mix (full cohort, status_relapse):")
print(event_mix_full)

stopifnot(
  as.integer(event_mix_full["0"]) == 1092L,
  as.integer(event_mix_full["1"]) == 120L,
  as.integer(event_mix_full["2"]) == 97L
)
if (sum(dat_cr$status_relapse == 1L) == 0L)
  stop("No cause-1 (relapse) events found -- cannot proceed.")


#################### Config spec ####################

# p = NA means "all remaining genes"; it is resolved to P_MAX after the
# variance filtering below, because P_MAX can still change at that point.
config_spec <- list(
  list(tag = "n1309-p300",   n = 1309L, p = 300L),
  list(tag = "n1309-p1000",  n = 1309L, p = 1000L),
  list(tag = "n1309-p5000",  n = 1309L, p = 5000L),
  list(tag = "n1309-p10000", n = 1309L, p = 10000L),
  list(tag = "n1309-pall",   n = 1309L, p = NA_integer_)
)

# Smoke-test config: tiny p at the FULL cohort, emitted only when TIMING_SMOKE=1.
# It deliberately uses all patients so that enabling smoke mode cannot change the
# patient-subset variance filter below, and therefore cannot alter the real
# configs' P_MAX.
if (Sys.getenv("TIMING_SMOKE", unset = "0") == "1") {
  config_spec <- c(config_spec, list(list(tag = "smoke", n = 1309L, p = 50L)))
  message("TIMING_SMOKE=1: smoke config (n=1309, p=50) added.")
}

config_ns <- vapply(config_spec, `[[`, integer(1), "n")
stopifnot(all(config_ns >= 1L), all(config_ns <= N_patients))
min_n <- min(config_ns)


#################### Draw nested orderings ####################

message("\nDrawing nested orderings (seed = ", seed, ")...")
set.seed(seed)

# Stratified patient ordering.
# Within each stratum of status_relapse (0/1/2), randomly permute the patients,
# then give each a rank in [0, 1) equal to its within-stratum position divided by
# the stratum size. Sorting by that rank interleaves the strata proportionally, so
# any prefix preserves the full-cohort event mix approximately.
#
# With the current ladder every config takes all N_patients, so this is just a
# permutation and the stratification does nothing. It is retained so a smaller-n
# config can be added later without redesigning this step.
strata_labels <- dat_cr$status_relapse
strata_rank   <- numeric(N_patients)
for (s in c(0L, 1L, 2L)) {
  idx_s <- which(strata_labels == s)
  perm  <- sample(length(idx_s))           # random permutation of 1..ns
  strata_rank[idx_s] <- (perm - 1L) / length(idx_s)
}
patient_order <- order(strata_rank)        # integer indices into dat_cr rows

if (min_n < N_patients) {
  prefix_mix <- table(factor(strata_labels[patient_order[seq_len(min_n)]],
                             levels = c(0L, 1L, 2L)))
  message(min_n, "-patient prefix event mix:")
  message("  ", paste(names(prefix_mix), as.integer(prefix_mix),
                      sep = "=", collapse = "  "))
}

# Gene ordering: single random permutation of all P_MAX columns. Every smaller-p
# subset uses the first p entries, so gene sets are strictly nested across configs.
gene_order <- sample(P_MAX)               # integer indices into columns of X_all

# Second variance filter. A gene that varies across all N_patients can still be
# constant within a smaller patient subset, which would make standardize_cols()
# error at materialize time. This only matters when some config uses fewer than
# all patients; with the current ladder it is skipped and P_MAX is unchanged.
if (min_n < N_patients) {
  var_in_min <- apply(X_all[patient_order[seq_len(min_n)], ], 2L, var)
  bad_in_min <- !is.finite(var_in_min) | var_in_min <= 0
  n_drop2    <- sum(bad_in_min)
  if (n_drop2 > 0) {
    keep_mask  <- !bad_in_min
    # Remap gene_order: drop entries pointing at removed columns, then reindex.
    old_to_new <- cumsum(keep_mask)         # old column index -> new column index
    gene_order <- old_to_new[gene_order[keep_mask[gene_order]]]
    X_all      <- X_all[, keep_mask, drop = FALSE]
    P_MAX      <- ncol(X_all)
    bad_names  <- names(var_in_min)[bad_in_min]
    message("Additional ", n_drop2, " gene(s) with zero variance in the ",
            min_n, "-patient subset dropped: ",
            paste(bad_names[seq_len(min(5L, n_drop2))], collapse = ", "),
            if (n_drop2 > 5L) paste0(" ... and ", n_drop2 - 5L, " more") else "")
    message("  P_MAX updated to ", P_MAX, ".")
  }
} else {
  message("All configs use the full cohort -- patient-subset variance filter skipped.")
}


#################### Resolve configs ####################

configs <- lapply(config_spec, function(cfg) {
  if (is.na(cfg$p)) cfg$p <- P_MAX
  cfg
})

bad_p <- vapply(configs, function(cfg) cfg$p > P_MAX, logical(1))
if (any(bad_p))
  stop("Config(s) request more predictors than exist (P_MAX = ", P_MAX, "): ",
       paste(vapply(configs[bad_p], `[[`, character(1), "tag"), collapse = ", "))

message("\nLadder: ",
        paste(vapply(configs,
                     function(cfg) sprintf("%s (n=%d, p=%d)", cfg$tag, cfg$n, cfg$p),
                     character(1)),
              collapse = "; "))


#################### Standardize helper ####################

standardize_cols <- function(x, tag) {
  mu  <- colMeans(x)
  sds <- apply(x, 2L, sd)
  bad <- !is.finite(sds) | sds <= 0
  if (any(bad))
    stop("[", tag, "] Zero-variance column(s) after subsetting (first 5): ",
         paste(colnames(x)[bad][seq_len(min(5L, sum(bad)))], collapse = ", "),
         if (sum(bad) > 5L) paste0(" ... and ", sum(bad) - 5L, " more") else "")
  sweep(sweep(x, 2L, mu, "-"), 2L, sds, "/")
}


#################### Materialize ####################

message("\nMaterializing configs...")
manifest_rows        <- vector("list", length(configs))
configs_materialized <- vector("list", length(configs))
names(configs_materialized) <- vapply(configs, `[[`, character(1), "tag")

for (i in seq_along(configs)) {
  cfg <- configs[[i]]
  tag <- cfg$tag
  n   <- cfg$n
  p   <- cfg$p

  pat_idx  <- patient_order[seq_len(n)]
  gene_idx <- gene_order[seq_len(p)]

  x_sub <- X_all[pat_idx, gene_idx, drop = FALSE]
  y_sub <- y_full[pat_idx, , drop = FALSE]

  stopifnot(identical(rownames(x_sub), rownames(y_sub)))

  x_std <- standardize_cols(x_sub, tag)

  event_mix_sub <- table(factor(y_sub[, "status_relapse"], levels = c(0L, 1L, 2L)))

  out <- list(
    x    = x_std,
    y    = y_sub,
    meta = list(
      tag         = tag,
      n           = n,
      p           = p,
      seed        = seed,
      event_mix   = event_mix_sub,
      gene_names  = colnames(x_std),
      patient_ids = rownames(x_std),
      created     = Sys.time()
    )
  )

  out_file <- file.path(timing_data_dir, paste0("timing_", tag, ".rds"))
  saveRDS(out, out_file)
  file_size <- file.info(out_file)$size

  message(sprintf("  [%s]  n=%-4d  p=%-6d  mix=%s  %.1f MB",
                  tag, n, p,
                  paste(as.integer(event_mix_sub), collapse = "/"),
                  file_size / 1e6))

  configs_materialized[[tag]] <- out
  manifest_rows[[i]] <- data.frame(
    tag             = tag,
    n               = n,
    p               = p,
    n_cause0        = as.integer(event_mix_sub["0"]),
    n_cause1        = as.integer(event_mix_sub["1"]),
    n_cause2        = as.integer(event_mix_sub["2"]),
    file            = out_file,
    file_size_bytes = file_size,
    stringsAsFactors = FALSE
  )
}


#################### Nesting assertions ####################

message("\nChecking nesting...")

# Derived from `configs` rather than hardcoded tags, so changing the ladder does
# not silently bypass these checks. The smoke config is excluded.
all_tags  <- vapply(configs, `[[`, character(1), "tag")
real_cfgs <- configs[all_tags != "smoke"]

# Patient sets nested by increasing n
ord_n <- order(vapply(real_cfgs, `[[`, integer(1), "n"))
for (k in seq_len(length(ord_n) - 1L)) {
  a <- real_cfgs[[ord_n[k]]]$tag
  b <- real_cfgs[[ord_n[k + 1L]]]$tag
  if (!all(configs_materialized[[a]]$meta$patient_ids %in%
           configs_materialized[[b]]$meta$patient_ids))
    stop("Patient sets are not nested: ", a, " is not a subset of ", b, ".")
}

# Gene sets nested by increasing p
ord_p <- order(vapply(real_cfgs, `[[`, integer(1), "p"))
for (k in seq_len(length(ord_p) - 1L)) {
  a <- real_cfgs[[ord_p[k]]]$tag
  b <- real_cfgs[[ord_p[k + 1L]]]$tag
  if (!all(configs_materialized[[a]]$meta$gene_names %in%
           configs_materialized[[b]]$meta$gene_names))
    stop("Gene sets are not nested: ", a, " is not a subset of ", b, ".")
}

# Every config has the dimensions it claims, and x/y rows align
for (cfg in configs) {
  obj <- configs_materialized[[cfg$tag]]
  stopifnot(
    nrow(obj$x) == cfg$n,
    ncol(obj$x) == cfg$p,
    identical(rownames(obj$x), rownames(obj$y))
  )
}

message("  All nesting assertions passed (", length(real_cfgs), " configs checked).")


#################### Manifest ####################

manifest      <- do.call(rbind, manifest_rows)
manifest_file <- file.path(timing_data_dir, "timing_manifest.csv")
write.csv(manifest, manifest_file, row.names = FALSE)

message("\nManifest:")
print(manifest[, c("tag", "n", "p", "n_cause0", "n_cause1", "n_cause2",
                   "file_size_bytes")])
message("\nDone. Manifest written to: ", manifest_file)
