# sim/timing/collect_timing.r
#
# Collect all per-task res_*.rds files into one CSV.
# Missing (p, rep) cells are filled from sacct output (passed as a file).
# A TIMEOUT row is treated as a censored observation at the time limit.
# If a task ran more than once (retry), the latest successful row wins.
#
# Environment variables:
#   TIMING_OUT_DIR   directory containing res_p*_rep*.rds files  (required)
#   TIMING_TAG       label written into the output filename       (required)
#   GRID_FILE        path to grid.tsv                            (required)
#   REPS             number of replicates per p value            (required)
#   SACCT_FILE       path to sacct output file; empty = skip     (optional)
#   ARRAY_JOB_ID     numeric SLURM job id (informational only)   (optional)
#
# Output:
#   <TIMING_OUT_DIR>/clintall_timing_<TIMING_TAG>.csv

.req_env <- function(name) {
  v <- Sys.getenv(name, unset = "")
  if (!nzchar(v)) stop("Required environment variable not set: ", name)
  v
}
.opt_env <- function(name, default = "") {
  v <- Sys.getenv(name, unset = "")
  if (nzchar(v)) v else default
}

# ---------------------------------------------------------------------------
# Helpers (defined before use)
# ---------------------------------------------------------------------------

# Parse "D-HH:MM:SS", "HH:MM:SS", or "MM:SS" to seconds.
.parse_hms <- function(s) {
  if (is.null(s) || is.na(s) || !nzchar(s)) return(NA_real_)
  s <- trimws(s)
  parts_day <- strsplit(s, "-")[[1]]
  days <- 0
  if (length(parts_day) == 2L) {
    days <- as.numeric(parts_day[1])
    s    <- parts_day[2]
  }
  parts <- suppressWarnings(as.numeric(strsplit(s, ":")[[1]]))
  if (length(parts) == 3L)      days * 86400 + parts[1]*3600 + parts[2]*60 + parts[3]
  else if (length(parts) == 2L) days * 86400 + parts[1]*60   + parts[2]
  else NA_real_
}

# Parse sacct MaxRSS (e.g. "512K", "1.2G") to MB.
.parse_maxrss <- function(s) {
  if (is.null(s) || is.na(s) || !nzchar(s) || s == "0") return(NA_real_)
  s    <- trimws(s)
  mult <- if (grepl("[Kk]$", s)) 1/1024
          else if (grepl("[Mm]$", s)) 1
          else if (grepl("[Gg]$", s)) 1024
          else 1/1024    # assume kB if no suffix
  val  <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", s)))
  val * mult
}

# rbind two data frames, adding missing columns as NA.
.bind_frames <- function(a, b) {
  if (nrow(a) == 0L) return(b)
  if (nrow(b) == 0L) return(a)
  all_cols <- union(names(a), names(b))
  for (col in all_cols) {
    if (!col %in% names(a)) a[[col]] <- NA
    if (!col %in% names(b)) b[[col]] <- NA
  }
  rbind(a[, all_cols, drop = FALSE], b[, all_cols, drop = FALSE])
}

# ---------------------------------------------------------------------------
# 0. Read env
# ---------------------------------------------------------------------------
timing_out_dir <- .req_env("TIMING_OUT_DIR")
timing_tag     <- .req_env("TIMING_TAG")
grid_file      <- .req_env("GRID_FILE")
reps           <- as.integer(.req_env("REPS"))
sacct_file     <- .opt_env("SACCT_FILE",   "")
array_job_id   <- .opt_env("ARRAY_JOB_ID", "")

cat(sprintf("[collect] tag=%s  dir=%s\n", timing_tag, timing_out_dir))

# ---------------------------------------------------------------------------
# 1. Read grid to build the expected (p, rep) set
# ---------------------------------------------------------------------------
grid <- read.table(grid_file, header = TRUE, sep = "\t")
G    <- nrow(grid)
expected <- expand.grid(p_index = seq_len(G), rep = seq_len(reps),
                        stringsAsFactors = FALSE)
expected$p <- grid$p[expected$p_index]
expected   <- expected[, c("p", "rep")]

cat(sprintf("[collect] Expecting %d result files (G=%d, reps=%d)\n",
            nrow(expected), G, reps))

# ---------------------------------------------------------------------------
# 2. Read all res_*.rds files
# ---------------------------------------------------------------------------
res_files <- list.files(timing_out_dir,
                        pattern = "^res_p[0-9]+_rep[0-9]+\\.rds$",
                        full.names = TRUE)
cat(sprintf("[collect] Found %d result files\n", length(res_files)))

read_one <- function(path) {
  tryCatch({
    r <- readRDS(path)
    as.data.frame(lapply(r, function(v) {
      if (length(v) != 1L) NA
      else if (inherits(v, "POSIXct")) format(v, "%Y-%m-%d %H:%M:%S")
      else v
    }), stringsAsFactors = FALSE)
  }, error = function(e) {
    message("WARNING: could not read ", path, ": ", conditionMessage(e))
    NULL
  })
}

rows <- lapply(res_files, read_one)
rows <- Filter(Negate(is.null), rows)

results <- if (length(rows) > 0L) do.call(rbind, rows) else data.frame()

# ---------------------------------------------------------------------------
# 3. Retry deduplication: keep latest successful row per (p, rep)
# ---------------------------------------------------------------------------
if (nrow(results) > 0L) {
  results$timestamp <- as.POSIXct(results$timestamp, format = "%Y-%m-%d %H:%M:%S")
  results <- results[order(results$p, results$rep, results$timestamp), ]

  .is_success <- function(err_col) {
    is.na(err_col) | err_col == "NA" | err_col == "none"
  }
  success <- .is_success(results$err)

  keep <- logical(nrow(results))
  for (i in seq_len(nrow(results))) {
    pi  <- results$p[i]
    ri  <- results$rep[i]
    dup <- which(keep &
                 !is.na(results$p)   & results$p   == pi &
                 !is.na(results$rep) & results$rep == ri)
    if (length(dup) == 0L) {
      keep[i] <- TRUE
    } else if (success[i]) {
      keep[dup] <- FALSE
      keep[i]   <- TRUE
    }
  }
  results <- results[keep, , drop = FALSE]
  cat(sprintf("[collect] After deduplication: %d rows\n", nrow(results)))
}

# ---------------------------------------------------------------------------
# 4. Fill missing (p, rep) cells from sacct
# ---------------------------------------------------------------------------
present <- if (nrow(results) > 0L) {
  paste(results$p, results$rep, sep = "_")
} else {
  character(0)
}

missing_grid <- expected[
  !paste(expected$p, expected$rep, sep = "_") %in% present, ,
  drop = FALSE]

if (nrow(missing_grid) > 0L) {
  cat(sprintf("[collect] %d (p,rep) cells missing from result files\n",
              nrow(missing_grid)))

  sacct_df <- NULL

  if (nzchar(sacct_file) && file.exists(sacct_file)) {
    cat(sprintf("[collect] Parsing sacct file: %s\n", sacct_file))
    sacct_raw <- tryCatch(
      read.table(sacct_file, header = TRUE, sep = "|",
                 stringsAsFactors = FALSE, fill = TRUE, comment.char = ""),
      error = function(e) {
        message("WARNING: could not parse sacct file: ", conditionMessage(e))
        NULL
      }
    )

    if (!is.null(sacct_raw) && nrow(sacct_raw) > 0L) {
      sacct_raw$task_id <- suppressWarnings(
        as.integer(sub(".*_([0-9]+)$", "\\1", sacct_raw[, 1L]))
      )
      sacct_raw <- sacct_raw[!is.na(sacct_raw$task_id), , drop = FALSE]

      fill_rows <- vector("list", nrow(missing_grid))
      for (i in seq_len(nrow(missing_grid))) {
        mp       <- missing_grid$p[i]
        mrep     <- missing_grid$rep[i]
        p_idx    <- which(grid$p == mp)
        tid_exp  <- (mrep - 1L) * G + (p_idx - 1L)

        s <- sacct_raw[sacct_raw$task_id == tid_exp, , drop = FALSE]
        if (nrow(s) > 0L) s <- s[nrow(s), , drop = FALSE]   # last attempt

        state_val     <- if (nrow(s) > 0L && "State"     %in% names(s)) s$State     else "MISSING"
        elapsed_val   <- if (nrow(s) > 0L && "Elapsed"   %in% names(s)) s$Elapsed   else NA_character_
        timelimit_val <- if (nrow(s) > 0L && "Timelimit" %in% names(s)) s$Timelimit else NA_character_
        maxrss_val    <- if (nrow(s) > 0L && "MaxRSS"    %in% names(s)) s$MaxRSS    else NA_character_

        timeout_secs  <- if (identical(state_val, "TIMEOUT")) .parse_hms(timelimit_val) else NA_real_

        fill_rows[[i]] <- data.frame(
          tag                = timing_tag,
          task_id            = tid_exp,
          p                  = mp,
          rep                = mrep,
          n                  = NA_integer_,
          n_events1          = NA_integer_,
          init_method        = NA_character_,
          s0                 = NA_real_,
          s1                 = NA_real_,
          t_read             = NA_real_,
          t_init             = NA_real_,
          t_em               = NA_real_,
          t_total            = timeout_secs,
          iterations         = NA_integer_,
          conv               = NA,
          nnz                = NA_integer_,
          init_lambda        = NA_real_,
          init_df            = NA_integer_,
          err                = state_val,
          peak_rss_mb        = .parse_maxrss(maxrss_val),
          host               = NA_character_,
          cpu_model          = NA_character_,
          slurm_job_id       = NA_character_,
          time_limit         = timelimit_val,
          R_version          = NA_character_,
          fastcmprsk_version = NA_character_,
          bhCRR_git_hash     = NA_character_,
          timestamp          = NA_character_,
          sacct_state        = state_val,
          sacct_elapsed      = elapsed_val,
          stringsAsFactors   = FALSE
        )
      }

      sacct_df <- do.call(rbind, Filter(Negate(is.null), fill_rows))
      if (!is.null(sacct_df) && nrow(sacct_df) > 0L) {
        results <- .bind_frames(results, sacct_df)
        cat(sprintf("[collect] Added %d sacct rows\n", nrow(sacct_df)))
      }
    }
  } else {
    cat("[collect] No sacct file provided or file not found; missing cells left as absent.\n")
  }
}

# ---------------------------------------------------------------------------
# 5. Write CSV and print summary
# ---------------------------------------------------------------------------
out_csv <- file.path(timing_out_dir,
                     sprintf("clintall_timing_%s.csv", timing_tag))
write.csv(results, out_csv, row.names = FALSE)
cat(sprintf("[collect] Wrote %d rows to: %s\n", nrow(results), out_csv))

if (nrow(results) > 0L && "p" %in% names(results)) {
  cat("\nSummary by p (t_total seconds):\n")
  p_vals <- sort(unique(results$p[!is.na(results$p)]))
  for (pv in p_vals) {
    sub_p <- results[!is.na(results$p) & results$p == pv, , drop = FALSE]
    ttot  <- sub_p$t_total[!is.na(sub_p$t_total)]
    n_err <- sum(!is.na(sub_p$err) & sub_p$err != "NA" & sub_p$err != "none",
                 na.rm = TRUE)
    cat(sprintf("  p=%6d  n=%2d  t_total: min=%.1f  med=%.1f  max=%.1f  errors=%d\n",
                pv, nrow(sub_p),
                if (length(ttot)) min(ttot)    else NA,
                if (length(ttot)) median(ttot) else NA,
                if (length(ttot)) max(ttot)    else NA,
                n_err))
  }
}

cat("[collect] Done.\n")
