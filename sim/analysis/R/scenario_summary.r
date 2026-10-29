# sim/analysis/R/scenario_summary.r -- describe the scenario and the runs
#
# OWNS
#   scenario_design_notes(res)  the DESIGN section of the scenario file header,
#                               as a character(1) paragraph (NULL if absent)
#   scenario_overview(res)      two-column Setting / Value table of the DGP
#   block_design_table(res)     one row per correlation block
#   run_provenance(res, root)   when and from which commit the results came
#   run_health(res, models)     per model: runs, converged, usable, timing
#   tuning_summary(res, models) per model with tuning columns (s0, s1,
#                               n_attempts, all_pass); today only bhCRR has
#                               them, but nothing here names bhCRR
#
# All functions take the list returned by read_scenario_results(). None of
# them names a model id; model display goes through the registry in models.r.

# ---------------------------------------------------------------------------
# Scenario description
# ---------------------------------------------------------------------------

scenario_design_notes <- function(res) {
  if (is.na(res$design_file)) return(NULL)
  lines <- readLines(res$design_file, warn = FALSE)
  # Header comment block only (up to the first non-comment, non-blank line)
  first_code <- which(!grepl("^\\s*(#|$)", lines))[1L]
  if (is.na(first_code)) first_code <- length(lines) + 1L
  hdr <- lines[seq_len(first_code - 1L)]
  # Section headings sit one space after the "#" ("# DESIGN"); body text is
  # indented further ("#   High dimensional ..."). The DESIGN section runs to
  # the next heading.
  is_head <- grepl("^#\\s?\\S", hdr)
  start <- grep("^#\\s*DESIGN\\s*$", hdr)[1L]
  if (is.na(start)) return(NULL)
  nxt <- which(is_head)
  end <- min(c(nxt[nxt > start], length(hdr) + 1L)) - 1L
  body <- trimws(sub("^#", "", hdr[(start + 1L):end]))
  paste(body[nzchar(body)], collapse = " ")
}

scenario_overview <- function(res) {
  d <- res$design
  s <- res$manifest$settings
  if (is.null(d) && is.null(s)) stop("No scenario file and no scenario.env.")
  num <- function(x) as.numeric(strsplit(x, ",")[[1L]])

  beta1 <- if (!is.null(d)) d$beta1_active else num(s$BETA1_ACTIVE)
  beta2 <- if (!is.null(d)) d$beta2_active else s$BETA2_ACTIVE
  get <- function(key) {
    if (!is.null(d[[key]])) return(d[[key]])
    s[[toupper(key)]]
  }
  beta2_txt <- if (identical(beta2, "neg_beta1")) "beta2 = -beta1" else
    sprintf("%d nonzero, |beta2| in [%.2f, %.2f]",
            sum(beta2 != 0), min(abs(beta2[beta2 != 0])),
            max(abs(beta2[beta2 != 0])))

  perf <- res$performance
  n_runs <- if (!is.null(perf)) length(unique(perf$run_number)) else NA
  eval_t <- if (!is.null(perf)) unique(perf$eval_time) else NA

  data.frame(
    Setting = c("Scenario", "Predictors (p)", "Observations (n)",
                "Completed runs", "Active predictors (cause 1)",
                "Cause 1 effects",
                "Cause 2 effects", "Censoring", "Observation window (u)",
                "Latent type", "Zero gap target (autotune beta_min)",
                "C-index horizon (eval_time)"),
    Value = c(
      res$scenario_name,
      format(res$p),
      format(get("nobs")),
      format(n_runs),
      sprintf("%d (%d positive, %d negative)", sum(beta1 != 0),
              sum(beta1 > 0), sum(beta1 < 0)),
      sprintf("|beta1| in [%.2f, %.2f]", min(abs(beta1[beta1 != 0])),
              max(abs(beta1[beta1 != 0]))),
      beta2_txt,
      sprintf("Exponential, rate %s", format(get("cens_rate"))),
      sprintf("[%s, %s]", format(get("u_min")), format(get("u_max"))),
      format(get("latent_type")),
      format(get("zero_gap_target")),
      paste(format(eval_t), collapse = ", ")
    ),
    stringsAsFactors = FALSE
  )
}

block_design_table <- function(res) {
  d <- res$design
  if (is.null(d)) {
    s <- res$manifest$settings
    kv <- function(x) {
      parts <- strsplit(strsplit(x, ",")[[1L]], "=")
      stats::setNames(as.numeric(vapply(parts, `[`, "", 2L)),
                      vapply(parts, `[`, "", 1L))
    }
    d <- list(block_props = kv(s$BLOCK_PROPS), block_rho = kv(s$BLOCK_RHO),
              latent_q = kv(s$LATENT_Q),
              beta1_active = as.numeric(strsplit(s$BETA1_ACTIVE, ",")[[1L]]),
              active_block = strsplit(s$ACTIVE_BLOCK, ",")[[1L]])
  }
  blocks <- names(d$block_props)

  # Realised block sizes come from the betas table (one run is enough: the
  # block map is deterministic within a p-group).
  n_pred <- if (!is.null(res$betas)) {
    b1 <- res$betas[res$betas$run_number == res$betas$run_number[1L] &
                      res$betas$Model == res$betas$Model[1L], ]
    as.integer(table(factor(b1$Block, levels = blocks)))
  } else round(d$block_props * res$p)

  act <- lapply(blocks, function(b) d$beta1_active[d$active_block == b])
  data.frame(
    Block         = blocks,
    `Within-block rho` = unname(d$block_rho[blocks]),
    `Latent q`    = unname(d$latent_q[blocks]),
    Proportion    = unname(d$block_props[blocks]),
    Predictors    = n_pred,
    Active        = vapply(act, length, integer(1)),
    `Active |beta1| range` = vapply(act, function(a)
      if (length(a)) sprintf("%.2f to %.2f", min(abs(a)), max(abs(a))) else "",
      character(1)),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Provenance and run health
# ---------------------------------------------------------------------------

run_provenance <- function(res, root) {
  m <- res$manifest$meta
  out <- data.frame(
    Item  = c("Scenario manifest created", "Manifest git commit",
              "Parsed tables last written", "Run numbers"),
    Value = c(m$created %||% NA, m$git %||% NA,
              format(max(file.mtime(res$files[!grepl("scenarios",
                                                      res$files)])),
                     "%Y-%m-%d %H:%M"),
              .run_range(res$performance$run_number)),
    stringsAsFactors = FALSE
  )
  rel <- file.path("sim", "results", res$scenario_name, "parsed")
  log <- tryCatch(
    suppressWarnings(system2("git", shQuote(c("-C", root, "log",
                                              "--date=format:%Y-%m-%d %H:%M",
                                              "--pretty=format:%h|%ad|%s",
                                              "--", rel)),
                             stdout = TRUE, stderr = FALSE)),
    error = function(e) character(0))
  commits <- if (length(log) && is.null(attr(log, "status"))) {
    parts <- strsplit(log, "|", fixed = TRUE)
    data.frame(Commit  = vapply(parts, `[`, "", 1L),
               Date    = vapply(parts, `[`, "", 2L),
               Message = vapply(parts, function(x)
                 paste(x[-(1:2)], collapse = "|"), ""),
               stringsAsFactors = FALSE)
  } else NULL
  list(summary = out, commits = commits)
}

.run_range <- function(r) {
  r <- sort(unique(r))
  if (!length(r)) return(NA_character_)
  gaps <- setdiff(seq(min(r), max(r)), r)
  sprintf("%d to %d (%d runs%s)", min(r), max(r), length(r),
          if (length(gaps)) paste0("; missing ", paste(gaps, collapse = ", "))
          else "")
}

# A model-run is "usable" when it converged and returned coefficients. Failed
# fits come through parse.r as converged == FALSE with NA estimates and an
# all-zero confusion matrix, so they must be removed before averaging.
.usable_runs <- function(perf) {
  conv <- is.na(perf$converged) | as.logical(perf$converged)
  counted <- (perf$tp + perf$fp + perf$tn + perf$fn) > 0
  conv & counted
}

run_health <- function(res, models) {
  perf <- res$performance
  models <- resolve_models(models, perf$Model)
  perf$usable <- .usable_runs(perf)
  perf$selected <- perf$tp + perf$fp
  setup <- ifelse(is.na(perf$setup_seconds), 0, perf$setup_seconds)
  perf$total_seconds <- perf$fit_seconds + setup

  out <- do.call(rbind, lapply(split(perf, perf$Model), function(d) {
    u <- d[d$usable, ]
    data.frame(
      Model           = d$Model[1L],
      Runs            = nrow(d),
      Converged       = sum(as.logical(d$converged) %in% TRUE),
      Usable          = nrow(u),
      `Median selected` = if (nrow(u)) stats::median(u$selected) else NA,
      `Median fit (min)` = round(stats::median(d$total_seconds) / 60, 2),
      check.names = FALSE, stringsAsFactors = FALSE)
  }))
  out <- label_models(out, models)
  out <- out[order(out$Model), ]
  rownames(out) <- NULL
  out
}

tuning_summary <- function(res, models) {
  perf <- res$performance
  tune_cols <- intersect(c("s0", "s1", "n_attempts", "all_pass", "iterations"),
                         names(perf))
  has <- rowSums(!is.na(perf[, tune_cols, drop = FALSE])) > 0
  perf <- perf[has, ]
  if (!nrow(perf)) return(NULL)
  models <- resolve_models(models, perf$Model)
  rng <- function(x, d = 4) if (all(is.na(x))) NA_character_ else
    sprintf("%s (%s to %s)", signif(stats::median(x, na.rm = TRUE), d),
            signif(min(x, na.rm = TRUE), d), signif(max(x, na.rm = TRUE), d))
  tab <- function(x) paste(sprintf("%s: %d", names(table(x)), table(x)),
                           collapse = ", ")
  out <- do.call(rbind, lapply(split(perf, perf$Model), function(d) data.frame(
    Model                = d$Model[1L],
    `s0 median (range)`  = rng(d$s0),
    `s1 median (range)`  = rng(d$s1),
    `Autotune attempts`  = tab(d$n_attempts),
    `Grid checks passed` = sprintf("%d / %d", sum(as.logical(d$all_pass) %in% TRUE),
                                   nrow(d)),
    `EM iterations median (range)` = rng(d$iterations, 3),
    check.names = FALSE, stringsAsFactors = FALSE)))
  out <- label_models(out, models)
  rownames(out) <- NULL
  out
}
