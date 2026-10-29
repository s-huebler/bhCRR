# sim/analysis/R/read_results.r -- load one scenario's parsed results
#
# OWNS
#   read_scenario_results(scenario, root, p = NULL)
#                       read every parsed table for one scenario and one p,
#                       plus the scenario.env manifest and the scenario
#                       definition file, into one list
#   read_manifest(path) parse scenario.env into $meta (the "# key: value"
#                       header lines) and $settings (the KEY=VALUE lines)
#
# INPUTS
#   scenario  character(1)  scenario slug, e.g. "equalblocks.01"
#   root      character(1)  bhCRR repository root
#   p         integer(1) or NULL
#                           which p-group to read. NULL picks the only one
#                           present and errors (listing the choices) if
#                           there are several. Per-predictor tables are never
#                           pooled across p because the block map depends on p.
#
# RETURNS  list of class "sim_results"
#   $scenario_name, $p, $parsed_dir
#   $betas            one row per predictor x model x run (parse.r C1)
#   $performance      one row per model x run (parse.r C2)
#   $by_block         one row per model x block x run
#   $false_positives  pooled per model x block
#   $manifest         list(meta, settings) from scenario.env, or NULL
#   $design           the `scenario` list from sim/scenarios/<scenario>.R,
#                     or NULL when the file is absent
#   $design_file      path to that file (NA when absent)
#   $files            paths read
#
# DEPENDS ON  base R only. bhCRR is NOT loaded: analysis works on the parsed
#             tables alone.

read_scenario_results <- function(scenario, root, p = NULL) {
  parsed_dir <- file.path(root, "sim", "results", scenario, "parsed")
  if (!dir.exists(parsed_dir))
    stop("No parsed results directory: ", parsed_dir,
         "\nRun parse_results() for this scenario first.")

  scen_rx <- gsub(".", "\\.", scenario, fixed = TRUE)
  perf_files <- list.files(parsed_dir,
                           pattern = sprintf("^model_performance_p(\\d+)_%s\\.csv$",
                                             scen_rx))
  p_avail <- sort(as.integer(sub("^model_performance_p(\\d+)_.*$", "\\1",
                                 perf_files)))
  if (!length(p_avail))
    stop("No model_performance_p<P>_", scenario, ".csv in ", parsed_dir)
  if (is.null(p)) {
    if (length(p_avail) > 1L)
      stop("Several p-groups in ", parsed_dir, ": ",
           paste(p_avail, collapse = ", "), ". Pass p = one of them.")
    p <- p_avail
  } else if (!p %in% p_avail) {
    stop("p = ", p, " not found; available: ", paste(p_avail, collapse = ", "))
  }

  f <- function(stem, ext) file.path(parsed_dir,
                                     sprintf("%s_p%d_%s.%s", stem, p, scenario, ext))
  files <- c(betas_rdata     = f("betas", "Rdata"),
             betas_csv       = f("betas", "csv"),
             performance     = f("model_performance", "csv"),
             by_block        = f("by_block", "csv"),
             false_positives = f("false_positives", "csv"),
             manifest        = file.path(parsed_dir, "scenario.env"),
             design          = file.path(root, "sim", "scenarios",
                                         paste0(scenario, ".R")))

  read_csv_or_null <- function(path) {
    if (!file.exists(path)) {
      warning("Missing parsed file: ", basename(path))
      return(NULL)
    }
    utils::read.csv(path, stringsAsFactors = FALSE)
  }

  # Betas: the .Rdata is ~40x smaller than the CSV, so prefer it.
  betas <- NULL
  if (file.exists(files[["betas_rdata"]])) {
    e <- new.env(parent = emptyenv())
    nm <- load(files[["betas_rdata"]], envir = e)
    betas <- e[[nm[1L]]]
  } else if (file.exists(files[["betas_csv"]])) {
    betas <- utils::read.csv(files[["betas_csv"]], stringsAsFactors = FALSE)
  } else {
    warning("No betas table for p = ", p)
  }

  design <- NULL
  if (file.exists(files[["design"]])) {
    e <- new.env()
    sys.source(files[["design"]], envir = e)
    design <- e$scenario
  }

  out <- list(
    scenario_name   = scenario,
    p               = p,
    parsed_dir      = parsed_dir,
    betas           = betas,
    performance     = read_csv_or_null(files[["performance"]]),
    by_block        = read_csv_or_null(files[["by_block"]]),
    false_positives = read_csv_or_null(files[["false_positives"]]),
    manifest        = if (file.exists(files[["manifest"]]))
                        read_manifest(files[["manifest"]]) else NULL,
    design          = design,
    design_file     = if (file.exists(files[["design"]])) files[["design"]] else NA,
    files           = files[file.exists(files)]
  )
  class(out) <- c("sim_results", "list")
  out
}

read_manifest <- function(path) {
  lines <- readLines(path, warn = FALSE)
  hdr <- grep("^#\\s*[A-Za-z_]+:", lines, value = TRUE)
  meta <- stats::setNames(
    as.list(trimws(sub("^#\\s*[A-Za-z_]+:\\s*", "", hdr))),
    sub("^#\\s*([A-Za-z_]+):.*$", "\\1", hdr))
  kv <- grep("^[A-Z0-9_]+=", lines, value = TRUE)
  settings <- stats::setNames(as.list(sub("^[^=]+=", "", kv)),
                              sub("=.*$", "", kv))
  list(meta = meta, settings = settings)
}

print.sim_results <- function(x, ...) {
  cat(sprintf("<sim_results> %s  p = %d\n", x$scenario_name, x$p))
  if (!is.null(x$performance))
    cat(sprintf("  runs: %d   models: %s\n",
                length(unique(x$performance$run_number)),
                paste(unique(x$performance$Model), collapse = ", ")))
  cat("  parsed_dir:", x$parsed_dir, "\n")
  invisible(x)
}
