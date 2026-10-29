# sim/analysis/R/load_analysis.r -- bring the analysis modules into scope
#
# OWNS
#   find_repo_root(start = getwd())
#                       walk up from `start` to the directory whose DESCRIPTION
#                       says Package: bhCRR. Works whether a .qmd is rendered
#                       (working directory = the .qmd's folder) or run chunk by
#                       chunk in Positron (working directory = repo root).
#   load_analysis(root = find_repo_root(), envir = globalenv())
#                       check the required packages, then source every module
#                       in sim/analysis/R/ in dependency order
#
# This is the only file in sim/analysis/ that contains source() paths. The
# analysis layer does NOT load bhCRR or sim/R/: it reads parsed tables only,
# so a report renders in seconds and on a machine without a compiled package.
#
# MODULES (dependency order)
#   models.r            model registry: ids, labels, colours
#   read_results.r      read_scenario_results()
#   scenario_summary.r  scenario description, provenance, run health
#   bias.r              estimate / bias / percent bias summaries and plots
#   selection.r         confusion matrix and selection metrics

find_repo_root <- function(start = getwd()) {
  env <- Sys.getenv("REPO_ROOT")
  if (nzchar(env)) return(normalizePath(env))
  dir <- normalizePath(start)
  repeat {
    desc <- file.path(dir, "DESCRIPTION")
    if (file.exists(desc) &&
        any(grepl("^Package:\\s*bhCRR\\s*$", readLines(desc, warn = FALSE))))
      return(dir)
    parent <- dirname(dir)
    if (identical(parent, dir))
      stop("Could not find the bhCRR repository root above ", start,
           ". Set REPO_ROOT or open the bhCRR project.")
    dir <- parent
  }
}

load_analysis <- function(root = find_repo_root(), envir = globalenv()) {
  need <- c("dplyr", "ggplot2", "stringr", "knitr")
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss))
    stop("Missing packages for sim/analysis: ", paste(miss, collapse = ", "),
         ". Install with renv::install().")

  if (!exists("%||%", envir = baseenv(), inherits = FALSE))
    assign("%||%", function(a, b) if (is.null(a)) b else a, envir = envir)

  dir <- file.path(root, "sim", "analysis", "R")
  modules <- c("models.r", "read_results.r", "scenario_summary.r",
               "bias.r", "selection.r")
  paths <- file.path(dir, modules)
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stop("Missing analysis module(s): ",
                            paste(basename(missing), collapse = ", "))
  for (p in paths) sys.source(p, envir = envir)
  invisible(paths)
}
