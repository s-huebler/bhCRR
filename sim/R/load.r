# sim/R/load.r -- bring bhCRR and all sim/R modules into scope
#
# This is the ONLY file in sim/ permitted to contain source() paths into the
# package or into sim/R/. Every other sim/ file assumes bhCRR is already
# attached and the sim modules are already sourced.
#
# OWNS
#   repo_root()
#                       locate the bhCRR repository root; resolution order:
#                       (a) REPO_ROOT env var if set and non-empty,
#                       (b) walk up from getwd() looking for DESCRIPTION with
#                           Package: bhCRR,
#                       (c) stop() naming both remedies
#                       Validates every candidate: a DESCRIPTION that exists
#                       but names a different package is an error, not a miss.
#
#   check_compiled(root)
#                       internal; verify src/bhCRR<dynlib.ext> exists under
#                       root; called by load_bhcrr() in cluster mode; errors
#                       with the pkgbuild::compile_dll() remedy if missing
#
#   load_bhcrr(mode = Sys.getenv("BHCRR_LOAD", "dev"), force = FALSE)
#                       attach bhCRR via pkgload::load_all(); idempotent:
#                       returns FALSE invisibly without reloading if
#                       "package:bhCRR" is already on the search path and
#                       force = FALSE; TRUE if it loaded
#                       mode = "dev"     compile = TRUE
#                       mode = "cluster" compile = FALSE (SO pre-built;
#                                        check_compiled() is called first)
#                       export_all = FALSE in both modes (see note below)
#
#   load_sim(mode = Sys.getenv("BHCRR_LOAD", "dev"), force = FALSE)
#                       call load_bhcrr(mode, force), then source the seven
#                       sim/R modules in dependency order into the global
#                       environment; return invisibly the character vector of
#                       file paths sourced
#
# INPUTS
#   mode   character(1)  "dev" or "cluster"
#   force  logical(1)    TRUE forces a reload even if bhCRR is already attached
#   root   character(1)  repo root path, passed to check_compiled()
#
# RETURNS
#   repo_root():    character(1) absolute path to the repository root
#   check_compiled(): called for its side effect; returns invisible(NULL)
#   load_bhcrr():   invisible(logical(1)) -- TRUE if bhCRR was (re)loaded,
#                   FALSE if it was already attached and force = FALSE
#   load_sim():     invisible(character(7)) -- paths of the sim/R modules
#                   sourced, in dependency order
#
# NOTE ON export_all = FALSE
#   pkgload::load_all() defaults to export_all = TRUE, which exposes every
#   internal object (functions with dot-prefixed names, un-@exported helpers,
#   etc.) on the package environment. library() does not. Using export_all =
#   FALSE in both modes keeps the sim layer honest against the package's actual
#   public API, so it cannot silently come to depend on something like
#   .cv_fold_path that a real install would never export. Nothing in the planned
#   sim layer needs an internal: fit.r uses fit_ssl_psdh, bhcrr_autotune, and
#   bhcrr_cv_control; score.r uses wolbers_c; comparators use fastcmprsk directly.
#
# NOTE ON renv
#   renv is activated by source("renv/activate.R") in the repo .Rprofile.
#   sim/chpc/run_array.slurm cds to REPO_ROOT before invoking Rscript, so
#   renv is already active before this file is sourced. Do not call
#   renv::load() or renv::activate() here.
#
# NOTE ON BHCRR_LOAD and cluster mode
#   BHCRR_LOAD is not yet exported by sim/chpc/config.sh, so cluster mode is
#   currently reachable only by passing mode = "cluster" explicitly. Migration
#   step 12 (sim/chpc wrappers) wires up the env var and adds the pre-array
#   compile step to sim/chpc/run_sim.sh.
#
# NOTE ON force = TRUE
#   After editing anything in R/ or src/, call load_bhcrr(force = TRUE) (or
#   equivalently load_sim(force = TRUE)) to pick up the changes. The idempotence
#   guard skips reload when bhCRR is already on the search path, so edits made
#   during an interactive session are invisible until a forced reload.
#
# DEPENDS ON
#   pkgload  (must be installed; load.r checks and errors if absent)
#   pkgbuild (checked only in cluster mode, for check_compiled())
#
# REPLACES
#   The source() calls scattered across sim/_legacy/code/*.R (every legacy
#   script sourced individual R/ files by relative path; those paths broke
#   whenever the working directory changed).

# ---------------------------------------------------------------------------
# repo_root -- internal
# ---------------------------------------------------------------------------

repo_root <- function() {

  .validate <- function(path) {
    desc <- file.path(path, "DESCRIPTION")
    if (!file.exists(desc)) return(FALSE)
    lines <- readLines(desc, warn = FALSE)
    pkg_line <- grep("^Package:", lines, value = TRUE)
    length(pkg_line) > 0 && trimws(sub("^Package:\\s*", "", pkg_line[1])) == "bhCRR"
  }

  # (a) environment variable
  env_root <- Sys.getenv("REPO_ROOT")
  if (nzchar(env_root)) {
    if (.validate(env_root))
      return(normalizePath(env_root, mustWork = TRUE))
    stop(
      "REPO_ROOT is set to '", env_root, "' but that directory does not contain ",
      "a DESCRIPTION file with Package: bhCRR. ",
      "Fix REPO_ROOT or unset it and set your working directory inside the repo."
    )
  }

  # (b) walk up from getwd()
  path <- normalizePath(getwd(), mustWork = TRUE)
  repeat {
    if (.validate(path)) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) break
    path <- parent
  }

  # (c) give up
  stop(
    "Cannot locate the bhCRR repository root. Either:\n",
    "  - export REPO_ROOT=/path/to/bhCRR  (sim/chpc/config.sh does this on CHPC)\n",
    "  - setwd() to any directory inside the bhCRR repository before sourcing load.r"
  )
}

# ---------------------------------------------------------------------------
# check_compiled -- internal; called by load_bhcrr() in cluster mode
# ---------------------------------------------------------------------------

check_compiled <- function(root) {
  # Array tasks share one REPO_ROOT on a shared filesystem. If compile = TRUE
  # were used in cluster mode, concurrent tasks would race on the same
  # src/*.o and src/bhCRR.so, corrupting the shared object. Cluster mode
  # therefore always uses a pre-built SO and never compiles.
  so_path <- file.path(root, "src", paste0("bhCRR", .Platform$dynlib.ext))
  if (!file.exists(so_path)) {
    stop(
      "Compiled shared object not found: ", so_path, "\n",
      "The shared object must be built once before submitting the array. ",
      "src/*.o and src/*.so are gitignored, so a fresh pull never has one.\n",
      "Fix: from the repo root, run:\n",
      "  Rscript -e 'pkgbuild::compile_dll()'\n",
      "Migration step 12 adds this as an automated pre-array step in ",
      "sim/chpc/run_sim.sh."
    )
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# load_bhcrr
# ---------------------------------------------------------------------------

load_bhcrr <- function(mode = Sys.getenv("BHCRR_LOAD", "dev"),
                       force = FALSE) {

  if (!requireNamespace("pkgload", quietly = TRUE))
    stop(
      "pkgload is required but not installed. ",
      "Install it with: install.packages(\"pkgload\")"
    )

  # Validate mode before the idempotence check so an invalid mode always errors.
  if (!mode %in% c("dev", "cluster"))
    stop("load_bhcrr(): 'mode' must be \"dev\" or \"cluster\", got \"", mode,
         "\". Set BHCRR_LOAD or pass mode explicitly.")

  # Idempotence: skip reload if already attached (unless forced).
  if (!force && "package:bhCRR" %in% search())
    return(invisible(FALSE))

  root <- repo_root()

  compile <- switch(mode,
    dev     = TRUE,
    cluster = { check_compiled(root); FALSE }
  )

  pkgload::load_all(
    root,
    compile    = compile,
    quiet      = TRUE,
    export_all = FALSE   # match library() behaviour; see header NOTE
  )

  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# load_sim
# ---------------------------------------------------------------------------

load_sim <- function(mode = Sys.getenv("BHCRR_LOAD", "dev"),
                     force = FALSE) {

  load_bhcrr(mode, force)

  root <- repo_root()

  # Explicit list in dependency order. Do NOT replace with dir() globbing:
  # order matters and globbing would silently pick up scratch files.
  # Files that are currently comment-only stubs source as harmless no-ops;
  # this lets drivers be written before all modules are implemented.
  modules <- file.path(root, "sim", "R", c(
    "scenario.r",
    "dgp.r",
    "fit.r",
    "comparators.r",
    "score.r",
    "run_one.r",
    "parse.r"
  ))

  missing <- modules[!file.exists(modules)]
  if (length(missing))
    stop("load_sim(): missing module file(s):\n",
         paste(" ", missing, collapse = "\n"))

  for (f in modules)
    source(f, local = FALSE)  # local = FALSE -> global environment

  invisible(modules)
}
