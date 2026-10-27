# sim/R/scenario.r -- scenario spec: construct, validate, and persist
#
# OWNS
#   scenario(name, nobs, zero_gap_target, beta1_active, beta2_active,
#            active_block, block_props, block_rho, latent_type, latent_q,
#            cens_rate, u_min, u_max)
#                     construct and validate a spec from the 12 DGP keys plus
#                     name; return a list of class "bhcrr_scenario" with
#                     $blocks appended as derived per-block structure
#
#   scenario_from_file(path)
#                     source a scenarios/*.R file into an isolated environment,
#                     verify it assigned exactly `scenario` and nothing else,
#                     enforce the name-filename contract, and validate via
#                     scenario(); refuse any file whose basename begins with "_"
#
#   new_scenario(name, dir = "sim/scenarios")
#                     copy sim/scenarios/_template.R to dir/<name>.R with the
#                     name field filled in; validate the slug; refuse to
#                     overwrite; return the new path invisibly
#
#   write_manifest(spec_or_path, dir)
#                     write scenario.env to dir; spec_or_path can be an
#                     already-constructed bhcrr_scenario or a file path to
#                     load via scenario_from_file(); return invisible(path)
#
#   read_manifest(path)
#                     parse scenario.env back into a named list of raw strings,
#                     one entry per KEY=value data line; splits on the FIRST
#                     "=" only (mirrors bash's cut -d= -f2-), so named-vector
#                     values like "indep=0.5,weak=0.25" are preserved intact
#
#   check_manifest(spec_or_path, dir)
#                     compare spec against the manifest in dir; comparison is
#                     NUMERIC for numeric keys (not string) so "0.40" and "0.4"
#                     do not cause spurious mismatches; named-vector keys are
#                     compared as sets by name, not by written order; stop with
#                     a message naming the key, recorded value, and current
#                     value on any mismatch; return invisible(TRUE) on pass
#
# INPUTS
#   name, nobs, ...:  see the scenario() signature below and _template.R
#   path:             character(1) path to a scenarios/*.R file
#   spec:             a bhcrr_scenario object (output of scenario())
#   spec_or_path:     character(1) path OR bhcrr_scenario; both accepted so
#                     write_manifest and check_manifest are callable from a
#                     one-line Rscript -e invocation without constructing a
#                     spec object first:
#                       Rscript -e 'source("sim/R/load.r"); load_sim(); \
#                         check_manifest("sim/scenarios/myScenario.R", "sim/results/myScenario")'
#   dir:              character(1) directory where scenario.env lives
#
# RETURNS
#   scenario():           list of class "bhcrr_scenario"; has the 12 DGP keys,
#                         $name, $beta2_active_supplied (the form given by the
#                         caller: "neg_beta1" or a numeric vector), and $blocks
#                         (named list, one entry per block label, holding label,
#                         proportion, rho, latent_type, latent_q -- the
#                         structure dgp.r consumes; predictor assignment is
#                         NOT done here, it depends on p and belongs to dgp.r)
#   scenario_from_file(): same as scenario()
#   new_scenario():       invisible(character(1)) -- path to the new file
#   write_manifest():     invisible(character(1)) -- path to scenario.env
#   read_manifest():      named list of raw strings; keys uppercase; data lines
#                         only (comment header lines are excluded)
#   check_manifest():     invisible(TRUE) on pass; stop() on mismatch
#
# SIGNATURE KEYS
#   The 12 keys that define a scenario (identical to _sig_keys in
#   sim/chpc/run_sim.sh) in the exact order that file uses:
#     ACTIVE_BLOCK  BETA1_ACTIVE  BETA2_ACTIVE  BLOCK_PROPS  BLOCK_RHO
#     CENS_RATE     LATENT_Q      LATENT_TYPE   NOBS         U_MAX
#     U_MIN         ZERO_GAP_TARGET
#   Changing any of these requires a new scenario name.
#   NOT signature: NPREDICTORS, RUN_START, RUN_END, RUNS_PER_TASK.
#
# DUAL WRITER NOTE (resolved at step 12b)
#   write_manifest() is now the SOLE writer of scenario.env.
#   run_sim.sh calls write_manifest() via a short Rscript before any sbatch,
#   so on CHPC the manifest is written once on the login node before the array
#   starts. batch_sweep() on CHPC array tasks passes manifest = "check_only"
#   and never writes. The old bash writer (_sig_keys array + string comparison
#   loop) was deleted; it compared values as raw strings, causing spurious
#   mismatches between "0.4" and "0.40". check_manifest() uses numeric
#   comparison and is immune to that; it is the sole reader and checker.
#   write_manifest() uses file.rename() for atomicity so a concurrent or
#   interrupted write cannot leave a partial file.
#
# DEPENDS ON
#   sim/R/  (none; scenario.r has no intra-sim dependencies)
#   bhCRR:: (none)
#
# REPLACES
#   The scenario manifest logic in sim/chpc/run_sim.sh (bash) lines 96-123
#   The scenario spec is the single definition of the DGP; the scenarios/
#   directory is its only authoritative home

# ---------------------------------------------------------------------------
# Internal constants
# ---------------------------------------------------------------------------

# 12 signature keys in the exact order from run_sim.sh's _sig_keys array
.sig_keys <- c("ACTIVE_BLOCK", "BETA1_ACTIVE", "BETA2_ACTIVE", "BLOCK_PROPS",
               "BLOCK_RHO", "CENS_RATE", "LATENT_Q", "LATENT_TYPE",
               "NOBS", "U_MAX", "U_MIN", "ZERO_GAP_TARGET")

# 13 user-settable keys (12 DGP keys + name)
.scenario_user_keys <- c("name", "nobs", "zero_gap_target",
                          "beta1_active", "beta2_active", "active_block",
                          "block_props", "block_rho", "latent_type",
                          "latent_q", "cens_rate", "u_min", "u_max")

# Forbidden patterns in scenario files (non-comment lines); checked as text
# before sourcing so a library() in a comment does not trigger the check
.forbidden_patterns <- c("library\\s*\\(", "require\\s*\\(",
                          "source\\s*\\(", "set\\.seed\\s*\\(")

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

.check_slug <- function(name) {
  if (!grepl("^[A-Za-z0-9][A-Za-z0-9.-]*$", name))
    stop(
      "'name' must match ^[A-Za-z0-9][A-Za-z0-9.-]*$ (no underscores).\n",
      "  Underscores are field delimiters in the p<P>_<scenario>_run<N> ",
      "filename format used by run_sim.sh and parse.r; a name containing '_' ",
      "would break filename parsing. Use '.' or '-' as word separators instead.\n",
      "  Got: \"", name, "\""
    )
}

.ser_num_vec <- function(x) paste(x, collapse = ",")

.ser_named_num <- function(x)
  paste(paste0(names(x), "=", x), collapse = ",")

.parse_named_num <- function(s) {
  parts <- strsplit(s, ",")[[1]]
  vals <- numeric(length(parts))
  nms  <- character(length(parts))
  for (i in seq_along(parts)) {
    eq <- regexpr("=", parts[i], fixed = TRUE)[1]
    nms[i]  <- substr(parts[i], 1, eq - 1)
    vals[i] <- as.numeric(substr(parts[i], eq + 1, nchar(parts[i])))
  }
  stats::setNames(vals, nms)
}

.num_eq <- function(a, b, tol = 1e-9) isTRUE(all(abs(a - b) < tol))

.named_num_eq <- function(a, b, tol = 1e-9) {
  if (!setequal(names(a), names(b))) return(FALSE)
  .num_eq(a[names(a)], b[names(a)], tol)
}

.git_sha <- function() {
  tryCatch(
    trimws(system2("git", c("rev-parse", "--short", "HEAD"),
                   stdout = TRUE, stderr = FALSE)[1]),
    error   = function(e) "unknown",
    warning = function(w) "unknown"
  )
}

.load_spec <- function(spec_or_path) {
  if (is.character(spec_or_path)) scenario_from_file(spec_or_path) else spec_or_path
}

# ---------------------------------------------------------------------------
# scenario() -- constructor and validator
# ---------------------------------------------------------------------------

scenario <- function(name,
                     nobs,
                     zero_gap_target,
                     beta1_active,
                     beta2_active,
                     active_block,
                     block_props,
                     block_rho,
                     latent_type = "gaussian",
                     latent_q,
                     cens_rate,
                     u_min,
                     u_max) {

  # ---- name ----------------------------------------------------------------
  .check_slug(name)

  # ---- beta2_active: resolve "neg_beta1" sentinel -------------------------
  beta2_active_supplied <- beta2_active
  if (is.character(beta2_active)) {
    if (!identical(beta2_active, "neg_beta1"))
      stop("'beta2_active' must be a numeric vector or the string \"neg_beta1\".")
    beta2_active <- -beta1_active
  }

  # ---- length consistency --------------------------------------------------
  n_active <- length(beta1_active)
  if (length(active_block) != n_active)
    stop("'active_block' must have the same length as 'beta1_active' (",
         n_active, "); got length ", length(active_block), ".")
  if (is.numeric(beta2_active_supplied) &&
      length(beta2_active_supplied) != n_active)
    stop("'beta2_active' (when a vector) must have the same length as ",
         "'beta1_active' (", n_active, "); got length ",
         length(beta2_active_supplied), ".")

  # ---- block labels --------------------------------------------------------
  label_set <- names(block_props)
  if (is.null(label_set) || length(label_set) == 0)
    stop("'block_props' must be a NAMED numeric vector (block labels as names).")
  if (is.null(names(block_rho)) || !setequal(names(block_rho), label_set))
    stop("'block_rho' must have the same label set as 'block_props'.\n",
         "  block_props labels: ", paste(sort(label_set), collapse = ", "), "\n",
         "  block_rho labels:   ",
         paste(sort(names(block_rho)), collapse = ", "))
  if (is.null(names(latent_q)) || !setequal(names(latent_q), label_set))
    stop("'latent_q' must have the same label set as 'block_props'.\n",
         "  block_props labels: ", paste(sort(label_set), collapse = ", "), "\n",
         "  latent_q labels:    ",
         paste(sort(names(latent_q)), collapse = ", "))
  bad_blocks <- setdiff(active_block, label_set)
  if (length(bad_blocks))
    stop("'active_block' entries not in block_props labels: ",
         paste(bad_blocks, collapse = ", "), ".\n",
         "  Known labels: ", paste(label_set, collapse = ", "))

  # ---- range checks --------------------------------------------------------
  if (!is.numeric(nobs) || length(nobs) != 1 || nobs < 1)
    stop("'nobs' must be a positive integer; got ", nobs, ".")
  tol <- sqrt(.Machine$double.eps)
  if (abs(sum(block_props) - 1) > tol)
    stop("'block_props' must sum to 1; got ", sum(block_props), ".")
  if (any(block_rho < 0) || any(block_rho >= 1))
    stop("'block_rho' values must be in [0, 1); got: ",
         paste(block_rho, collapse = ", "))
  if (any(latent_q <= 0) || any(latent_q >= 1))
    stop("'latent_q' values must be in (0, 1); got: ",
         paste(latent_q, collapse = ", "))
  if (!is.numeric(cens_rate) || length(cens_rate) != 1 || cens_rate < 0)
    stop("'cens_rate' must be a non-negative number; got ", cens_rate, ".")
  if (u_min > u_max)
    stop("'u_min' (", u_min, ") must be <= 'u_max' (", u_max, ").")
  if (!latent_type %in% c("gaussian", "bernoulli"))
    stop("'latent_type' must be \"gaussian\" or \"bernoulli\"; got \"",
         latent_type, "\".")

  # ---- build $blocks: per-label derived structure for dgp.r ---------------
  blocks <- lapply(label_set, function(lbl) {
    list(
      label       = lbl,
      proportion  = unname(block_props[lbl]),
      rho         = unname(block_rho[lbl]),
      latent_type = latent_type,
      latent_q    = unname(latent_q[lbl])
    )
  })
  names(blocks) <- label_set

  structure(
    list(
      name                  = name,
      nobs                  = as.integer(nobs),
      zero_gap_target       = zero_gap_target,
      beta1_active          = beta1_active,
      beta2_active          = beta2_active,
      beta2_active_supplied = beta2_active_supplied,
      active_block          = active_block,
      block_props           = block_props,
      block_rho             = block_rho,
      latent_type           = latent_type,
      latent_q              = latent_q,
      cens_rate             = cens_rate,
      u_min                 = u_min,
      u_max                 = u_max,
      blocks                = blocks
    ),
    class = c("bhcrr_scenario", "list")
  )
}

# ---------------------------------------------------------------------------
# scenario_from_file()
# ---------------------------------------------------------------------------

scenario_from_file <- function(path) {
  path <- normalizePath(path, mustWork = FALSE)
  bn   <- basename(path)

  # Refuse template files (leading underscore)
  if (startsWith(bn, "_"))
    stop("scenario_from_file() refuses files whose basename begins with '_'.\n",
         "  '", bn, "' is a template. ",
         "Copy it with new_scenario() to create a real scenario.")

  if (!file.exists(path))
    stop("Scenario file not found: '", path, "'")

  # Check for forbidden calls before sourcing (text scan, not execution)
  lines        <- readLines(path, warn = FALSE)
  non_comment  <- lines[!grepl("^\\s*#", lines)]
  for (pat in .forbidden_patterns) {
    if (any(grepl(pat, non_comment)))
      stop("Scenario file '", bn, "' contains a forbidden call (pattern: '",
           sub("\\\\s\\*\\\\\\(", "()", pat),
           "').\n  Scenario files must be pure list constructors; no ",
           "library(), require(), source(), or set.seed() calls.")
  }

  # Source into an isolated environment so side effects are contained
  e <- new.env(parent = baseenv())
  tryCatch(
    eval(parse(file = path), envir = e),
    error = function(err)
      stop("Error parsing/sourcing '", bn, "':\n  ", conditionMessage(err))
  )

  # Must assign exactly one object named "scenario"
  objs <- ls(envir = e, all.names = TRUE)
  if (!identical(objs, "scenario"))
    stop("Scenario file '", bn,
         "' must assign exactly one object named 'scenario' and nothing else.\n",
         "  Found: ", paste(sort(objs), collapse = ", "))

  spec_raw <- get("scenario", envir = e)
  if (!is.list(spec_raw))
    stop("'scenario' in '", bn, "' must be a list; got ",
         class(spec_raw)[1], ".")
  if (is.null(spec_raw$name))
    stop("'scenario$name' is missing in '", bn, "'.")

  # name field must match filename (without extension)
  name_base <- sub("\\.[Rr]$", "", bn)
  if (!identical(spec_raw$name, name_base))
    stop("scenario$name (\"", spec_raw$name, "\") does not match the filename ",
         "(\"", bn, "\").\n",
         "  A mismatch would silently write results into the wrong directory.")

  # Filter to known user keys (ignore 'blocks' if erroneously included)
  user_keys <- intersect(.scenario_user_keys, names(spec_raw))
  do.call(scenario, spec_raw[user_keys])
}

# ---------------------------------------------------------------------------
# new_scenario()
# ---------------------------------------------------------------------------

new_scenario <- function(name, dir = "sim/scenarios") {
  .check_slug(name)

  out_path <- file.path(dir, paste0(name, ".R"))
  if (file.exists(out_path))
    stop("File already exists: '", out_path, "'.\n",
         "  Choose a different name or delete the existing file first.")

  template_path <- file.path(repo_root(), "sim", "scenarios", "_template.R")
  if (!file.exists(template_path))
    stop("Template not found at: '", template_path, "'")

  lines <- readLines(template_path, warn = FALSE)
  # Replace the name placeholder with the new name
  lines <- sub('"FILL\\.IN\\.NAME"', sprintf('"%s"', name), lines)

  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  writeLines(lines, out_path)
  invisible(out_path)
}

# ---------------------------------------------------------------------------
# write_manifest()
# ---------------------------------------------------------------------------

write_manifest <- function(spec_or_path, dir) {
  spec <- .load_spec(spec_or_path)

  # BETA2_ACTIVE: write the SUPPLIED form to match what bash would write
  b2_str <- if (identical(spec$beta2_active_supplied, "neg_beta1")) {
    "neg_beta1"
  } else {
    .ser_num_vec(spec$beta2_active_supplied)
  }

  # Keys in the exact _sig_keys order from run_sim.sh
  lines <- c(
    sprintf("# scenario:  %s", spec$name),
    sprintf("# created:   %s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S")),
    sprintf("# git:       %s", .git_sha()),
    sprintf("ACTIVE_BLOCK=%s",    .ser_num_vec(spec$active_block)),
    sprintf("BETA1_ACTIVE=%s",    .ser_num_vec(spec$beta1_active)),
    sprintf("BETA2_ACTIVE=%s",    b2_str),
    sprintf("BLOCK_PROPS=%s",     .ser_named_num(spec$block_props)),
    sprintf("BLOCK_RHO=%s",       .ser_named_num(spec$block_rho)),
    sprintf("CENS_RATE=%s",       spec$cens_rate),
    sprintf("LATENT_Q=%s",        .ser_named_num(spec$latent_q)),
    sprintf("LATENT_TYPE=%s",     spec$latent_type),
    sprintf("NOBS=%s",            spec$nobs),
    sprintf("U_MAX=%s",           spec$u_max),
    sprintf("U_MIN=%s",           spec$u_min),
    sprintf("ZERO_GAP_TARGET=%s", spec$zero_gap_target)
  )

  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  path <- file.path(dir, "scenario.env")
  tmp  <- tempfile(tmpdir = dir, fileext = ".env.tmp")
  writeLines(lines, tmp)
  file.rename(tmp, path)
  invisible(path)
}

# ---------------------------------------------------------------------------
# read_manifest()
# ---------------------------------------------------------------------------

read_manifest <- function(path) {
  if (!file.exists(path))
    stop("Manifest not found: '", path, "'")
  lines <- readLines(path, warn = FALSE)
  # Data lines are uppercase-key lines; comment header lines start with "#"
  data_lines <- lines[grepl("^[A-Z]", lines)]
  result <- vector("list", length(data_lines))
  keys   <- character(length(data_lines))
  for (i in seq_along(data_lines)) {
    eq     <- regexpr("=", data_lines[i], fixed = TRUE)[1]
    keys[i]     <- substr(data_lines[i], 1, eq - 1)
    result[[i]] <- substr(data_lines[i], eq + 1, nchar(data_lines[i]))
  }
  stats::setNames(result, keys)
}

# ---------------------------------------------------------------------------
# check_manifest()
# ---------------------------------------------------------------------------

check_manifest <- function(spec_or_path, dir) {
  spec <- .load_spec(spec_or_path)

  env_path <- file.path(dir, "scenario.env")
  if (!file.exists(env_path))
    stop("No scenario.env found at: '", env_path, "'")

  rec <- read_manifest(env_path)
  tol <- 1e-9

  # Helper: report mismatch for a key
  .mismatch <- function(key, current, recorded_str) {
    stop(sprintf(
      "Scenario mismatch -- %s:\n  recorded : %s\n  current  : %s\n%s",
      key, recorded_str, paste(current, collapse = ","),
      "Fix the mismatch or choose a new scenario name."
    ))
  }

  # ACTIVE_BLOCK -- character vector comparison
  if (!is.null(rec$ACTIVE_BLOCK)) {
    rec_ab <- strsplit(rec$ACTIVE_BLOCK, ",")[[1]]
    if (!identical(rec_ab, spec$active_block))
      .mismatch("ACTIVE_BLOCK", spec$active_block, rec$ACTIVE_BLOCK)
  }

  # BETA1_ACTIVE -- numeric vector comparison (immune to "0.40" vs "0.4")
  if (!is.null(rec$BETA1_ACTIVE)) {
    rec_b1 <- as.numeric(strsplit(rec$BETA1_ACTIVE, ",")[[1]])
    if (!.num_eq(rec_b1, spec$beta1_active, tol))
      .mismatch("BETA1_ACTIVE", spec$beta1_active, rec$BETA1_ACTIVE)
  }

  # BETA2_ACTIVE -- resolve "neg_beta1" before comparing numerically
  if (!is.null(rec$BETA2_ACTIVE)) {
    rec_b2 <- if (identical(rec$BETA2_ACTIVE, "neg_beta1")) {
      -spec$beta1_active
    } else {
      as.numeric(strsplit(rec$BETA2_ACTIVE, ",")[[1]])
    }
    if (!.num_eq(rec_b2, spec$beta2_active, tol))
      .mismatch("BETA2_ACTIVE", spec$beta2_active, rec$BETA2_ACTIVE)
  }

  # BLOCK_PROPS, BLOCK_RHO, LATENT_Q -- named numeric sets by name
  for (key in c("BLOCK_PROPS", "BLOCK_RHO", "LATENT_Q")) {
    if (!is.null(rec[[key]])) {
      rec_val  <- .parse_named_num(rec[[key]])
      spec_val <- switch(key,
        BLOCK_PROPS = spec$block_props,
        BLOCK_RHO   = spec$block_rho,
        LATENT_Q    = spec$latent_q
      )
      if (!.named_num_eq(rec_val, spec_val, tol))
        .mismatch(key, spec_val, rec[[key]])
    }
  }

  # CENS_RATE, NOBS, U_MAX, U_MIN, ZERO_GAP_TARGET -- single numeric
  for (key in c("CENS_RATE", "NOBS", "U_MAX", "U_MIN", "ZERO_GAP_TARGET")) {
    if (!is.null(rec[[key]])) {
      rec_val  <- as.numeric(rec[[key]])
      spec_val <- switch(key,
        CENS_RATE       = spec$cens_rate,
        NOBS            = as.numeric(spec$nobs),
        U_MAX           = spec$u_max,
        U_MIN           = spec$u_min,
        ZERO_GAP_TARGET = spec$zero_gap_target
      )
      if (!.num_eq(rec_val, spec_val, tol))
        .mismatch(key, spec_val, rec[[key]])
    }
  }

  # LATENT_TYPE -- string comparison
  if (!is.null(rec$LATENT_TYPE)) {
    if (!identical(rec$LATENT_TYPE, spec$latent_type))
      .mismatch("LATENT_TYPE", spec$latent_type, rec$LATENT_TYPE)
  }

  invisible(TRUE)
}
