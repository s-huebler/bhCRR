# sim/R/run_one.r -- scenario + run index -> one complete result list
#
# OWNS
#   run_one(scenario, p, run_id,
#           load_mode = Sys.getenv("BHCRR_LOAD", "dev"),
#           seed = TRUE,
#           fit_control = NULL, s0_seq = NULL, s1_seq = NULL,
#           max_widen = 1L)
#                       the single unit of work; draws data, fits bhCRR,
#                       fits comparators, scores all models, returns a result
#                       list; sim/local/batch.R and sim/chpc/R/run_array.r
#                       both call this function and contain no workflow logic
#                       of their own -- all simulation logic lives here
#
# INPUTS
#   scenario    list or character(1)   validated scenario spec (list of class
#                                      "bhcrr_scenario") or path to a
#                                      scenarios/*.R file (character); when a
#                                      path is given, scenario_from_file() is
#                                      called internally; an unvalidated plain
#                                      list is an error -- pass it through
#                                      scenario() first
#   p           integer(1)             number of predictors for this run
#   run_id      integer(1)             run index; used for seeding and for the
#                                      $run_id field of the result
#   load_mode   character(1)           passed to load_bhcrr(); "dev" or "cluster"
#   seed        logical(1) or integer(1)
#                                      TRUE (default): derive a deterministic
#                                      integer from (scenario$name, p, run_id)
#                                      and call set.seed() once before sim_data().
#                                      FALSE: do not call set.seed() at all,
#                                      reproducing legacy behaviour.
#                                      Any other integer: use that value directly.
#   fit_control bhcrr_cv_control or NULL
#                                      passed to fit_bhcrr(); NULL builds the
#                                      default; MUST NOT carry a non-NULL $seed
#                                      field when seed = TRUE (see B5 guard)
#   s0_seq      numeric or NULL        override the s0 grid in fit_bhcrr()
#   s1_seq      numeric or NULL        override the s1 grid in fit_bhcrr()
#   max_widen   integer(1)             max grid-widening attempts (default 1L)
#
# RETURNS
#   list with fields:
#     $run_id       integer(1)
#     $p            integer(1)
#     $scenario     list          the validated bhcrr_scenario spec used
#
#     $data         list(x, y, beta1, beta2)
#                     x:     matrix n x p  design matrix, colnames X_1..X_p
#                     y:     matrix n x 2  (time, status), status in {0,1,2}
#                     beta1: numeric p     full cause-1 coefficient vector;
#                            active values at their named positions, zeros
#                            elsewhere
#                     beta2: numeric p     full cause-2 coefficient vector
#                   block_map, block_signature and design_settings from
#                   sim_data() are stored in $meta, not $data; one canonical
#                   home avoids drift between $data and $meta
#
#     $fit          list     output of fit_bhcrr() -- see fit.r for contract
#     $comparators  list     output of fit_comparators() -- see comparators.r;
#                            retains top-level setup_seconds and standardization
#                            fields; do NOT merge this with $fg_oracle (see C4)
#     $fg_oracle    list     output of fit_fg_oracle() -- see comparators.r;
#                            stored separately so setup_seconds and
#                            standardization are not mistaken for model outputs
#                            by any downstream code that iterates over models
#
#     $scores       list(overall = data.frame, by_block = data.frame)
#                   output of score_run(); see score.r for column contracts
#
#     $meta         list
#       $meta$run_id           integer(1)
#       $meta$p                integer(1)
#       $meta$scenario_name    character(1)
#       $meta$seed             integer(1) or NA_integer_
#                                seed passed to set.seed(); NA when seed = FALSE
#       $meta$git_sha          character(1)  git SHA of HEAD at run time;
#                                "unknown" on any git failure
#       $meta$bhcrr_version    character(1)  utils::packageVersion("bhCRR")
#       $meta$r_version        character(1)  R.version.string
#       $meta$elapsed_seconds  numeric(1)    wall seconds for the entire
#                                run_one() call; numeric, never difftime
#       $meta$block_map        data.frame  one row per predictor:
#                                name       character  "X_1", "X_2", ...
#                                block      character  block label
#                                block_rho  numeric    within-block rho
#                                designed   logical    in designed active set
#                                active1    logical    beta1 != 0
#                                active2    logical    beta2 != 0
#       $meta$block_signature  character(1)  fingerprint of the realized
#                                column-to-block map (see dgp.r SIGNATURE note)
#       $meta$design_settings  list  the generator's campaign-fixed knobs plus
#                                the realized event mix for this run:
#                                  nobs, npredictors, block_props, latent,
#                                  latent_q, latent_a, latent_balanced,
#                                  unit_variance, layout, event_mix
#                                These are NOT among the 12 scenario signature
#                                keys, so this list is their only record in the
#                                stored result. Changing any of them requires
#                                a code edit (not a new scenario name), and that
#                                change would be silent in the manifest. Storing
#                                them here makes it detectable in old vs. new
#                                run files.
#
# BLOCK MAP PROVENANCE
#   The map uses NO RNG. It is a pure function of p, block_props, block_rho,
#   active_block and layout. The RNG sets column values, never block assignment.
#   Stored rather than recomputed so that a later edit to block_props or layout
#   cannot silently relabel results already on disk.
#
# SEEDING POLICY
#   Legacy (sim/_legacy/code/batch_run.R) explicitly deleted .Random.seed from
#   the global environment between runs to prevent a leaked seed from making
#   consecutive runs correlated. The hazard: if run 2 is always seeded from run
#   1's state, runs across (p, scenario) combinations share a statistical
#   relationship that biases variance estimates.
#
#   This implementation seeds on all three of (scenario$name, p, run_id) using
#   a pure-base-R hash. A seed derived from run_id alone would reintroduce the
#   correlation legacy feared: run 1 at p=25 and run 1 at p=300 would get the
#   same stream. The triple-key seed avoids that while making any single run
#   reproducible from $meta$seed without re-deriving the hash.
#
# ERROR POLICY
#   run_one() does NOT catch errors. It lets them propagate so the driver
#   (sim/local/batch.R for local runs; sim/chpc/R/run_array.r for CHPC array
#   jobs) owns all failure handling: keeping the sweep alive, writing failure
#   records, and logging. A single caller owns recovery; tryCatch here would
#   shadow errors and complicate debugging.
#
# DEPENDS ON
#   sim/R/load.r         repo_root(), load_bhcrr()
#   sim/R/scenario.r     scenario_from_file(), scenario()
#   sim/R/dgp.r          sim_data()
#   sim/R/fit.r          fit_bhcrr()
#   sim/R/comparators.r  fit_comparators(), fit_fg_oracle()
#   sim/R/score.r        score_run()
#
# NOTE ON load_bhcrr() vs load_sim()
#   run_one() calls load_bhcrr(load_mode), NOT load_sim(). load_sim() sources
#   run_one.r itself into the global environment; calling it from inside
#   run_one() would re-source this file during its own execution, which would
#   shadow the running function with a new definition mid-call. load_bhcrr() is
#   idempotent: when a driver has already called load_sim(), this is a no-op.
#
# REPLACES
#   sim/_legacy/code/single_run_complex.R     (entire file; one-run case)
#   sim/_legacy/code/batch_run_complex.R      (the inner loop body)
#   sim/chpc/R/run_batch.R becomes sim/chpc/R/run_array.r, a thin wrapper
#     that reads env vars, resolves the scenario, and calls run_one() with
#     no workflow logic of its own; that file is not absorbed here
#   Manifest write and check are driver-level concerns: write_manifest() is
#     called by sim/local/batch.R (step 10) and the check by
#     sim/chpc/R/run_array.r (step 12). run_one() receives no output directory
#     and therefore cannot write or check a manifest.
#
# WHAT GETS SAVED
#   For now the result keeps as much as possible -- fitted objects, the
#   autotune object, comparator fits, x and y -- so that any diagnostic can
#   be run after the fact while the method is still being troubleshooted.
#   This is a known cost in file size and parse time, accepted on purpose.
#   The revisit trigger is the point where the presented summary statistics
#   are settled. Do NOT trim anything during the current migration pass.

# ---------------------------------------------------------------------------
# .derive_seed -- pure-base-R hash of (scenario_name, p, run_id)
# ---------------------------------------------------------------------------

.derive_seed <- function(scenario_name, p, run_id) {
  # Construct a canonical string from all three identifiers so that no two
  # (scenario, p, run_id) triples share a seed. A separator that cannot appear
  # in a valid scenario name (which must match ^[A-Za-z0-9][A-Za-z0-9.-]*$)
  # is used to prevent ("A", 12, 3) and ("A1", 2, 3) from colliding.
  s <- paste0(scenario_name, "|p=", as.integer(p), "|run=", as.integer(run_id))
  # djb2-style hash over UTF-8 bytes. Use floating-point doubles throughout
  # and reduce modulo 2^31 - 1 (Mersenne prime) at every step to avoid R
  # integer overflow (R's L integers silently produce NA on overflow, which
  # would crash set.seed()). Result is always in [0, 2^31 - 2], a valid
  # set.seed() argument. No package dependencies.
  bytes <- utf8ToInt(s)
  m <- 2147483647.0   # 2^31 - 1
  h <- 5381.0
  for (b in bytes)
    h <- (h * 33.0 + b) %% m
  as.integer(h)
}

# ---------------------------------------------------------------------------
# run_one
# ---------------------------------------------------------------------------

run_one <- function(scenario,
                    p,
                    run_id,
                    load_mode   = Sys.getenv("BHCRR_LOAD", "dev"),
                    seed        = TRUE,
                    fit_control = NULL,
                    s0_seq      = NULL,
                    s1_seq      = NULL,
                    max_widen   = 1L) {

  t_start <- proc.time()[["elapsed"]]

  # ---- C2: load_bhcrr, NOT load_sim (see NOTE ON load_bhcrr vs load_sim) ----
  load_bhcrr(load_mode)

  # ---- C1: resolve scenario ----
  if (is.character(scenario)) {
    scenario <- scenario_from_file(scenario)
  } else if (is.list(scenario)) {
    if (!inherits(scenario, "bhcrr_scenario"))
      stop("run_one(): 'scenario' is a plain list but not of class ",
           "\"bhcrr_scenario\". Pass it through scenario() for validation ",
           "before calling run_one(), or supply a file path instead.")
  } else {
    stop("run_one(): 'scenario' must be a bhcrr_scenario list or a character ",
         "path to a scenarios/*.R file. Got ", class(scenario)[1L], ".")
  }

  # ---- B5: guard against fit_control$seed overriding our seeding ----
  if (isTRUE(seed) && !is.null(fit_control) && !is.null(fit_control$seed))
    stop("run_one(): seed = TRUE but fit_control$seed is non-NULL ",
         "(fit_control$seed = ", fit_control$seed, "). ",
         "bhcrr_cv_control$seed would reset the RNG stream partway through ",
         "the run and break the one-seed-per-run reproducibility property. ",
         "Either set seed = FALSE, remove fit_control$seed, or pass an integer ",
         "seed directly to run_one(seed = <integer>).")

  # ---- B2 / B3: seeding ----
  if (isFALSE(seed)) {
    meta_seed <- NA_integer_
  } else if (isTRUE(seed)) {
    meta_seed <- .derive_seed(scenario$name, p, run_id)
    set.seed(meta_seed)
  } else {
    meta_seed <- as.integer(seed)
    set.seed(meta_seed)
  }

  # ---- C3: pipeline ----
  data   <- sim_data(scenario, p)

  fit    <- fit_bhcrr(data$x, data$y, scenario, run_id,
                      control   = fit_control,
                      s0_seq    = s0_seq,
                      s1_seq    = s1_seq,
                      max_widen = max_widen)

  comps  <- fit_comparators(data$x, data$y)

  oracle <- fit_fg_oracle(data$x, data$y, data$beta1)

  # C3: score with the established merged-comparator pattern (see score.r C3 note)
  scores <- score_run(fit,
                      c(comps, list(fg_oracle = oracle)),
                      data$x, data$y, scenario,
                      data$beta1, data$block_map)

  # ---- C5: meta ----
  git_sha <- tryCatch({
    root <- repo_root()
    trimws(system2("git", c("-C", root, "rev-parse", "--short", "HEAD"),
                   stdout = TRUE, stderr = FALSE))
  }, error   = function(e) "unknown",
     warning = function(e) "unknown")
  if (!nzchar(git_sha)) git_sha <- "unknown"

  elapsed <- as.numeric(proc.time()[["elapsed"]] - t_start)

  meta <- list(
    run_id          = as.integer(run_id),
    p               = as.integer(p),
    scenario_name   = scenario$name,
    seed            = meta_seed,
    git_sha         = git_sha,
    bhcrr_version   = as.character(utils::packageVersion("bhCRR")),
    r_version       = R.version.string,
    elapsed_seconds = elapsed,
    block_map       = data$block_map,
    block_signature = data$block_signature,
    design_settings = data$design_settings
  )

  # ---- C6: $data holds x, y, beta1, beta2 only ----
  list(
    run_id      = as.integer(run_id),
    p           = as.integer(p),
    scenario    = scenario,
    data        = list(x     = data$x,
                       y     = data$y,
                       beta1 = data$beta1,
                       beta2 = data$beta2),
    fit         = fit,
    comparators = comps,
    fg_oracle   = oracle,
    scores      = scores,
    meta        = meta
  )
}
