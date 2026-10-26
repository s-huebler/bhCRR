# sim/R/dgp.r -- block-correlated covariates and competing-risks outcomes
#
# OWNS
#   .generate_block_design(nobs, npredictors, beta1_active, beta2_active,
#                          active_block, block_props, block_rho, latent,
#                          latent_q, latent_balanced, unit_variance, layout,
#                          block_id, predictor_prefix)
#                    internal; full port of generate_block_X.R; returns a
#                    named list with X, block_id, block_size, block_rho,
#                    signature, beta1, beta2, active_idx, meta, settings
#
#   sim_covariates(n, p, blocks, latent_balanced = TRUE, unit_variance = TRUE,
#                  layout = "active_first", block_id = NULL,
#                  predictor_prefix = "X_")
#                    public wrapper; rebuilds named-vector args from blocks and
#                    calls .generate_block_design with no active predictors;
#                    returns the design matrix only (column names X_1..X_p)
#
#   check_design_unchanged(supplied, returned)
#                    internal; stop() if the matrix returned by
#                    simulateTwoCauseFineGrayModel differs from the supplied
#                    one beyond 1e-12 tolerance; extracted from the call flow
#                    so it is independently unit-testable; note: this catches
#                    substitution but not internal rescaling of X before event
#                    generation (that rescaling, if present, would be silent)
#
#   sim_outcomes(x, scenario)
#                    draw (time, status); assumes active_first layout, so the
#                    first length(scenario$beta1_active) columns of x are the
#                    active predictors; reconstructs full-length beta vectors
#                    from that assumption; returns list(y, event_mix)
#
#   sim_data(scenario, p)
#                    convenience wrapper; calls .generate_block_design then
#                    sim_outcomes, assembles all contract fields; returns
#                    list(x, y, beta1, beta2, block_map, block_signature,
#                         design_settings)
#
# INPUTS
#   n         integer(1)   number of observations
#   p         integer(1)   total number of predictors
#   blocks    list         scenario$blocks; one element per block label, each
#                          holding: label, proportion, rho, latent_type,
#                          latent_q; NO per-predictor index information
#                          (predictor-to-block assignment depends on p and is
#                          computed here from proportions -- see A2 fix)
#   x         matrix n x p column names X_1..X_p; from sim_covariates()
#   scenario  list         validated bhcrr_scenario from scenario.r
#   supplied  matrix       design matrix given to simulateTwoCauseFineGrayModel
#   returned  matrix       $X returned by simulateTwoCauseFineGrayModel
#
# RETURNS
#   sim_covariates:   numeric matrix n x p, column names X_1..X_p
#                     (A1 fix: legacy uses predictor_prefix = "X_", so columns
#                     are X_1..X_p, not V1..Vp; downstream code depends on it:
#                     comparators build paste0("X_", seq_len(p)) and parse.r
#                     joins the block map by predictor name)
#
#   check_design_unchanged: invisible(NULL); stop() on mismatch
#
#   sim_outcomes:     list(
#                       y         matrix n x 2, colnames c("time","status"),
#                                 status integer in {0L,1L,2L}
#                       event_mix named numeric(3): censored, cause1, cause2
#                                 proportions; run_one records these in $meta
#                     )
#
#   sim_data:         list(
#                       x               matrix n x p
#                       y               matrix n x 2
#                       beta1           numeric p  full cause-1 coef vector;
#                                       active values in columns 1..k under
#                                       active_first layout, zeros elsewhere
#                       beta2           numeric p  full cause-2 coef vector
#                       block_map       data.frame, one row per predictor:
#                                         name      character  "X_1", ...
#                                         block     character  block label
#                                         block_rho numeric    within-block rho
#                                                   (A3 fix: column is
#                                                   block_rho, NOT rho;
#                                                   parse.r must read this name)
#                                         designed  logical
#                                         active1   logical  beta1 != 0
#                                         active2   logical  beta2 != 0
#                       block_signature character(1) fingerprint of the column-
#                                       to-block map; see SIGNATURE note below
#                       design_settings list  all generator knobs plus event
#                                       mix, for provenance (see CAMPAIGN-FIXED
#                                       note); run_one.r (step 9) must carry
#                                       this into $meta
#                     )
#
# SIGNATURE (extended from generate_block_X.R -- PART C change)
#   Legacy signature omits latent, latent_balanced and unit_variance, all of
#   which change the generated distribution. Since the signature's purpose is
#   to identify runs that may be pooled, that is a hole. Extended format:
#     p=<p>|blocks=<names>|size=<sizes>|rho=<rhos>|layout=<layout>|
#     active=<labels>|pinned=<bool>|latent=<latent>|lbal=<T/F>|uvar=<T/F>
#   The first seven fields are identical to generate_block_X.R and unchanged.
#   Three new fields follow "|pinned=...". This is a deliberate divergence;
#   all existing run files are being deleted so there is no compatibility cost.
#
# CAMPAIGN-FIXED DEFAULTS
#   latent_balanced, unit_variance, layout, block_id and predictor_prefix are
#   NOT among the 12 scenario signature keys and cannot be set from a scenario
#   file. sim_data() always passes the campaign defaults below; changing one
#   requires a code edit. All three are embedded in block_signature (see above)
#   so the change is detectable in stored results.
#     latent_balanced  = TRUE
#     unit_variance    = TRUE
#     layout           = "active_first"
#     predictor_prefix = "X_"
#
# NOTE ON block_corr_summary()
#   generate_block_X.R (lines 419-460) includes block_corr_summary(), a
#   validation diagnostic that computes realized within- vs cross-block
#   correlations. It is NOT ported here. It belongs in sim/R/diagnostics.r
#   (migration step TBD) since it is a development tool, not part of the DGP.
#
# DEPENDS ON
#   sim/R/scenario.r                            (bhcrr_scenario fields)
#   fastcmprsk::simulateTwoCauseFineGrayModel   (outcome generation)
#   stats::rexp                                 (exponential censoring)
#   stats::rnorm                                (covariate draws)
#
# REPLACES
#   sim/_legacy/code/generate_block_X.R  (entire file)
#   sim/_legacy/code/batch_run_complex.R lines 104-160 (outcome generation)

# ---------------------------------------------------------------------------
# .generate_block_design -- internal full port of generate_block_X.R
# ---------------------------------------------------------------------------

.generate_block_design <- function(
    nobs, npredictors,
    beta1_active     = numeric(0),
    beta2_active     = numeric(0),
    active_block     = character(0),
    block_props,
    block_rho,
    latent           = c("gaussian", "bernoulli"),
    latent_q,
    latent_balanced  = TRUE,
    unit_variance    = TRUE,
    layout           = c("active_first", "contiguous"),
    block_id         = NULL,
    predictor_prefix = "X_") {

  latent <- match.arg(latent)
  layout <- match.arg(layout)
  block_id_in <- block_id
  block_id <- NULL

  ## 1. Scalar validation
  if (length(nobs) != 1L || !is.finite(nobs) || nobs < 2L)
    stop("nobs must be a single finite integer of at least 2.")
  if (length(npredictors) != 1L || !is.finite(npredictors) || npredictors < 1L)
    stop("npredictors must be a single finite positive integer.")
  nobs <- as.integer(nobs)
  npredictors <- as.integer(npredictors)

  ## 2. Block specification
  block_names <- names(block_props)
  if (is.null(block_names) || any(!nzchar(block_names)) ||
      anyDuplicated(block_names) > 0L)
    stop("block_props must be a named vector with unique, non-empty block names.")
  if (any(!is.finite(block_props)) || any(block_props < 0) ||
      sum(block_props) <= 0)
    stop("block_props must be finite, non-negative, and not all zero.")
  if (!setequal(names(block_rho), block_names))
    stop("block_rho must be named with exactly the same blocks as block_props: ",
         paste(block_names, collapse = ", "), ".")
  block_rho <- block_rho[block_names]
  if (any(!is.finite(block_rho)) || any(block_rho < 0) || any(block_rho >= 1))
    stop("block_rho values must be finite and in [0, 1).")
  if (latent == "bernoulli") {
    if (!setequal(names(latent_q), block_names))
      stop("latent_q must be named with exactly the same blocks as block_props.")
    latent_q <- latent_q[block_names]
    if (any(!is.finite(latent_q)) || any(latent_q <= 0) || any(latent_q >= 1))
      stop("latent_q values must be finite and strictly between 0 and 1.")
  }

  ## 3. Active-predictor specification
  n_active <- length(beta1_active)
  if (length(beta2_active) != n_active || length(active_block) != n_active)
    stop("beta1_active, beta2_active and active_block must all have the same ",
         "length (got ", n_active, ", ", length(beta2_active), ", ",
         length(active_block), ").")
  if (n_active > npredictors)
    stop("More active predictors (", n_active, ") than predictors (",
         npredictors, ").")
  if (any(!is.finite(c(beta1_active, beta2_active))))
    stop("beta1_active and beta2_active must be finite.")
  active_block <- as.character(active_block)
  unknown_block <- setdiff(active_block, block_names)
  if (length(unknown_block))
    stop("active_block names blocks that do not exist: ",
         paste(unique(unknown_block), collapse = ", "), ".")

  ## 4. Block sizes (Hamilton apportionment; deterministic, no RNG)
  min_size <- ifelse(block_rho > 0, 2L, 1L)
  names(min_size) <- block_names
  if (npredictors < sum(min_size))
    stop("npredictors = ", npredictors, " is too small for this block ",
         "specification; needs at least ", sum(min_size), ".")

  if (is.null(block_id_in)) {
    target     <- block_props / sum(block_props) * npredictors
    block_size <- floor(target)
    shortfall  <- npredictors - sum(block_size)
    if (shortfall > 0L) {
      bump <- order(target - floor(target), decreasing = TRUE)[seq_len(shortfall)]
      block_size[bump] <- block_size[bump] + 1L
    }
    names(block_size) <- block_names
    adjusted <- FALSE
    repeat {
      short <- which(block_size < min_size)
      if (!length(short)) break
      donor <- which.max(block_size - min_size)
      if (block_size[donor] - min_size[donor] <= 0L)
        stop("Cannot satisfy block minimum sizes at npredictors = ", npredictors,
             "; raise npredictors or change block_props.")
      block_size[short[1L]] <- block_size[short[1L]] + 1L
      block_size[donor]     <- block_size[donor] - 1L
      adjusted <- TRUE
    }
    if (adjusted)
      warning("block_props adjusted so every block meets its minimum size. ",
              "Realized sizes: ",
              paste(block_names, block_size, sep = "=", collapse = ", "), ".")
  } else {
    if (length(block_id_in) != npredictors)
      stop("block_id must have one entry per predictor (expected ", npredictors,
           ", got ", length(block_id_in), ").")
    block_id_in <- as.character(block_id_in)
    unknown_pin <- setdiff(block_id_in, block_names)
    if (length(unknown_pin))
      stop("block_id names blocks not in block_props: ",
           paste(unique(unknown_pin), collapse = ", "), ".")
    block_size <- as.integer(table(factor(block_id_in, levels = block_names)))
    names(block_size) <- block_names
    too_small <- which(block_size < min_size)
    if (length(too_small))
      stop("Pinned block_id leaves block(s) ",
           paste(block_names[too_small], collapse = ", "),
           " below minimum size.")
  }

  active_per_block <- table(factor(active_block, levels = block_names))
  overfull <- which(as.integer(active_per_block) > block_size)
  if (length(overfull))
    stop("Block(s) ", paste(block_names[overfull], collapse = ", "),
         " assigned more active predictors than columns.")

  ## 5. Column layout
  block_id   <- character(npredictors)
  active_idx <- integer(n_active)

  if (!is.null(block_id_in)) {
    block_id <- block_id_in
    if (layout == "active_first") {
      if (n_active > 0L &&
          !identical(block_id[seq_len(n_active)], active_block))
        stop("Pinned block_id must agree with active_block in its first ",
             n_active, " entries under layout = \"active_first\".")
      active_idx <- seq_len(n_active)
    } else {
      for (k in seq_len(n_active)) {
        b    <- active_block[k]
        free <- setdiff(which(block_id == b), active_idx[seq_len(k - 1L)])
        active_idx[k] <- free[1L]
      }
    }
  } else if (layout == "active_first") {
    if (n_active > 0L) {
      block_id[seq_len(n_active)] <- active_block
      active_idx <- seq_len(n_active)
    }
    filler_counts <- block_size - as.integer(active_per_block)
    if (npredictors > n_active)
      block_id[(n_active + 1L):npredictors] <-
        rep(block_names, times = filler_counts)
  } else {
    block_id   <- rep(block_names, times = block_size)
    next_free  <- cumsum(c(0L, block_size[-length(block_size)])) + 1L
    names(next_free) <- block_names
    for (k in seq_len(n_active)) {
      b              <- active_block[k]
      active_idx[k]  <- next_free[[b]]
      next_free[[b]] <- next_free[[b]] + 1L
    }
  }

  ## 6. Draw covariates (shift-plus-noise scaffold)
  X        <- matrix(0.0, nrow = nobs, ncol = npredictors)
  shift_sd <- stats::setNames(sqrt(block_rho / (1 - block_rho)), block_names)
  latent_a <- stats::setNames(rep(NA_real_, length(block_names)), block_names)

  for (b in block_names) {
    cols <- which(block_id == b)
    if (!length(cols)) next

    v <- shift_sd[[b]]^2
    if (v <= 0) {
      shift        <- rep(0.0, nobs)
      latent_a[[b]] <- 0
    } else if (latent == "gaussian") {
      shift <- shift_sd[[b]] * stats::rnorm(nobs)
      # latent_a stays NA for gaussian; the shift scale IS shift_sd[[b]] but
      # the a-inversion formula is specific to the Bernoulli construction
    } else {
      q             <- latent_q[[b]]
      a             <- sqrt(v / (q * (1 - q)))
      latent_a[[b]] <- a
      ind <- if (latent_balanced) {
        n1 <- round(q * nobs)
        sample(rep(c(1.0, 0.0), times = c(n1, nobs - n1)))
      } else {
        as.numeric(stats::runif(nobs) < q)
      }
      shift <- a * (ind - q)
    }

    eps <- matrix(stats::rnorm(nobs * length(cols)), nrow = nobs)
    Z   <- shift + eps

    ## 7. Unit-variance rescale
    if (unit_variance)
      Z <- Z * sqrt(1 - block_rho[[b]])

    X[, cols] <- Z
  }

  ## 8. Scatter betas
  predictor_names <- paste0(predictor_prefix, seq_len(npredictors))
  colnames(X) <- predictor_names

  beta1 <- numeric(npredictors)
  beta2 <- numeric(npredictors)
  if (n_active > 0L) {
    beta1[active_idx] <- as.numeric(beta1_active)
    beta2[active_idx] <- as.numeric(beta2_active)
  }
  names(beta1) <- predictor_names
  names(beta2) <- predictor_names

  ## 9. Per-predictor metadata
  meta <- data.frame(
    index    = seq_len(npredictors),
    name     = predictor_names,
    block    = block_id,
    rho      = as.numeric(block_rho[block_id]),   # raw name; sim_data renames
    beta1    = unname(beta1),
    beta2    = unname(beta2),
    designed = seq_len(npredictors) %in% active_idx,
    active1  = unname(beta1 != 0),
    active2  = unname(beta2 != 0),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  ## 10. Block signature (PART C extension)
  ## The three new fields at the end (latent, lbal, uvar) extend the legacy
  ## format from generate_block_X.R. They were absent from the original even
  ## though all three change the generated distribution. The first seven pipe-
  ## delimited fields are identical and in the same order as the original.
  signature <- paste0(
    "p=", npredictors,
    "|blocks=", paste(block_names, collapse = ","),
    "|size=", paste(as.integer(block_size), collapse = ","),
    "|rho=", paste(format(as.numeric(block_rho), digits = 6), collapse = ","),
    "|layout=", layout,
    "|active=", paste(active_block, collapse = ","),
    "|pinned=", !is.null(block_id_in),
    "|latent=", latent,
    "|lbal=", latent_balanced,
    "|uvar=", unit_variance
  )

  list(
    X          = X,
    block_id   = block_id,
    block_size = block_size,
    block_rho  = block_rho,
    signature  = signature,
    beta1      = beta1,
    beta2      = beta2,
    active_idx = active_idx,
    meta       = meta,
    settings   = list(
      nobs            = nobs,
      npredictors     = npredictors,
      block_props     = block_props,
      latent          = latent,
      latent_q        = if (latent == "bernoulli") latent_q else NULL,
      latent_a        = latent_a,
      latent_balanced = latent_balanced,
      unit_variance   = unit_variance,
      layout          = layout
    )
  )
}

# ---------------------------------------------------------------------------
# sim_covariates -- public; rebuilds legacy args from blocks list
# ---------------------------------------------------------------------------

sim_covariates <- function(n, p, blocks,
                           latent_balanced  = TRUE,
                           unit_variance    = TRUE,
                           layout           = "active_first",
                           block_id         = NULL,
                           predictor_prefix = "X_") {
  block_names <- names(blocks)
  block_props <- stats::setNames(
    vapply(blocks, `[[`, numeric(1), "proportion"), block_names)
  block_rho_v <- stats::setNames(
    vapply(blocks, `[[`, numeric(1), "rho"),        block_names)
  latent      <- blocks[[1L]]$latent_type   # scenario.r validates single type
  latent_q    <- stats::setNames(
    vapply(blocks, `[[`, numeric(1), "latent_q"),   block_names)

  design <- .generate_block_design(
    nobs             = n,
    npredictors      = p,
    beta1_active     = numeric(0),
    beta2_active     = numeric(0),
    active_block     = character(0),
    block_props      = block_props,
    block_rho        = block_rho_v,
    latent           = latent,
    latent_q         = latent_q,
    latent_balanced  = latent_balanced,
    unit_variance    = unit_variance,
    layout           = layout,
    block_id         = block_id,
    predictor_prefix = predictor_prefix
  )
  design$X
}

# ---------------------------------------------------------------------------
# check_design_unchanged -- internal
# ---------------------------------------------------------------------------

check_design_unchanged <- function(supplied, returned) {
  if (!isTRUE(all.equal(unname(as.matrix(returned)),
                        unname(supplied),
                        tolerance = 1e-12)))
    stop("simulateTwoCauseFineGrayModel() did not return the supplied design ",
         "matrix unchanged; the event times were not generated from X_design.")
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# sim_outcomes
# ---------------------------------------------------------------------------

sim_outcomes <- function(x, scenario) {
  if (!requireNamespace("fastcmprsk", quietly = TRUE))
    stop("fastcmprsk is required for sim_outcomes(); install it first.")

  n     <- nrow(x)
  p     <- ncol(x)
  k     <- length(scenario$beta1_active)

  if (!is.numeric(scenario$cens_rate) || length(scenario$cens_rate) != 1 ||
      !is.finite(scenario$cens_rate) || scenario$cens_rate < 0)
    stop("scenario$cens_rate must be a finite non-negative number; got ",
         scenario$cens_rate, ".")

  # Reconstruct full-length beta vectors under active_first layout.
  # Under active_first, columns 1..k are the active predictors in order.
  beta1 <- c(scenario$beta1_active, rep(0.0, p - k))
  beta2 <- c(scenario$beta2_active, rep(0.0, p - k))

  sim <- fastcmprsk::simulateTwoCauseFineGrayModel(
    nobs    = n,
    beta1   = beta1,
    beta2   = beta2,
    X       = x,
    u.min   = scenario$u_min,
    u.max   = scenario$u_max,
    p       = 0.5,
    returnX = TRUE
  )

  check_design_unchanged(supplied = x, returned = sim$X)

  # Exponential censoring
  ftime   <- sim$ftime
  fstatus <- sim$fstatus
  if (scenario$cens_rate > 0) {
    cens_time  <- stats::rexp(n, rate = scenario$cens_rate)
    is_censored <- cens_time < ftime
    ftime[is_censored]   <- cens_time[is_censored]
    fstatus[is_censored] <- 0L
  }

  event_mix <- c(
    censored = mean(fstatus == 0),
    cause1   = mean(fstatus == 1),
    cause2   = mean(fstatus == 2)
  )

  y <- matrix(c(ftime, as.integer(fstatus)), ncol = 2,
              dimnames = list(NULL, c("time", "status")))

  list(y = y, event_mix = event_mix)
}

# ---------------------------------------------------------------------------
# sim_data
# ---------------------------------------------------------------------------

sim_data <- function(scenario, p) {
  blocks      <- scenario$blocks
  block_names <- names(blocks)

  block_props <- stats::setNames(
    vapply(blocks, `[[`, numeric(1), "proportion"), block_names)
  block_rho_v <- stats::setNames(
    vapply(blocks, `[[`, numeric(1), "rho"),        block_names)
  latent      <- blocks[[1L]]$latent_type
  latent_q    <- stats::setNames(
    vapply(blocks, `[[`, numeric(1), "latent_q"),   block_names)

  # Campaign-fixed defaults (see header)
  design <- .generate_block_design(
    nobs             = scenario$nobs,
    npredictors      = p,
    beta1_active     = scenario$beta1_active,
    beta2_active     = scenario$beta2_active,
    active_block     = scenario$active_block,
    block_props      = block_props,
    block_rho        = block_rho_v,
    latent           = latent,
    latent_q         = latent_q,
    latent_balanced  = TRUE,
    unit_variance    = TRUE,
    layout           = "active_first",
    block_id         = NULL,
    predictor_prefix = "X_"
  )

  outcomes <- sim_outcomes(design$X, scenario)

  # Build block_map: rename meta$rho -> block_rho, drop index/beta1/beta2
  m <- design$meta
  block_map <- data.frame(
    name      = m$name,
    block     = m$block,
    block_rho = m$rho,    # A3 fix: the column is block_rho, not rho
    designed  = m$designed,
    active1   = m$active1,
    active2   = m$active2,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  design_settings <- c(
    design$settings,
    list(event_mix = outcomes$event_mix)
  )

  list(
    x               = design$X,
    y               = outcomes$y,
    beta1           = design$beta1,
    beta2           = design$beta2,
    block_map       = block_map,
    block_signature = design$signature,
    design_settings = design_settings
  )
}
