###############################################################################
## generate_block_X()
##
## Block-correlated design matrix for the competing-risks simulation, following
## the covariate structure of
##
##   Binder H, Schumacher M (2008). Adapting prediction error estimates for
##   biased complexity selection in high-dimensional bootstrap samples.
##   Stat Appl Genet Mol Biol 7(1), Article 12. Section 3.1.
##
## as reused for competing risks in
##
##   Binder H, Allignol A, Schumacher M, Beyersmann J (2009). Boosting for
##   high-dimensional time-to-event data with competing risks.
##   Bioinformatics 25(7), 890-896. Section 3.2.1.
##
## STRUCTURE. Predictors are partitioned into named blocks. Within a block every
## predictor shares one per-subject latent draw and then adds its own N(0, 1)
## noise. That induces exchangeable (compound-symmetric) correlation inside a
## block and exactly zero correlation between blocks, since each block draws its
## own independent latent:
##
##     z_ij = shift_ib + eps_ij,   Var(shift_ib) = v_b,   eps_ij ~ N(0, 1)
##     rho_b = v_b / (1 + v_b)      <=>      v_b = rho_b / (1 - rho_b)
##
## SETTING RHO DIRECTLY. Both papers hard-code the shared shift as
## a_b * I(u_ib < q_b) and then tabulate the correlation it happens to produce.
## Inverting that relation,
##
##     a_b = sqrt( v_b / (q_b * (1 - q_b)) ),
##
## makes rho the parameter the caller sets rather than a downstream consequence.
## The inversion reproduces the published constants exactly, which is the check
## that this really is the papers' generator with a different front end:
##
##     rho 0.50, q 0.5  ->  a = 2.000   (paper: 2, the -1 / +1 split)
##     rho 0.35, q 0.4  ->  a = 1.498   (paper: 1.5)
##     rho 0.05, q 0.7  ->  a = 0.501   (paper: 0.5)
##     rho 0.32, q 0.3  ->  a = 1.499   (paper: 1.5)
##
## LATENT TYPE. latent = "bernoulli" is the papers' construction and yields the
## bimodal, microarray-like marginals they intended. latent = "gaussian"
## replaces the Bernoulli shift with sqrt(v_b) * N(0, 1), which hits rho exactly
## in every sample rather than in expectation and needs no special case at
## rho = 0. Both run through the identical shift-plus-noise scaffold below.
##
## VARIANCE SCALING. The papers' raw covariates have marginal variance
## 1 + v_b, so a coefficient of 0.5 carries about 1.4x the signal in a rho = 0.5
## block as in an uncorrelated one, purely from scale. unit_variance = TRUE
## rescales each block to unit marginal variance, which leaves the correlations
## untouched (dividing a block by a constant cannot change them) and puts every
## beta on a common scale so bias is comparable across blocks. This is a
## deliberate departure from Binder & Schumacher.
##
## RELATION TO THE EVENT-TIME MODEL. This function supplies covariates only.
## The papers generate cause-specific hazards from a Cox-exponential model,
## whereas simulateTwoCauseFineGrayModel() generates from the Fine-Gray
## subdistribution. beta1 here stays a subdistribution coefficient and is NOT
## the same object as Binder's alpha_1. Same design matrix, different mechanism.
##
## RNG. No seed is set here. Set one in the calling script if a scenario needs
## to be reproducible.
###############################################################################

generate_block_X <- function(nobs,
                             npredictors,
                             beta1_active     = numeric(0),
                             beta2_active     = numeric(0),
                             active_block     = character(0),
                             block_props      = c(indep = 1/3, weak = 1/3, strong = 1/3),
                             block_rho        = c(indep = 0, weak = 0.3, strong = 0.6),
                             latent           = c("gaussian", "bernoulli"),
                             latent_q         = c(indep = 0.5, weak = 0.4, strong = 0.5),
                             latent_balanced  = TRUE,
                             unit_variance    = TRUE,
                             layout           = c("active_first", "contiguous"),
                             predictor_prefix = "X_") {

  latent <- match.arg(latent)
  layout <- match.arg(layout)

  ## 1. Scalar inputs -------------------------------------------------------
  if (length(nobs) != 1L || !is.finite(nobs) || nobs < 2L) {
    stop("nobs must be a single finite integer of at least 2.")
  }
  if (length(npredictors) != 1L || !is.finite(npredictors) || npredictors < 1L) {
    stop("npredictors must be a single finite positive integer.")
  }
  nobs <- as.integer(nobs)
  npredictors <- as.integer(npredictors)

  ## 2. Block specification -------------------------------------------------
  ## block_props, block_rho and (for the Bernoulli latent) latent_q must all be
  ## named with the same block names. Everything downstream is indexed by name,
  ## never by position, so the caller cannot silently misalign rho with a block.
  block_names <- names(block_props)
  if (is.null(block_names) || any(!nzchar(block_names)) ||
      anyDuplicated(block_names) > 0L) {
    stop("block_props must be a named vector with unique, non-empty block names.")
  }
  if (any(!is.finite(block_props)) || any(block_props < 0) || sum(block_props) <= 0) {
    stop("block_props must be finite, non-negative, and not all zero.")
  }
  if (!setequal(names(block_rho), block_names)) {
    stop("block_rho must be named with exactly the same blocks as block_props: ",
         paste(block_names, collapse = ", "), ".")
  }
  block_rho <- block_rho[block_names]
  if (any(!is.finite(block_rho)) || any(block_rho < 0) || any(block_rho >= 1)) {
    stop("block_rho values must be finite and in [0, 1).")
  }
  if (latent == "bernoulli") {
    if (!setequal(names(latent_q), block_names)) {
      stop("latent_q must be named with exactly the same blocks as block_props.")
    }
    latent_q <- latent_q[block_names]
    if (any(!is.finite(latent_q)) || any(latent_q <= 0) || any(latent_q >= 1)) {
      stop("latent_q values must be finite and strictly between 0 and 1.")
    }
  }

  ## 3. Active-predictor specification --------------------------------------
  ## beta1_active, beta2_active and active_block are aligned elementwise: entry
  ## k is one designed predictor, with its cause-1 effect, its cause-2 effect,
  ## and the block it is planted in. An entry stays "designed active" even where
  ## one of its two betas is zero, which is how Binder et al. place covariates
  ## that act on one cause only.
  n_active <- length(beta1_active)
  if (length(beta2_active) != n_active || length(active_block) != n_active) {
    stop("beta1_active, beta2_active and active_block must all have the same ",
         "length (got ", n_active, ", ", length(beta2_active), ", ",
         length(active_block), ").")
  }
  if (n_active > npredictors) {
    stop("More active predictors (", n_active, ") than predictors (",
         npredictors, ").")
  }
  if (any(!is.finite(c(beta1_active, beta2_active)))) {
    stop("beta1_active and beta2_active must be finite.")
  }
  active_block <- as.character(active_block)
  unknown_block <- setdiff(active_block, block_names)
  if (length(unknown_block)) {
    stop("active_block names blocks that do not exist: ",
         paste(unique(unknown_block), collapse = ", "), ".")
  }

  ## 4. Block sizes ---------------------------------------------------------
  ## Largest-remainder (Hamilton) apportionment so the block sizes sum to
  ## exactly npredictors. A block carrying rho > 0 needs at least two members
  ## for a within-block correlation to exist at all; an uncorrelated block needs
  ## one. With the default three blocks that floor is 5 predictors.
  target <- block_props / sum(block_props) * npredictors
  block_size <- floor(target)
  shortfall <- npredictors - sum(block_size)
  if (shortfall > 0L) {
    bump <- order(target - floor(target), decreasing = TRUE)[seq_len(shortfall)]
    block_size[bump] <- block_size[bump] + 1L
  }
  names(block_size) <- block_names

  min_size <- ifelse(block_rho > 0, 2L, 1L)
  names(min_size) <- block_names
  if (npredictors < sum(min_size)) {
    stop("npredictors = ", npredictors, " is too small for this block ",
         "specification; it needs at least ", sum(min_size),
         " (2 per correlated block, 1 per uncorrelated block).")
  }
  adjusted <- FALSE
  repeat {
    short <- which(block_size < min_size)
    if (!length(short)) break
    donor <- which.max(block_size - min_size)
    if (block_size[donor] - min_size[donor] <= 0L) {
      stop("Cannot satisfy the minimum size of every block at npredictors = ",
           npredictors, "; raise npredictors or change block_props.")
    }
    block_size[short[1L]] <- block_size[short[1L]] + 1L
    block_size[donor] <- block_size[donor] - 1L
    adjusted <- TRUE
  }
  if (adjusted) {
    warning("block_props were adjusted so every block meets its minimum size. ",
            "Realized sizes: ",
            paste(block_names, block_size, sep = "=", collapse = ", "), ".")
  }

  active_per_block <- table(factor(active_block, levels = block_names))
  overfull <- which(as.integer(active_per_block) > block_size)
  if (length(overfull)) {
    stop("Block(s) ", paste(block_names[overfull], collapse = ", "),
         " were assigned more active predictors than they have columns.")
  }

  ## 5. Column layout -------------------------------------------------------
  ## "active_first" puts the designed actives in columns 1..n_active, in the
  ## order the caller wrote them, so beta1[seq_along(beta1_active)] lines up
  ## with beta1_active and downstream code that assumes leading actives keeps
  ## working. Blocks are then non-contiguous in column order, which costs
  ## nothing: correlation is defined by block membership, not adjacency.
  ## "contiguous" reproduces the papers' column ordering instead.
  block_id <- character(npredictors)
  active_idx <- integer(n_active)

  if (layout == "active_first") {
    if (n_active > 0L) {
      block_id[seq_len(n_active)] <- active_block
      active_idx <- seq_len(n_active)
    }
    filler_counts <- block_size - as.integer(active_per_block)
    if (npredictors > n_active) {
      block_id[(n_active + 1L):npredictors] <-
        rep(block_names, times = filler_counts)
    }
  } else {
    block_id <- rep(block_names, times = block_size)
    next_free <- cumsum(c(0L, block_size[-length(block_size)])) + 1L
    names(next_free) <- block_names
    for (k in seq_len(n_active)) {
      b <- active_block[k]
      active_idx[k] <- next_free[[b]]
      next_free[[b]] <- next_free[[b]] + 1L
    }
  }

  ## 6. Draw the covariates -------------------------------------------------
  ## One shared latent per block per subject, plus independent N(0, 1) noise per
  ## cell. The shift vector recycles down the columns of the block, which is
  ## what makes the correlation exchangeable.
  X <- matrix(0.0, nrow = nobs, ncol = npredictors)
  shift_sd <- setNames(sqrt(block_rho / (1 - block_rho)), block_names)
  latent_a <- setNames(rep(NA_real_, length(block_names)), block_names)

  for (b in block_names) {
    cols <- which(block_id == b)
    if (!length(cols)) next

    v <- shift_sd[[b]]^2
    if (v <= 0) {
      shift <- rep(0.0, nobs)
      latent_a[[b]] <- 0
    } else if (latent == "gaussian") {
      shift <- shift_sd[[b]] * stats::rnorm(nobs)
    } else {
      q <- latent_q[[b]]
      a <- sqrt(v / (q * (1 - q)))
      latent_a[[b]] <- a
      ## latent_balanced fixes the number of shifted subjects at round(q * nobs)
      ## rather than drawing each independently, matching the deterministic
      ## i <= 0.5n split Binder & Schumacher use in their first block and making
      ## rho exact in the sample rather than in expectation.
      ind <- if (latent_balanced) {
        n1 <- round(q * nobs)
        sample(rep(c(1.0, 0.0), times = c(n1, nobs - n1)))
      } else {
        as.numeric(stats::runif(nobs) < q)
      }
      ## Centre on the nominal q so the block has population mean zero. Any
      ## residual mean is absorbed by the baseline hazard downstream.
      shift <- a * (ind - q)
    }

    eps <- matrix(stats::rnorm(nobs * length(cols)), nrow = nobs)
    Z <- shift + eps

    ## 7. Optional rescale to unit marginal variance ------------------------
    ## Var(z) = 1 + v = 1 / (1 - rho), so multiplying by sqrt(1 - rho) gives
    ## sd 1 without touching the correlation.
    if (unit_variance) {
      Z <- Z * sqrt(1 - block_rho[[b]])
    }
    X[, cols] <- Z
  }

  ## 8. Scatter the betas to their columns ----------------------------------
  ## Built by index rather than by c(beta_active, rep(0, ...)) so it stays
  ## correct under either layout. which(beta1 != 0) downstream still recovers
  ## the cause-1 oracle set, which is a subset of the designed actives.
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

  ## 9. Per-predictor metadata ----------------------------------------------
  ## Carried out of the function so results can be summarised by correlation
  ## stratum, which is the whole point of the block design.
  meta <- data.frame(
    index    = seq_len(npredictors),
    name     = predictor_names,
    block    = block_id,
    rho      = as.numeric(block_rho[block_id]),
    beta1    = unname(beta1),
    beta2    = unname(beta2),
    designed = seq_len(npredictors) %in% active_idx,
    active1  = unname(beta1 != 0),
    active2  = unname(beta2 != 0),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  ## 10. Return -------------------------------------------------------------
  list(
    X          = X,
    block_id   = block_id,
    block_size = block_size,
    block_rho  = block_rho,
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


###############################################################################
## block_corr_summary()
##
## Validation helper. Reports the realized mean within-block correlation against
## the target and the largest absolute cross-block correlation, which should sit
## near zero. Run it at large nobs when checking a new block specification.
###############################################################################

block_corr_summary <- function(X, block_id, block_rho = NULL) {
  if (ncol(X) != length(block_id)) {
    stop("block_id must have one entry per column of X.")
  }
  R <- stats::cor(X)
  blocks <- unique(block_id)

  within <- vapply(blocks, function(b) {
    cols <- which(block_id == b)
    if (length(cols) < 2L) return(NA_real_)
    mean(R[cols, cols][upper.tri(diag(length(cols)))])
  }, numeric(1))

  cross_max <- vapply(blocks, function(b) {
    cols <- which(block_id == b)
    other <- setdiff(seq_along(block_id), cols)
    if (!length(other)) return(NA_real_)
    max(abs(R[cols, other]))
  }, numeric(1))

  out <- data.frame(
    block        = blocks,
    n_predictors = as.integer(table(factor(block_id, levels = blocks))),
    target_rho   = if (is.null(block_rho)) NA_real_ else as.numeric(block_rho[blocks]),
    mean_within  = as.numeric(within),
    max_cross    = as.numeric(cross_max),
    row.names    = NULL,
    stringsAsFactors = FALSE
  )
  out$sd_check <- vapply(blocks, function(b) {
    mean(apply(X[, block_id == b, drop = FALSE], 2L, stats::sd))
  }, numeric(1))
  out
}
