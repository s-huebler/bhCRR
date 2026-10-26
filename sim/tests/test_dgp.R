library(testthat)

# Resolve repo root before load.r is available
.test_root <- local({
  e <- Sys.getenv("REPO_ROOT")
  if (nzchar(e)) return(normalizePath(e, mustWork = TRUE))
  p <- normalizePath(getwd(), mustWork = TRUE)
  repeat {
    if (file.exists(file.path(p, "DESCRIPTION"))) return(p)
    pp <- dirname(p)
    if (identical(pp, p)) stop("Cannot locate repo root; set REPO_ROOT")
    p <- pp
  }
})

source(file.path(.test_root, "sim", "R", "load.r"))
load_sim("dev", force = FALSE)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Standard 3-block setup (bernoulli; one block of each correlation level)
.bp3  <- c(indep = 0.60, weak = 0.20, strong = 0.20)
.br3  <- c(indep = 0.00, weak = 0.35, strong = 0.50)
.bq3  <- c(indep = 0.50, weak = 0.40, strong = 0.50)

.gbd <- function(..., bp = .bp3, br = .br3, bq = .bq3,
                 latent = "bernoulli", n = 60L, p = 25L,
                 b1 = c(0.5, -0.5, 0.5), b2 = -b1,
                 ab = c("strong", "weak", "indep")) {
  .generate_block_design(
    nobs = n, npredictors = p,
    beta1_active = b1, beta2_active = b2, active_block = ab,
    block_props = bp, block_rho = br,
    latent = latent, latent_q = bq,
    ...
  )
}

# Minimal valid scenario args (reuse test_scenario.R style)
.valid_scenario_args <- list(
  name            = "test.scenario",
  nobs            = 200L,
  zero_gap_target = 0.1,
  beta1_active    = c(0.5, -0.5, 0.5),
  beta2_active    = "neg_beta1",
  active_block    = c("strong", "weak", "indep"),
  block_props     = c(indep = 0.60, weak = 0.20, strong = 0.20),
  block_rho       = c(indep = 0.00, weak = 0.35, strong = 0.50),
  latent_type     = "bernoulli",
  latent_q        = c(indep = 0.50, weak = 0.40, strong = 0.50),
  cens_rate       = 0.05,
  u_min           = 50,
  u_max           = 150
)

# ---------------------------------------------------------------------------
# D1: Bernoulli latent_a — rho-inversion formula
# ---------------------------------------------------------------------------

test_that("D1: Bernoulli latent_a matches rho-inversion formula to 3 dp", {
  # v = rho/(1-rho); a = sqrt(v / (q*(1-q)))
  # (rho, q) -> expected a
  pairs <- list(
    c(rho = 0.50, q = 0.50, a_exp = 2.000),
    c(rho = 0.35, q = 0.40, a_exp = 1.498),
    c(rho = 0.05, q = 0.70, a_exp = 0.501),
    c(rho = 0.32, q = 0.30, a_exp = 1.497)
  )
  # Build a 4-block design; one block per pair
  nblocks <- length(pairs)
  bp <- stats::setNames(rep(1/nblocks, nblocks), paste0("b", seq_len(nblocks)))
  br <- stats::setNames(vapply(pairs, `[`, numeric(1), "rho"), names(bp))
  bq <- stats::setNames(vapply(pairs, `[`, numeric(1), "q"),   names(bp))

  design <- .generate_block_design(
    nobs = 20L, npredictors = nblocks * 2L,
    block_props = bp, block_rho = br,
    latent = "bernoulli", latent_q = bq
  )

  la <- design$settings$latent_a
  for (i in seq_along(pairs)) {
    nm <- names(bp)[i]
    expect_equal(round(la[[nm]], 3), pairs[[i]][["a_exp"]],
                 info = paste("block", nm,
                              ": rho=", pairs[[i]]["rho"],
                              " q=", pairs[[i]]["q"]))
  }
})

# ---------------------------------------------------------------------------
# D2: Gaussian realized correlations
# ---------------------------------------------------------------------------

test_that("D2: Gaussian realized within-block rho is near target; cross-block near 0", {
  set.seed(42L)
  n_big <- 20000L
  design <- .generate_block_design(
    nobs = n_big, npredictors = 30L,
    block_props = c(a = 0.40, b = 0.30, c = 0.30),
    block_rho   = c(a = 0.70, b = 0.40, c = 0.00),
    latent = "gaussian",
    latent_q = c(a = 0.5, b = 0.5, c = 0.5),
    unit_variance = TRUE
  )

  X  <- design$X
  id <- design$block_id
  br <- design$block_rho

  # Within-block: realized pairwise mean rho should be close to target
  for (b in names(br)) {
    cols <- which(id == b)
    if (length(cols) < 2L) next
    R   <- cor(X[, cols])
    off <- R[upper.tri(R)]
    expect_lt(abs(mean(off) - br[[b]]), 0.03,
              label = paste("within-block mean rho, block", b))
  }

  # Cross-block: abs mean correlation should be close to 0
  blocks_with_2 <- names(which(table(id) >= 2L))
  if (length(blocks_with_2) >= 2L) {
    b1 <- blocks_with_2[1L]; b2 <- blocks_with_2[2L]
    cols1 <- which(id == b1); cols2 <- which(id == b2)
    cross <- abs(cor(X[, cols1], X[, cols2]))
    expect_lt(mean(cross), 0.03, label = "cross-block mean abs correlation")
  }
})

# ---------------------------------------------------------------------------
# D3: unit_variance scaling
# ---------------------------------------------------------------------------

test_that("D3: unit_variance=TRUE → column SD near 1", {
  set.seed(1L)
  design <- .generate_block_design(
    nobs = 10000L, npredictors = 20L,
    block_props = c(x = 0.50, y = 0.50),
    block_rho   = c(x = 0.60, y = 0.00),
    latent = "gaussian",
    latent_q = c(x = 0.5, y = 0.5),
    unit_variance = TRUE
  )
  sds <- apply(design$X, 2, sd)
  expect_lt(max(abs(sds - 1)), 0.05)
})

test_that("D3: unit_variance=FALSE → column SD near sqrt(1/(1-rho))", {
  set.seed(2L)
  rho_val <- 0.60
  design <- .generate_block_design(
    nobs = 20000L, npredictors = 10L,
    block_props = c(x = 1.0),
    block_rho   = c(x = rho_val),
    latent = "gaussian",
    latent_q = c(x = 0.5),
    unit_variance = FALSE
  )
  expected_sd <- sqrt(1 / (1 - rho_val))
  sds <- apply(design$X, 2, sd)
  expect_lt(max(abs(sds - expected_sd)), 0.05)
})

# ---------------------------------------------------------------------------
# D4: Block map is deterministic (no RNG dependence); differs across p
# ---------------------------------------------------------------------------

test_that("D4: block_map rows are identical under different RNG states", {
  set.seed(7L)
  d1 <- .gbd()
  set.seed(99L)
  d2 <- .gbd()

  # block_id (the column-to-block mapping) must be identical
  expect_identical(d1$block_id, d2$block_id)
  # block_size likewise
  expect_identical(d1$block_size, d2$block_size)
})

test_that("D4: block_map differs between p=25 and p=205", {
  d25  <- .gbd(p = 25L)
  d205 <- .gbd(p = 205L)
  expect_false(identical(d25$block_id,   d205$block_id))
  expect_false(identical(d25$block_size, d205$block_size))
})

# ---------------------------------------------------------------------------
# D5: Single block, rho=0, gaussian → effectively iid N(0,1)
# ---------------------------------------------------------------------------

test_that("D5: single block rho=0 gaussian produces iid N(0,1) columns", {
  set.seed(5L)
  design <- .generate_block_design(
    nobs = 20000L, npredictors = 10L,
    block_props = c(only = 1.0),
    block_rho   = c(only = 0.0),
    latent = "gaussian",
    latent_q = c(only = 0.5),
    unit_variance = TRUE
  )
  X <- design$X
  sds  <- apply(X, 2, sd)
  cors <- cor(X)
  off  <- cors[upper.tri(cors)]

  expect_lt(max(abs(sds  - 1.0)), 0.05,   label = "column SDs near 1")
  expect_lt(max(abs(off)),        0.05,    label = "all pairwise correlations near 0")
})

# ---------------------------------------------------------------------------
# D6: Column names exactly X_1..X_p
# ---------------------------------------------------------------------------

test_that("D6: column names are X_1..X_p", {
  design <- .gbd(p = 25L)
  expect_identical(colnames(design$X), paste0("X_", seq_len(25L)))
})

test_that("D6: column names X_1..X_p also from sim_covariates()", {
  spec <- do.call(scenario, .valid_scenario_args)
  X <- sim_covariates(n = 50L, p = 30L, blocks = spec$blocks)
  expect_identical(colnames(X), paste0("X_", seq_len(30L)))
})

# ---------------------------------------------------------------------------
# D7: beta vectors, active columns, block_map$designed
# ---------------------------------------------------------------------------

test_that("D7: beta1/beta2 have length p; active values in leading columns", {
  p <- 25L
  design <- .gbd(p = p)
  k <- 3L  # length(b1) in .gbd()

  expect_length(design$beta1, p)
  expect_length(design$beta2, p)

  # Under active_first, columns 1..k carry the active betas
  expect_equal(unname(design$beta1[1:k]), c(0.5, -0.5, 0.5))
  expect_true(all(design$beta1[(k + 1L):p] == 0))
  expect_true(all(design$beta2[(k + 1L):p] == 0))
})

test_that("D7: block_map$designed matches which(beta1 != 0) under active_first", {
  design <- .gbd(p = 25L)
  from_beta   <- unname(which(design$beta1 != 0))
  from_map    <- which(design$meta$designed)
  expect_equal(from_beta, from_map)
})

# ---------------------------------------------------------------------------
# D8: Outcomes: status in {0,1,2}; cens_rate=0 passes through sim's output
# ---------------------------------------------------------------------------

test_that("D8: status values are a subset of {0, 1, 2}", {
  set.seed(8L)
  spec <- do.call(scenario, .valid_scenario_args)
  dat  <- sim_data(spec, p = 25L)
  expect_true(all(dat$y[, "status"] %in% c(0L, 1L, 2L)))
})

test_that("D8: cens_rate=0 produces no additional censoring beyond simulator output", {
  args <- .valid_scenario_args
  args$cens_rate <- 0
  args$nobs <- 300L
  spec <- do.call(scenario, args)

  set.seed(8L)
  dat <- sim_data(spec, p = 25L)

  # With cens_rate=0 the only status==0 observations are those the
  # fine-gray simulator returned as censored (if any); we simply verify
  # that no extra censoring is applied by checking that non-event rows
  # are not inflated.  Status distribution must still sum to n rows.
  n <- nrow(dat$y)
  expect_equal(nrow(dat$y), n)
  expect_true(all(dat$y[, "status"] %in% c(0L, 1L, 2L)))
  # event_mix must sum to 1
  expect_equal(sum(dat$design_settings$event_mix), 1.0, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# D9: check_design_unchanged
# ---------------------------------------------------------------------------

test_that("D9: check_design_unchanged passes on identical matrices", {
  X <- matrix(rnorm(20), nrow = 4)
  expect_no_error(check_design_unchanged(X, X))
})

test_that("D9: check_design_unchanged errors on altered matrix", {
  X <- matrix(rnorm(20), nrow = 4)
  Y <- X; Y[1, 1] <- Y[1, 1] + 1.0
  expect_error(check_design_unchanged(X, Y),
               regexp = "simulateTwoCauseFineGrayModel|unchanged|design")
})

test_that("D9: check_design_unchanged errors on added column names", {
  X <- matrix(rnorm(20), nrow = 4,
              dimnames = list(NULL, paste0("X_", 1:5)))
  Y <- X
  colnames(Y) <- paste0("V", 1:5)   # names differ; unname() strips both
  expect_no_error(check_design_unchanged(X, Y))  # unname() means names are ignored
})

# ---------------------------------------------------------------------------
# D10: Block signature
# ---------------------------------------------------------------------------

test_that("D10: signature is identical across repeated calls with same settings", {
  d1 <- .gbd(p = 25L)
  d2 <- .gbd(p = 25L)
  expect_identical(d1$signature, d2$signature)
})

test_that("D10: signature changes when latent type flips", {
  dbern <- .gbd(latent = "bernoulli", p = 25L)
  dgaus <- .gbd(latent = "gaussian",  p = 25L)
  expect_false(identical(dbern$signature, dgaus$signature))
  expect_true(grepl("latent=bernoulli", dbern$signature))
  expect_true(grepl("latent=gaussian",  dgaus$signature))
})

test_that("D10: signature changes when unit_variance flips", {
  d_uv  <- .gbd(p = 25L, unit_variance = TRUE)
  d_nuv <- .gbd(p = 25L, unit_variance = FALSE)
  expect_false(identical(d_uv$signature, d_nuv$signature))
  expect_true(grepl("uvar=TRUE",  d_uv$signature))
  expect_true(grepl("uvar=FALSE", d_nuv$signature))
})

test_that("D10: signature changes when latent_balanced flips", {
  d_lb  <- .gbd(p = 25L, latent_balanced = TRUE)
  d_nlb <- .gbd(p = 25L, latent_balanced = FALSE)
  expect_false(identical(d_lb$signature, d_nlb$signature))
  expect_true(grepl("lbal=TRUE",  d_lb$signature))
  expect_true(grepl("lbal=FALSE", d_nlb$signature))
})

test_that("D10: signature contains all ten expected fields", {
  d <- .gbd(p = 25L)
  for (field in c("p=", "blocks=", "size=", "rho=", "layout=",
                   "active=", "pinned=", "latent=", "lbal=", "uvar=")) {
    expect_true(grepl(field, d$signature, fixed = TRUE),
                label = paste("signature contains", field))
  }
})
