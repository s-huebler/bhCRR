# sim/scenarios/equalblocks.01.R
#
# DESIGN
#   High dimensional scenario intended for p = 1000 (npredictors is NOT a
#   signature key; it is passed on the run_sim.sh command line). Three blocks of
#   equal size (indep / weak / strong), 30 designed actives spread evenly across
#   them (10 per block).
#
#   |beta1| ~ Uniform(0.25, 0.75) with a random sign, so each active effect is
#   in [0.25, 0.75] or [-0.75, -0.25]. Magnitudes were softened from an earlier
#   draft (50 actives, U(0.5, 1.5)) whose linear predictor had sd ~5.9, which
#   made event type nearly deterministic for half the sample. beta2 = -beta1 via the "neg_beta1" sentinel.
#
# WHY THE VALUES ARE LITERALS, NOT runif()/sample() CALLS
#   This file is re-sourced by scenario_from_file() on every run and every
#   manifest check, and run_one() sources it BEFORE it calls set.seed(). A
#   random draw here would give each run (and each array task) different "true"
#   betas, and check_manifest() would stop on BETA1_ACTIVE. set.seed() is also
#   rejected by the forbidden-call scan. So the draws were made ONCE, offline,
#   and pasted in. They reproduce exactly with R's default RNG
#   (Mersenne-Twister, Inversion, Rejection) using:
#
#     set.seed(20261027)
#     k   <- 30L
#     mag <- round(runif(k, 0.25, 0.75), 2)
#     sgn <- sample(c(-1, 1), k, replace = TRUE)
#     beta1_active <- sgn * mag
#     active_block <- sample(rep(c("indep", "weak", "strong"), length.out = k))
#
#   To redraw, run that snippet with a new seed and paste the results in UNDER A
#   NEW SCENARIO NAME (beta1_active and active_block are signature keys).
#
# CONTRACT (from _template.R)
#   Assign exactly ONE object named `scenario`, a plain list; no library(),
#   require(), source(), set.seed(); no other side effects. Deterministic base R
#   (c(), rep(), seq(), round()) would be allowed, but only if no helper object
#   is left behind and nothing touches the RNG.
#
# SIGNATURE KEYS (12 total)
#   Changing ANY key below requires a new scenario name so that old and new
#   results are never mixed in sim/results/.

scenario <- list(

  # ---------------------------------------------------------------------------
  # name  character(1)  SIGNATURE KEY
  #   Scenario slug used in all output filenames and the results directory.
  #   Must match ^[A-Za-z0-9][A-Za-z0-9.-]*$ (see NAME RULE above).
  #   Consumed by: scenario_from_file(), write_manifest(), parse_results()
  #   [PLACEHOLDER -- new_scenario() fills this in automatically]
  # ---------------------------------------------------------------------------
  name            = "equalblocks.01",

  # ---------------------------------------------------------------------------
  # nobs  integer(1)  SIGNATURE KEY  units: observations  range: >= 1
  #   Total sample size per simulated dataset.
  #   Consumed by: dgp.r::sim_data()
  #   [PLACEHOLDER -- 100 is intentionally small for smoke-testing]
  # ---------------------------------------------------------------------------
  nobs            = 500L,

  # ---------------------------------------------------------------------------
  # zero_gap_target  numeric(1)  SIGNATURE KEY  range: (0, Inf)
  #   Minimum clinically-relevant effect size (passed as beta_min to
  #   bhcrr_autotune). Autotune will reject tuning grids that cannot
  #   distinguish this gap from zero. Units: same scale as beta1_active
  #   coefficients (subdistribution log-hazard ratio units).
  #   Consumed by: fit.r::fit_bhcrr() -> bhcrr_autotune()
  #   [PLACEHOLDER]
  # ---------------------------------------------------------------------------
  zero_gap_target = 0.10,

  # ---------------------------------------------------------------------------
  # beta1_active  numeric(k)  SIGNATURE KEY  range: any nonzero real
  #   Cause-1 Fine-Gray SUBDISTRIBUTION coefficients for the k designed active
  #   predictors. NOT cause-specific coefficients. These are the true nonzero
  #   entries; dgp.r pads them with zeros to length p.
  #   Length must equal length(active_block).
  #   Consumed by: dgp.r::sim_data(), score.r::score_run() (via beta1)
  #   30 actives, |beta| ~ U(0.25, 0.75), random sign (19 negative, 11 positive).
  #   Generated offline; see the header for the seed and recipe. Entry i is
  #   column X_i under the active_first layout.
  # ---------------------------------------------------------------------------
  beta1_active    = c(
     0.36, -0.70, -0.30,  0.45, -0.29, -0.49,  0.69,  0.65, -0.65, -0.61,  #  1-10
     0.33, -0.47, -0.72,  0.33, -0.44,  0.48, -0.47, -0.29,  0.60,  0.70,  # 11-20
    -0.66, -0.46,  0.69, -0.51,  0.37, -0.38, -0.45, -0.32, -0.34, -0.44   # 21-30
  ),

  # ---------------------------------------------------------------------------
  # beta2_active  numeric(k) OR "neg_beta1"  SIGNATURE KEY
  #   Cause-2 Fine-Gray SUBDISTRIBUTION coefficients for the active predictors.
  #   The sentinel "neg_beta1" means use -beta1_active; this matches the typical
  #   antagonistic relationship between competing risks in the simulation
  #   literature and is the recommended default unless you specifically want to
  #   model a different cause-2 pattern.
  #   If given as a vector, must have the same length as beta1_active.
  #   Consumed by: dgp.r::sim_data()
  #   Sentinel: beta2 = -beta1, so every cause-2 effect is also in
  #   [0.25, 0.75] in absolute value with the opposite sign.
  # ---------------------------------------------------------------------------
  beta2_active    = "neg_beta1",

  # ---------------------------------------------------------------------------
  # active_block  character(k)  SIGNATURE KEY
  #   Block label for each designed active predictor, in the same order as
  #   beta1_active. Every entry must appear in block_props / block_rho /
  #   latent_q. This encodes the intended correlation structure for the
  #   active predictors.
  #   Length must equal length(beta1_active).
  #   Consumed by: dgp.r::sim_data() (for the block_map)
  #   Balanced random assignment (10 per block), shuffled.
  #   Generated offline; see the header for the seed and recipe.
  # ---------------------------------------------------------------------------
  active_block    = c(
    "strong", "weak",   "weak",   "weak",   "strong", "strong", "weak",   "indep",  "strong", "strong",  #  1-10
    "strong", "indep",  "indep",  "strong", "weak",   "weak",   "weak",   "indep",  "indep",  "strong",  # 11-20
    "indep",  "indep",  "strong", "indep",  "indep",  "weak",   "strong", "weak",   "indep",  "weak"     # 21-30
  ),

  # ---------------------------------------------------------------------------
  # block_props  named numeric  SIGNATURE KEY  sum must equal 1  range: (0, 1) per entry
  #   Proportion of the p predictors assigned to each block type. The names
  #   are the BLOCK LABELS: arbitrary strings chosen by the scenario author.
  #   The SAME label set must appear in block_props, block_rho, and latent_q,
  #   and every entry of active_block must be one of these labels.
  #   Consumed by: dgp.r::sim_data() (to partition p predictors across blocks)
  #   [PLACEHOLDER]
  # ---------------------------------------------------------------------------
  block_props     = c(indep = 0.34, weak = 0.33, strong = 0.33),

  # ---------------------------------------------------------------------------
  # block_rho  named numeric  SIGNATURE KEY  range: [0, 1) per entry
  #   Within-block Pearson correlation (after unit-variance rescaling). 0 = no
  #   correlation (independent block). The label set must match block_props.
  #   The rho-to-latent-shift conversion is a = sqrt((rho/(1-rho)) / (q(1-q)));
  #   see dgp.r's NOTES FOR THE IMPLEMENTER.
  #   Consumed by: dgp.r::sim_covariates()
  #   [PLACEHOLDER]
  # ---------------------------------------------------------------------------
  block_rho       = c(indep = 0.00, weak = 0.30, strong = 0.70),

  # ---------------------------------------------------------------------------
  # latent_type  character(1)  SIGNATURE KEY  values: "gaussian" | "bernoulli"
  #   Distribution of the latent variable used to induce block correlation.
  #   "gaussian" is the default and covers most use cases.
  #   "bernoulli" replicates the source papers' binary-covariate structure;
  #   use it only when paper-faithful comparison is required.
  #   Consumed by: dgp.r::sim_covariates()
  # ---------------------------------------------------------------------------
  latent_type     = "gaussian",

  # ---------------------------------------------------------------------------
  # latent_q  named numeric  SIGNATURE KEY  range: (0, 1) per entry
  #   Marginal quantile of the latent variable used as the copula threshold.
  #   Controls the effective correlation strength together with block_rho.
  #   The label set must match block_props.
  #   Consumed by: dgp.r::sim_covariates()
  #   [PLACEHOLDER]
  # ---------------------------------------------------------------------------
  latent_q        = c(indep = 0.50, weak = 0.40, strong = 0.50),

  # ---------------------------------------------------------------------------
  # cens_rate  numeric(1)  SIGNATURE KEY  range: >= 0  units: Exp rate
  #   Rate parameter for independent exponential censoring applied after event
  #   generation. 0 = no censoring. Higher values = heavier censoring.
  #   Consumed by: dgp.r::sim_outcomes() via stats::rexp
  #   [PLACEHOLDER]
  # ---------------------------------------------------------------------------
  cens_rate       = 0.05,

  # ---------------------------------------------------------------------------
  # u_min  numeric(1)  SIGNATURE KEY  range: > 0, u_min <= u_max
  #   Lower bound of the observation window used in outcome simulation.
  #   Consumed by: dgp.r::sim_outcomes() via fastcmprsk::simulateTwoCauseFineGrayModel
  #   [PLACEHOLDER]
  # ---------------------------------------------------------------------------
  u_min           = 50,

  # ---------------------------------------------------------------------------
  # u_max  numeric(1)  SIGNATURE KEY  range: >= u_min
  #   Upper bound of the observation window used in outcome simulation.
  #   Consumed by: dgp.r::sim_outcomes() via fastcmprsk::simulateTwoCauseFineGrayModel
  #   [PLACEHOLDER]
  # ---------------------------------------------------------------------------
  u_max           = 150

)
