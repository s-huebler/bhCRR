# sim/scenarios/_template.R -- copy this file with new_scenario() to start a new scenario
#
# USAGE
#   new_scenario("myScenario.01", dir = "sim/scenarios")
#   # Opens sim/scenarios/myScenario.01.R with name already filled in.
#   # Then replace EVERY value marked [PLACEHOLDER] with a real number.
#
# CONTRACT
#   This file must:
#     - assign exactly ONE object named `scenario`, a plain list
#     - contain NO calls to library(), require(), source(), or set.seed()
#     - produce NO other side effects
#   scenario_from_file() enforces all of these before any validation runs.
#
# NAME RULE
#   `name` must match ^[A-Za-z0-9][A-Za-z0-9.-]*$  (letters, digits, dots,
#   hyphens; NO underscores). Underscores are field delimiters in the
#   p<P>_<scenario>_run<N> filename format used by run_sim.sh and parse.r;
#   a name containing '_' would break filename parsing. Use '.' or '-' as
#   word separators instead.
#   `name` must EQUAL the file's basename without extension (enforced by
#   scenario_from_file). A mismatch would silently write results to the wrong
#   directory.
#
# SIGNATURE KEYS (12 total)
#   Changing ANY key below requires a new scenario name so that old and new
#   results are never mixed in sim/results/. Non-signature parameters
#   (npredictors, run range) are passed on the run_sim.sh command line
#   and are NOT stored here.
#
# PLACEHOLDER VALUES
#   The values below are deliberately small and round so the scenario runs fast
#   for smoke-testing. They are NOT scientifically motivated.
#   Replace every value marked [PLACEHOLDER] before real use.

scenario <- list(

  # ---------------------------------------------------------------------------
  # name  character(1)  SIGNATURE KEY
  #   Scenario slug used in all output filenames and the results directory.
  #   Must match ^[A-Za-z0-9][A-Za-z0-9.-]*$ (see NAME RULE above).
  #   Consumed by: scenario_from_file(), write_manifest(), parse_results()
  #   [PLACEHOLDER -- new_scenario() fills this in automatically]
  # ---------------------------------------------------------------------------
  name            = "troubleshoot.01",

  # ---------------------------------------------------------------------------
  # nobs  integer(1)  SIGNATURE KEY  units: observations  range: >= 1
  #   Total sample size per simulated dataset.
  #   Consumed by: dgp.r::sim_data()
  #   [PLACEHOLDER -- 100 is intentionally small for smoke-testing]
  # ---------------------------------------------------------------------------
  nobs            = 100L,

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
  #   [PLACEHOLDER -- 3 actives]
  # ---------------------------------------------------------------------------
  beta1_active    = c(0.5, -0.5, 0.5),

  # ---------------------------------------------------------------------------
  # beta2_active  numeric(k) OR "neg_beta1"  SIGNATURE KEY
  #   Cause-2 Fine-Gray SUBDISTRIBUTION coefficients for the active predictors.
  #   The sentinel "neg_beta1" means use -beta1_active; this matches the typical
  #   antagonistic relationship between competing risks in the simulation
  #   literature and is the recommended default unless you specifically want to
  #   model a different cause-2 pattern.
  #   If given as a vector, must have the same length as beta1_active.
  #   Consumed by: dgp.r::sim_data()
  #   [PLACEHOLDER -- using sentinel]
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
  #   [PLACEHOLDER]
  # ---------------------------------------------------------------------------
  active_block    = c("indep" , "indep", "indep"),

  # ---------------------------------------------------------------------------
  # block_props  named numeric  SIGNATURE KEY  sum must equal 1  range: (0, 1) per entry
  #   Proportion of the p predictors assigned to each block type. The names
  #   are the BLOCK LABELS: arbitrary strings chosen by the scenario author.
  #   The SAME label set must appear in block_props, block_rho, and latent_q,
  #   and every entry of active_block must be one of these labels.
  #   Consumed by: dgp.r::sim_data() (to partition p predictors across blocks)
  #   [PLACEHOLDER]
  # ---------------------------------------------------------------------------
  block_props     = c(indep = 1, weak = 0.00, strong = 0.00),

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
