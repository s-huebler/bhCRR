# sim/R/dgp.r -- block-correlated covariates and competing-risks outcomes
#
# OWNS
#   sim_covariates(n, p, blocks)   draw an n x p design matrix with
#                                  within-block Gaussian copula correlation;
#                                  a single all-independent block (rho = 0)
#                                  reproduces the simple fastcmprsk generator,
#                                  so no separate code path exists for the
#                                  uncorrelated case -- it is a configuration
#   sim_outcomes(x, scenario)      given covariates and a scenario spec, draw
#                                  (time, status) via
#                                  fastcmprsk::simulateTwoCauseFineGrayModel
#                                  with exponential censoring via stats::rexp
#                                  applied after event generation
#   sim_data(scenario, p)          convenience wrapper: covariates + outcomes
#                                  in one call; also constructs the block map
#                                  and block signature
#
# INPUTS
#   n         integer(1)   number of observations
#   p         integer(1)   total number of predictors
#   blocks    list         from scenario$blocks; each element holds rho,
#                          proportion, latent_type, latent_q, and the indices
#                          of the active predictors in that block
#   x         matrix n x p numeric design matrix from sim_covariates()
#   scenario  list         validated scenario spec from scenario.r
#
# RETURNS
#   sim_covariates: numeric matrix n x p, column names V1..Vp
#   sim_outcomes:   numeric matrix n x 2, colnames c("time", "status"),
#                   status in {0L, 1L, 2L}
#   sim_data:       list with elements
#                     x              matrix n x p
#                     y              matrix n x 2
#                     beta1          numeric p  full cause-1 coefficient vector
#                                    (active values at named columns, zeros
#                                    elsewhere)
#                     beta2          numeric p  full cause-2 coefficient vector
#                                    (same layout as beta1)
#                     block_map      data.frame one row per predictor:
#                                      name, block, block_rho, designed,
#                                      active1, active2
#                     block_signature character(1) fingerprint of the realized
#                                    column-to-block map; passed to $meta by
#                                    run_one.r
#
# NOTES FOR THE IMPLEMENTER
#   rho inversion: a = sqrt((rho/(1-rho)) / (q(1-q))) is what turns the
#   source papers' hard-coded latent shift into a settable rho. This formula
#   must be the FIRST thing tested if the block construction is ever rewritten.
#
#   unit-variance block rescaling is ON by default, a deliberate departure
#   from Binder & Schumacher, so that a given beta carries the same signal
#   in correlated and uncorrelated blocks.
#
#   the Gaussian latent is the default; Bernoulli exists for paper-faithful
#   replication of the source papers' covariate structure.
#
#   layout = "active_first" puts designed actives in the leading columns;
#   blocks are then non-contiguous, which costs nothing because correlation
#   is defined by block membership, not by column adjacency.
#
#   beta1 is a Fine-Gray SUBDISTRIBUTION coefficient, not a cause-specific
#   one. The papers' covariate structure is replicated; their hazard model
#   is not.
#
# DEPENDS ON
#   sim/R/scenario.r                            (scenario spec fields)
#   fastcmprsk::simulateTwoCauseFineGrayModel   (outcome generation; block-
#                                                designed X is passed in)
#   stats::rexp                                 (exponential censoring applied
#                                                after event generation)
#
# REPLACES
#   sim/_legacy/code/generate_block_X.R  (entire file)
#   inline DGP logic in sim/_legacy/code/batch_run_complex.R lines 1-~80
#
# NOT YET IMPLEMENTED
