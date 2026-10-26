# sim/R/dgp.r -- block-correlated covariates and competing-risks outcomes
#
# OWNS
#   sim_covariates(n, p, blocks)   draw an n x p design matrix with
#                                  within-block Gaussian copula correlation;
#                                  a single all-independent block (rho = 0)
#                                  reproduces the simple fastcmprsk generator,
#                                  so no separate code path exists for the
#                                  uncorrelated case -- it is a configuration
#   sim_outcomes(x, scenario)      given covariates and a scenario spec,
#                                  draw (time, status) via subdistribution
#                                  hazards with cause-1 and cause-2 coefficients
#   sim_data(scenario, p)          convenience wrapper: covariates + outcomes
#                                  in one call, returns list(x, y)
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
#   sim_data:       list(x = matrix n x p, y = matrix n x 2)
#
# DEPENDS ON
#   sim/R/scenario.r  (scenario spec fields)
#   bhCRR::           (none; DGP is self-contained)
#
# REPLACES
#   sim/_legacy/code/generate_block_X.R  (entire file)
#   inline DGP logic in sim/_legacy/code/batch_run_complex.R lines 1-~80
#
# NOT YET IMPLEMENTED
