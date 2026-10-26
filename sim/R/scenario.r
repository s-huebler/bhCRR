# sim/R/scenario.r -- scenario spec: construct, validate, and persist
#
# OWNS
#   scenario(...)            construct a validated scenario spec list from the
#                            12 DGP keys; errors on invalid combinations
#   scenario_from_file(path) source a scenarios/*.R file and return its
#                            validated spec; entry point for run_one.r
#   write_manifest(spec, dir) write scenario.env to a results directory;
#                              scenario.env is generated provenance, not a
#                              hand-maintained input -- it is always derived
#                              from the scenario spec, never edited by hand
#   read_manifest(path)      parse scenario.env back into a named list for
#                              comparison / audit
#   check_manifest(spec, dir) compare a spec against the existing manifest
#                              in dir; stop with a clear mismatch message if
#                              any of the 12 signature keys differ (mirrors the
#                              bash check in sim/chpc/run_sim.sh)
#
# INPUTS
#   ...      named arguments corresponding to the 12 DGP keys:
#              nobs, zero_gap_target,
#              beta1_active (numeric vector), beta2_active ("neg_beta1" or vector),
#              active_block (character vector, same length as beta1_active),
#              block_props  (named numeric; names: block labels),
#              block_rho    (named numeric; names: block labels),
#              latent_type  ("gaussian"),
#              latent_q     (named numeric; names: block labels),
#              cens_rate, u_min, u_max
#   path     character(1)  path to a scenarios/*.R file
#   spec     list          validated scenario spec (output of scenario())
#   dir      character(1)  directory where scenario.env is written/read
#
# RETURNS
#   scenario():            list of class "bhcrr_scenario" with the 12 DGP keys
#                          plus $blocks (derived block structure for dgp.r)
#   scenario_from_file():  same as scenario()
#   write_manifest():      invisible(path) to the written file
#   read_manifest():       named list of raw key=value strings
#   check_manifest():      invisible(TRUE) on match; stop() on mismatch
#
# DEPENDS ON
#   sim/R/  (none; scenario.r has no intra-sim dependencies)
#   bhCRR:: (none)
#
# REPLACES
#   The scenario manifest logic in sim/chpc/run_sim.sh (bash) lines 96-123
#   The scenario spec is the single definition of the DGP; the scenarios/
#   directory is its only authoritative home
#
# NOT YET IMPLEMENTED
