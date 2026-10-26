# sim/R/load.r -- bring bhCRR into scope for every sim/ script
#
# OWNS
#   load_bhcrr(mode)   switch between local dev (devtools::load_all) and
#                      cluster execution (library(bhCRR)), selected by the
#                      'mode' argument or the BHCRR_LOAD env var
#
# INPUTS
#   mode  character(1)  "dev" or "cluster"; defaults to Sys.getenv("BHCRR_LOAD",
#                       "dev"); "dev" calls devtools::load_all(repo_root()),
#                       "cluster" calls library(bhCRR)
#
# RETURNS
#   invisible(NULL); called for its side effect of attaching bhCRR
#
# DEPENDS ON
#   (none within sim/R; devtools for "dev" mode, bhCRR installed for "cluster")
#
# REPLACES
#   The source() calls scattered across sim/_legacy/code/*.R (every legacy
#   script sourced individual R/ files by relative path; those paths broke
#   whenever the working directory changed). This file is the ONE place that
#   loads bhCRR so no other sim/ file ever contains a source() path.
#
# NOT YET IMPLEMENTED
