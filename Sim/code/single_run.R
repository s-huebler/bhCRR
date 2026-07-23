#################### Set Up ################
library(tidyverse)
library(fastcmprsk, lib.loc = "/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library")
library(survival)

 source("/Users/sophiehuebler/Documents/bhCRR/R/helpers.r")
 source("/Users/sophiehuebler/Documents/bhCRR/R/utils.r")
 source("/Users/sophiehuebler/Documents/bhCRR/R/dedupe_warnings.r")

 source("/Users/sophiehuebler/Documents/bhCRR/R/update_betas.r")
 source("/Users/sophiehuebler/Documents/bhCRR/R/expected_inclusion_probs.r")
 source("/Users/sophiehuebler/Documents/bhCRR/R/expected_penalty_weights.r")
 source("/Users/sophiehuebler/Documents/bhCRR/R/update_mixture_prob.r")

 Rcpp::sourceCpp("/Users/sophiehuebler/Documents/bhCRR/src/cv_fastcrrp.cpp")
 Rcpp::sourceCpp("/Users/sophiehuebler/Documents/bhCRR/src/RcppExports.cpp")
 source("/Users/sophiehuebler/Documents/bhCRR/R/cv_fastCrrp_cpp.r")
 source("/Users/sophiehuebler/Documents/bhCRR/R/fit_ssl_psdh.r")

 source("/Users/sophiehuebler/Documents/bhCRR/R/generate_foldid.r")
 source("/Users/sophiehuebler/Documents/bhCRR/R/predict_from_ssl_psdh.r")
 source("/Users/sophiehuebler/Documents/bhCRR/R/wolbers_c.r")
 source("/Users/sophiehuebler/Documents/bhCRR/R/cv_ssl_psdh.r")
 source("/Users/sophiehuebler/Documents/bhCRR/R/tune_ssl_psdh.r")
 source("/Users/sophiehuebler/Documents/bhCRR/R/threshold.R")
 source("/Users/sophiehuebler/Documents/bhCRR/R/tuning_diagnostics.r")



#################### Data Generation ################

# To change per scenario
nobs = 200
npredictors = 25
beta1_active = c(0.40, - 0.50, 0.60, 0.75, - 0.80)
beta2_active = c(0, 0.3, 0, 0, -0.2)


# Calculated (fixed)
beta1 <- c(beta1_active, rep(0, npredictors-length(beta1_active)))
beta2 <- c(beta2_active, rep(0, npredictors-length(beta2_active)))

# Simulated
sim <- simulateTwoCauseFineGrayModel(nobs = nobs,
                                     beta1 = beta1,
                                     beta2 = beta2,
                                     X = NULL,
                                     u.min = 100,
                                     u.max = 100,
                                     p = 0.5,
                                     returnX = TRUE)

sim_data <- cbind(data.frame(ID = 1:nobs,
                             TTE = sim$ftime,
                             Status = as.integer(sim$fstatus)),
                  sim$X)

names(sim_data)<- c("ID", "TTE", "Status", paste0("X_", 1:(ncol(sim_data)-3)))

x <- as.matrix(sim_data %>%select(starts_with("X")))
y <- as.matrix(sim_data %>%
                 mutate(Status = as.numeric(Status))%>%
                 select(TTE, Status))




#################### SSL Model ################

ssl_psdh_time1 <- Sys.time()

#Initial model fit
mod <- fit_ssl_psdh(x, y,
                    ss=c(0.04, 0.6),
                    initial_sparsity = 0.05,
                    theta_a = 1,
                    theta_b = 1,
                    maxit = 50,
                    epsilon=1e-04,
                    init = NULL,
                    init_lam_path = 10^seq(log10(0.1),
                                                log10(0.001),
                                                length = 10),
                    inner_maxit_start = 1000)

#Tuning
#
# The manual pre-flight / post-flight checks are now automated in
# R/tuning_diagnostics.r (bhcrr_autotune). It: estimates theta from the untuned
# fit, clamps the reasonable s1 range into the feasible band, solves for the s0
# region around the clinical zero-gap target, recommends and validates a grid,
# runs tune_ssl_psdh(), and auto-widens + re-tunes once if the optimum lands on
# a grid edge. Everything runs unattended and reports flags for review.

zero_gap_target <- 0.1   # clinically relevant minimum treatment effect

autotune <- bhcrr_autotune(mod,
                           beta_min   = zero_gap_target,
                           beta_floor = 0.01,
                           nfolds     = 10,
                           ncv        = 2,
                           foldid     = NULL,
                           max_widen  = 1)

# Pre-flight diagnostics, the validation checks, and the chosen pair.
print(autotune$preflight)
print(autotune)
print(autotune$validation$checks)
plot_score_heatmap(autotune$tuning)

# Tuned scale parameters (spike s0, slab s1) selected by the autotuner.
best_s0 <- autotune$best$s0
best_s1 <- autotune$best$s1



final_mod <-fit_ssl_psdh(x, y,
                         ss=c(best_s0, best_s1),
                         initial_sparsity = 0.05,
                         theta_a = 1,
                         theta_b = 1,
                         maxit = 100,
                         epsilon=1e-04,
                         init = NULL,
                         init_lam_path = 10^seq(log10(0.1),
                                                log10(0.001),
                                                length = 10),
                         inner_maxit_start = 1000)



ssl_psdh_time2 <- Sys.time()
ssl_psdh_time = ssl_psdh_time2 - ssl_psdh_time1



