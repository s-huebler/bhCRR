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

#Find reasonable scale area for s0 based on what is reasonable/possible for s1 and what the clinically relevant treatment effects might be

#Fix theta as the estimated proportion of non-zero coefficients in the untuned model
theta_guess <- sum(mod$coefficients$Estimate != 0)/nrow(mod$coefficients)

#Find bounds for s1, make sure that reasonable range (0.5 - 1) falls within the bounds
s1_bounds <- s1_range(nobs, npredictors)
if(s1_bounds[1] > 0.5 | 1 > s1_bounds[2]){
  s1_guess1 <- s1_bounds[2]
  s1_guess2 <- s1_bounds[2]
}else{
  s1_guess1 <- 0.5
  s1_guess2 <- 1
}

#Supply clinically relevant minimum treatment effect (this value will set where the zero gap is)

zero_gap_target <- 0.1

#Make sure the tuning covers around this area where I solve for the upper bound of s0 across a couple of guesses of s1 and maybe a couple potential zero gap targets
s0_upper_bound1 <- s0_upper_range(nobs, s1_guess1)
solve_for_s0(beta_min = zero_gap_target, n = nobs, s1 = s1_guess1, theta = theta_guess, lower_bound = 0.0001, upper_bound = s0_upper_bound1)

s0_upper_bound2 <- s0_upper_range(nobs, s1_guess2)
solve_for_s0(beta_min = zero_gap_target, n = nobs, s1 = s1_guess2, theta = theta_guess, lower_bound = 0.0001, upper_bound = s0_upper_bound2)

#Absolute smallest s0 could be to have a zero gap, just checking that we aren't under that
solve_for_s0(beta_min = 0.01, n = nobs, s1 = s1_guess2, theta = theta_guess, lower_bound = 0.0001, upper_bound = s0_upper_bound2)




# Run the tuning
tuning1 <- tune_ssl_psdh(mod, seq(0.005, 0.034, 0.004),
                        c(0.5, 0.75, 1), nfolds=10, ncv=2, foldid=NULL)
plot_score_heatmap(tuning1)

# Zoom in a bit
tuning2 <-  tune_ssl_psdh(mod, seq(0.015, 0.025, 0.002),
                          c(0.45, 0.5, 0.55), nfolds=10, ncv=2, foldid=NULL)
plot_score_heatmap(tuning2)

# Make sure there is only 1 max value
tuning2 %>%
  filter(score_mean == max(tuning2$score_mean))

# Check that the zero gap is still reasonable
zero_gap(200, 0.023, 0.45, .2)



final_mod <-fit_ssl_psdh(x, y,
                         ss=c(0.023, 0.45),
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



