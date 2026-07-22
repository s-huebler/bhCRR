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

#Initial model fit
mod <- fit_ssl_psdh(x, y,
                    ss=c(0.04, 0.7),
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

#Check each part of tuning is working
fol <- generate_foldid(nobs=nobs, nfolds=10, foldid=NULL, ncv=2)

eval_time <- quantile(y[y[,2]==1,1], .5)

predicted <- predict_from_ssl_psdh(mod, x, eval_time)

cindex <- wolbers_c(y, predicted, eval_time)

cv_ssl_psdh(mod, fol$foldid, 0.04, 0.6, 2, .5)

test <- matrix(NA,
               nrow = length(seq(0.03, 0.08, 0.01)),
               ncol = length(seq(0.6, 1, .15)))

loops_time1 <- Sys.time()
for(i in seq(0.03, 0.08, 0.01)){
  for(j in seq(0.6, 1, .15)){
    test[i,j]<- cv_ssl_psdh(mod, fol$foldid, i, j, 2, .5)$measures[1]
  }
}
loops_time2 <- Sys.time()
print(loops_time2-loops_time1)


tune_time1 <- Sys.time()
tuning <- tune_ssl_psdh(mod, seq(0.03, 0.08, 0.01), seq(0.6, 1, .15), nfolds=10, ncv=2, foldid=NULL)
tune_time2 <- Sys.time()
print(tune_time2-tune_time1)
