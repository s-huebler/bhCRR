library(tidyverse)
library(fastcmprsk, lib.loc = "/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library")
library(survival)

source("/Users/sophiehuebler/Documents/bhCRR/R/cv_ssl_psdh.r")
source("/Users/sophiehuebler/Documents/bhCRR/R/fit_ssl_psdh.r")
source("/Users/sophiehuebler/Documents/bhCRR/R/update_betas.r")
source("/Users/sophiehuebler/Documents/bhCRR/R/expected_inclusion_probs.r")
source("/Users/sophiehuebler/Documents/bhCRR/R/expected_penalty_weights.r")
source("/Users/sophiehuebler/Documents/bhCRR/R/update_mixture_prob.r")
source("/Users/sophiehuebler/Documents/bhCRR/R/helpers.r")
source("/Users/sophiehuebler/Documents/bhCRR/R/tune_ssl_psdh.r")
source("/Users/sophiehuebler/Documents/bhCRR/R/wolbers_c.r")
source("/Users/sophiehuebler/Documents/bhCRR/R/predict_from_ssl_psdh.r")
source("/Users/sophiehuebler/Documents/bhCRR/R/generate_foldid.r")
source("/Users/sophiehuebler/Documents/bhCRR/R/dedupe_warnings.r")
source("/Users/sophiehuebler/Documents/bhCRR/R/cv_fastCrrp.r")


###########Generate Data Functions ##########################################

easy_sim <- function(nobs, beta1, beta2){
  sim <- simulateTwoCauseFineGrayModel(nobs = 50,
                                       beta1 = beta1,
                                       beta2 = beta2,
                                       X = NULL,
                                       u.min = 100,
                                       u.max = 100,
                                       p = 0.5,
                                       returnX = TRUE)

  sim_data <- cbind(data.frame(ID = nobs,
                               TTE = sim$ftime,
                               Status = as.integer(sim$fstatus)),
                    sim$X)

  names(sim_data)<- c("ID", "TTE", "Status", paste0("X_", 1:(ncol(sim_data)-3)))

  x <- as.matrix(sim_data %>%select(starts_with("X")))
  y <- as.matrix(sim_data %>%
                   mutate(Status = as.numeric(Status))%>%
                   select(TTE, Status))

  return(list("dataframe" = sim_data,
              "x" = x,
              "y" = y))
}

easy_sim_wrapper <- function(x){
data <- easy_sim(50,
                     c(0.40, - 0.50, 0.60, 0.75, - 0.80, rep(0,x-5)),
                     c(0, 0.3, 0, 0, -0.2, rep(0,x-5))
)

comparisons <- data.frame(Name = c(paste("B", 1:x, sep = "")),
                   True_beta = c(c(0.40, - 0.50, 0.60, 0.75, - 0.80),
                                 rep(0, x-5))
)

return(list(data = data,
            comparisons=comparisons))
}




#####################Outcomes Function####################################
calculate_metrics <- function(predicted, true_val) {

  # Replace NA values with 0 (assuming NA means the variable was not selected)
  predicted <- ifelse(is.na(predicted), 0, predicted)
  true_val <- ifelse(is.na(true_val), 0, true_val)

  # Determine selection status (non-zero coefficient means selected)
  pred_selected <- predicted != 0
  true_selected <- true_val != 0

  # Calculate Confusion Matrix elements
  TP <- sum(pred_selected & true_selected)
  TN <- sum(!pred_selected & !true_selected)
  FP <- sum(pred_selected & !true_selected)
  FN <- sum(!pred_selected & true_selected)

  # Calculate Metrics (with safeguards against division by zero)
  sensitivity <- ifelse((TP + FN) == 0, NA, TP / (TP + FN))
  specificity <- ifelse((TN + FP) == 0, NA, TN / (TN + FP))
  ppv <- ifelse((TP + FP) == 0, NA, TP / (TP + FP))

  # F1 Score: 2 * (Precision * Recall) / (Precision + Recall)
  # Note: PPV is Precision, Sensitivity is Recall
  f1 <- ifelse((ppv + sensitivity) == 0 | is.na(ppv) | is.na(sensitivity),
               NA,
               2 * (ppv * sensitivity) / (ppv + sensitivity))

  # Mean Squared Error
  mse <- mean((predicted - true_val)^2)

  # Return as a single-row dataframe
  return(data.frame(
    TP = TP,
    TN = TN,
    FP = FP,
    FN = FN,
    Sensitivity = sensitivity,
    Specificity = specificity,
    PPV = ppv,
    F1 = f1,
    MSE = mse
  ))
}

# Identify the columns that contain model predictions
# We exclude 'Name' and 'True_beta'
metrics_summary <- function(df){
  model_cols <- setdiff(names(df), c("Name", "True_beta"))

  # Apply the function to each model column and bind the results
  results_list <- lapply(model_cols, function(model) {
    metrics <- calculate_metrics(df[[model]], df$True_beta)
    metrics$Model <- model
    return(metrics)
  })

  # Combine everything into a clean summary table
  results_df <- bind_rows(results_list) %>%
    select(TP, TN, FP, FN, Model, Sensitivity, Specificity, PPV, F1, MSE) # Reorder columns for readability

  results_df
}


#####################SSL Wrap Function####################################

psdh_ssl_func <- function(x, y, init_sparsity, ncv_count,
                          tuning = "wolbers",
                          s0_seq = seq(0.005, 0.1, length.out = 5),
                          s1_seq = seq(0.5,   0.9, length.out = 5),
                          nfolds = 10,
                          initial_shrinkage = NULL,
                          fixed_global_shrinkage = NULL,
                          save_tunes = NULL) {

  start_time <- Sys.time()
  x <- as.matrix(x)
  y <- as.matrix(y)
  fit <- fit_ssl_psdh(x, y, initial_sparsity = init_sparsity, fixed_global_shrinkage = fixed_global_shrinkage, initial_shrinkage = NULL)

  if (tuning == "wolbers") {
    tunes <- tryCatch(
      tune_ssl_psdh(fit,
                    s0_seq,
                    s1_seq,
                    ncv              = ncv_count,
                    initial_sparsity = init_sparsity),
      error = function(msg) NA
    )
  } else if (tuning == "normal") {
    # Manual grid search using standard Harrell C-index (survival::concordance)
    # on the SSL linear predictor x_test %*% beta. Folds shared across all
    # (s0, s1) pairs and pooled across folds before averaging over reps.
    n      <- nrow(y)
    fol    <- generate_foldid(nobs = n, nfolds = nfolds,
                              foldid = NULL, ncv = ncv_count)
    foldid <- fol$foldid

    param_grid <- expand.grid(s0 = s0_seq, s1 = s1_seq)
    valid_grid <- param_grid[param_grid$s1 > param_grid$s0, ]
    if (nrow(valid_grid) == 0) {
      stop("No valid (s0, s1) pairs with s1 > s0.")
    }

    results_list <- lapply(seq_len(nrow(valid_grid)), function(p) {
      s0_p <- valid_grid$s0[p]
      s1_p <- valid_grid$s1[p]

      pooled <- rep(NA_real_, ncv_count)
      for (rep_idx in seq_len(ncv_count)) {
        lp_all <- rep(NA_real_, n)
        for (f in seq_len(nfolds)) {
          omit_idx <- which(foldid[, rep_idx] == f)
          y_train  <- as.matrix(y[-omit_idx, ])
          x_train  <- as.matrix(x[-omit_idx, , drop = FALSE])

          sub_fit <- try(suppressWarnings(
            fit_ssl_psdh(x = x_train, y = y_train,
                         ss = c(s0_p, s1_p),
                         initial_sparsity = init_sparsity,
                         maxit = 50, epsilon = 1e-04,
                         fixed_global_shrinkage = fixed_global_shrinkage)
          ), silent = TRUE)
          if (inherits(sub_fit, "try-error")) next

          x_test <- x[omit_idx, , drop = FALSE]
          lp_all[omit_idx] <- as.vector(
            x_test %*% as.vector(sub_fit$final_model_object$coef)
          )
        }

        pooled[rep_idx] <- tryCatch(
          survival::concordance(
            Surv(y[, 1], y[, 2] == 1) ~ lp_all
          )$concordance,
          error = function(e) NA_real_
        )
      }

      c(s0         = s0_p,
        s1         = s1_p,
        score_mean = mean(pooled, na.rm = TRUE),
        score_sd   = if (ncv_count == 1) NA_real_ else sd(pooled, na.rm = TRUE))
    })

    tunes <- as.data.frame(do.call(rbind, results_list))
  } else {
    stop("`tuning` must be 'normal' or 'wolbers'.")
  }

  #print(head(tunes))
  max_tune <- tunes %>%
    dplyr::filter(score_mean == max(tunes$score_mean, na.rm = TRUE))

  if(!is.null(save_tunes)){
    assign(x = save_tunes, value = tunes, envir = .GlobalEnv)
  }


  final_fit <- fit_ssl_psdh(x, y,
                            ss = c(max_tune$s0[1], max_tune$s1[1]),
                            initial_sparsity = init_sparsity,
                            fixed_global_shrinkage = fixed_global_shrinkage)

  end_time <-Sys.time()
  print(end_time - start_time)

  final_fit
}


#####################Run Generate Date####################################
one_sim <- function(x){
generate_data <- easy_sim_wrapper(x) #Replace with argument

sim_data <- generate_data$data

sim_comparisons <- generate_data$comparison

#####################Run Functions####################################

####################FINE AND GRAY ##########################################################
sim_data$dataframe$Status_factor <- factor(sim_data$dataframe$Status,
                                      levels = c(0, 1, 2),
                                      labels = c("censor", "event", "competing"))




sim_fg_data <- finegray(Surv(TTE, Status_factor) ~ ., data = sim_data$dataframe, etype = "event")




sim_fg_model <- coxph(Surv(fgstart, fgstop, fgstatus) ~ X_1 + X_2 + X_3 + X_4 + X_5,
                      weight = fgwt,
                      data = sim_fg_data)
sim_comparisons$Fg = c(sim_fg_model$coefficients, rep(0, x-5))

####################FASTCMPRSK

####################NORMAL
sim_fstcmp_normalTuned <- cv_fastCrrp(x = sim_data$x, time = sim_data$y[,1], status = sim_data$y[,2],
                                      k = 5, penalty = "LASSO",
                                      lambda = 10^seq(log10(0.1),log10(0.001),length = 25),
                                      tuning = "normal")

sim_fstcmp_normalTuned_best_lambda <- sim_fstcmp_normalTuned$lambda_min
sim_comparisons$fstcmp_normalTuned = as.vector(sim_fstcmp_normalTuned$full_model$coef[,
                                                                                      which(sim_fstcmp_normalTuned$lambda == sim_fstcmp_normalTuned_best_lambda)])


####################WOLBERS
sim_fstcmp_wolbersTuned <- cv_fastCrrp(x = sim_data$x, time = sim_data$y[,1], status = sim_data$y[,2],
                                       k = 5, penalty = "LASSO",
                                       lambda = 10^seq(log10(0.1),log10(0.001),length = 25),
                                       tuning = "wolbers")
sim_fstcmp_wolbersTuned_best_lambda <- sim_fstcmp_wolbersTuned$lambda_min
sim_comparisons$fstcmp_wolbersTuned = as.vector(sim_fstcmp_wolbersTuned$full_model$coef[,
                                                                                        which(sim_fstcmp_wolbersTuned$lambda == sim_fstcmp_wolbersTuned_best_lambda)])


####################SSL ####################################################


####################NORMAL
sim_ssl_normalTuned <- psdh_ssl_func( sim_data$x, sim_data$y, 0.5, 1,
                                      tuning = "normal",
                                      s0_seq = seq(0.005, 0.1, length.out = 10),
                                      s1_seq = seq(0.8,   1, length.out = 4),
                                      nfolds = 10,
                                      fixed_global_shrinkage = NULL,
                                      save_tunes = NULL)
sim_comparisons$ssl_normalTuned = as.vector(sim_ssl_normalTuned$final_model_object$coef)

####################WOLBERS
sim_ssl_wolbersTuned <- psdh_ssl_func( sim_data$x, sim_data$y, 0.5, 1,
                                       tuning = "wolbers",
                                       s0_seq = seq(0.005, 0.1, length.out =10),
                                       s1_seq = seq(0.8,   1, length.out = 4 ),
                                       nfolds = 10,
                                       fixed_global_shrinkage = NULL,
                                       save_tunes = NULL)
sim_comparisons$ssl_wolbersTuned = as.vector(sim_ssl_wolbersTuned$final_model_object$coef)

#####################Outcomes####################################

model_results <- metrics_summary(sim_comparisons)
model_results$Predictors <- x #Replace with arguement
#print(model_results)

sim_comparisons
}

for(i in 1:15){
  temp_file <- paste0("~/Documents/bhCRR/Sim/Local_Testing/Sim100_", i, ".Rdata")
  temp <- one_sim(100)
  save(temp, file = temp_file)
}

