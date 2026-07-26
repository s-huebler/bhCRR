plot_model_comparison <- function(df, active_only = TRUE) {

  # 1. Calculate Biases and Summarize
  summary_df <- df %>%
    mutate(
      Active = True_beta != 0,
      Bias_Est = Estimate - True_beta,
      Bias_FG = FG_True - True_beta
    ) %>%
    # Filter if requested
    {if (active_only) filter(., Active == TRUE) else .} %>%
    group_by(predictor, True_beta, Active) %>%
    summarise(
      n_runs       = n(),

      # New Model (SSL) Summaries
      Avg_Est      = mean(Estimate, na.rm = TRUE),
      SD_Est       = sd(Estimate, na.rm = TRUE),
      Avg_Bias_Est = mean(Bias_Est, na.rm = TRUE),
      SE_Bias_Est  = sd(Bias_Est, na.rm = TRUE) / sqrt(n()),

      # True FG Summaries
      Avg_FG       = mean(FG_True, na.rm = TRUE),
      SD_FG        = sd(FG_True, na.rm = TRUE),
      Avg_Bias_FG  = mean(Bias_FG, na.rm = TRUE),
      SE_Bias_FG   = sd(Bias_FG, na.rm = TRUE) / sqrt(n()),

      .groups = "drop"
    ) %>%
    # Sort predictors numerically (V1, V2, ... V10)
    mutate(predictor = factor(predictor, levels = stringr::str_sort(unique(predictor), numeric = TRUE)))

  # 2. Pivot for Plot 1: Average Estimates
  # We want New Model and FG Model side-by-side
  est_long <- summary_df %>%
    select(predictor, True_beta, Active, Avg_Est, SD_Est, Avg_FG, SD_FG) %>%
    pivot_longer(
      cols = c(Avg_Est, Avg_FG),
      names_to = "Model",
      values_to = "Mean_Estimate"
    ) %>%
    mutate(
      SD = ifelse(Model == "Avg_Est", SD_Est, SD_FG),
      Model = ifelse(Model == "Avg_Est", "SSL Model", "True FG")
    )

  # PLOT 1: Estimates Comparison
  p1 <- ggplot(est_long, aes(x = predictor, y = Mean_Estimate, color = Model)) +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbar(aes(ymin = Mean_Estimate - 1.96 * SD,
                      ymax = Mean_Estimate + 1.96 * SD),
                  position = position_dodge(width = 0.5), width = 0.3) +
    # Add True Beta marker dead center (no dodging)
    geom_point(aes(x = predictor, y = True_beta),
               color = "black", shape = 15, size = 2.5, inherit.aes = FALSE) +
    labs(
      x = "Predictor",
      y = "Average Estimate",
      title = "Estimator Comparison: SSL vs. True Fine-Gray",
      subtitle = "Black squares indicate the true simulated parameter value.",
      caption = "Error bars: Mean Estimate +/- 1.96 * SD(Estimate)"
    ) +
    scale_color_manual(values = c("SSL Model" = "#E69F00", "True FG" = "#56B4E9")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.position = "bottom")


  # 3. Pivot for Plot 2: Average Bias
  bias_long <- summary_df %>%
    select(predictor, True_beta, Active, Avg_Bias_Est, SE_Bias_Est, Avg_Bias_FG, SE_Bias_FG) %>%
    pivot_longer(
      cols = c(Avg_Bias_Est, Avg_Bias_FG),
      names_to = "Model",
      values_to = "Mean_Bias"
    ) %>%
    mutate(
      SE = ifelse(Model == "Avg_Bias_Est", SE_Bias_Est, SE_Bias_FG),
      Model = ifelse(Model == "Avg_Bias_Est", "SSL Model", "True FG")
    )

  # PLOT 2: Bias Comparison
  p2 <- ggplot(bias_long, aes(x = predictor, y = Mean_Bias, color = Model)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.8) +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbar(aes(ymin = Mean_Bias - SE,
                      ymax = Mean_Bias + SE),
                  position = position_dodge(width = 0.5), width = 0.3) +
    labs(
      x = "Predictor",
      y = "Average Bias (Estimate - True Beta)",
      title = "Bias Comparison: SSL vs. True Fine-Gray",
      caption = "Error bars: Mean Bias +/- 1 Standard Error of the Bias"
    ) +
    scale_color_manual(values = c("SSL Model" = "#E69F00", "True FG" = "#56B4E9")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.position = "bottom")

  # Print plots
  print(p1)
  print(p2)

  # Invisibly return the summarized data in case you need to inspect it
  invisible(summary_df)
}


plot_model_comparison(betas_p25, active_only = TRUE)

plot_model_comparison(betas_p5, active_only = TRUE)
