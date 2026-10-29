plot_model_comparison <- function(df, active_only = TRUE) {

  # Define the custom color palette for the models
  model_colors <- c(
    "SSL" = "#E31A1C",             # Bold Red
    "LASSO_BIC" = "#A6CEE3",       # Light Blue
    "LASSO_Wolbers" = "#1F78B4",   # Dark Blue
    "aLASSO_BIC" = "#B2DF8A",      # Light Green
    "aLASSO_Wolbers" = "#33A02C",  # Dark Green
    "SCAD_BIC" = "#CAB2D6",        # Light Purple
    "SCAD_Wolbers" = "#6A3D9A",     # Dark Purple
    "MCP_BIC" = "#FDBF6F",         # Light Orange
    "MCP_Wolbers" = "#FF7F00"      # Dark Orange
  )

  # 1. Calculate Biases and Summarize
  summary_df <- df %>%
    group_by(Model) %>%
    mutate(
      Active = True != 0,
      Bias_Est = Estimate - True,
      Bias_FG = FG - True,
      Pct_Bias_Est = 100 * (Estimate - True) / True,
      Pct_Bias_FG  = 100 * (FG - True) / True
    ) %>%
    # Filter if requested
    {if (active_only) filter(., Active == TRUE) else .} %>%
    group_by(predictor, Model, True, Active) %>%
    summarise(
      n_runs       = n(),

      # New Model (SSL) Summaries
      Avg_Est      = mean(Estimate, na.rm = TRUE),
      SD_Est       = sd(Estimate, na.rm = TRUE),
      Avg_Bias_Est = mean(Bias_Est, na.rm = TRUE),
      SD_Bias_Est  = sd(Bias_Est, na.rm = TRUE),
      SE_Bias_Est  = sd(Bias_Est, na.rm = TRUE) / sqrt(n()),
      Avg_Pct_Bias_Est = mean(Pct_Bias_Est, na.rm = TRUE),
      SD_Pct_Bias_Est  = sd(Pct_Bias_Est, na.rm = TRUE),
      SE_Pct_Bias_Est  = sd(Pct_Bias_Est, na.rm = TRUE) / sqrt(n()),

      # True FG Summaries
      Avg_FG       = mean(FG, na.rm = TRUE),
      SD_FG        = sd(FG, na.rm = TRUE),
      Avg_Bias_FG  = mean(Bias_FG, na.rm = TRUE),
      SD_Bias_FG   = sd(Bias_FG, na.rm = TRUE),
      SE_Bias_FG   = sd(Bias_FG, na.rm = TRUE) / sqrt(n()),

      .groups = "drop"
    ) %>%
    # Sort predictors numerically (V1, V2, ... V10)
    mutate(predictor = factor(predictor,
                              levels = stringr::str_sort(unique(predictor),
                                                         numeric = TRUE)))

  # PLOT 1: Estimates Comparison
  lvls <- if (is.factor(summary_df$predictor)) {
    levels(summary_df$predictor)
  } else {
    sort(unique(as.character(summary_df$predictor)))
  }

  true_df <- summary_df %>%
    distinct(predictor, True) %>%
    mutate(x = match(as.character(predictor), lvls))

  hw <- 0.45  # half-width of the true-value line; tweak to taste

  p1 <- ggplot(summary_df, aes(x = predictor, y = Avg_Est, color = Model)) +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbar(aes(ymin = Avg_Est -  SD_Est,
                      ymax = Avg_Est +  SD_Est),
                  position = position_dodge(width = 0.5), width = 0.3) +
    # True simulated value: dotted line spanning the whole group
    geom_segment(data = true_df,
                 aes(x = x - hw, xend = x + hw, y = True, yend = True),
                 color = "black", linetype = "dashed", linewidth = 0.8,
                 inherit.aes = FALSE) +
    scale_color_manual(values = model_colors) + # Added Color Mapping
    labs(
      x = "Predictor",
      y = "Average Estimate",
      title = "Estimator Comparison: Fine and Gray with Different Penalizations",
      subtitle = "Dotted line indicates the true coefficient used in simulation.",
      caption = "Error bars: Mean Estimate +/-  SD(Estimate)"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.position = "bottom")



  # PLOT 2: Bias Comparison
  p2 <- ggplot(summary_df, aes(x = predictor, y = Avg_Bias_Est, color = Model)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.8) +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbar(aes(ymin = Avg_Bias_Est - 1.96*SE_Bias_Est,
                      ymax = Avg_Bias_Est + 1.96*SE_Bias_Est),
                  position = position_dodge(width = 0.5), width = 0.3) +
    scale_color_manual(values = model_colors) + # Added Color Mapping
    labs(
      x = "Predictor",
      y = "Average Bias (Estimate - True Beta)",
      title = "Bias Comparison: Fine and Gray with Different Penalizations",
      caption = "Error bars: Mean Bias +/- 1.96* Standard Error of the Bias"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.position = "bottom")

  # PLOT 3: Percent Bias Comparison
  p3 <- ggplot(summary_df, aes(x = predictor, y = Avg_Pct_Bias_Est, color = Model)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.8) +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbar(aes(ymin = Avg_Pct_Bias_Est - 1.96*SE_Pct_Bias_Est,
                      ymax = Avg_Pct_Bias_Est + 1.96*SE_Pct_Bias_Est),
                  position = position_dodge(width = 0.5), width = 0.3) +
    scale_color_manual(values = model_colors) + # Added Color Mapping
    labs(
      x = "Predictor",
      y = "Average Percent Bias (%)  [(Estimate - True Beta) / True Beta]",
      title = "Percent Bias Comparison: Fine and Gray with Different Penalizations",
      caption = "Error bars: Mean Percent Bias +/- 1.96* Standard Error of the Percent Bias"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.position = "bottom")

  # Print plots
  print(p1)
  print(p2)
  print(p3)

  # Invisibly return the summarized data in case you need to inspect it
  invisible(summary_df)
}
