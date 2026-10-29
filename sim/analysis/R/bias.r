# sim/analysis/R/bias.r -- coefficient estimates, bias and percent bias
#
# OWNS
#   summarise_bias(betas, models, active_only = TRUE)
#                       Monte Carlo summary per predictor x model: mean and SD
#                       of the estimate, mean bias with its Monte Carlo SE,
#                       mean percent bias with its SE
#   bias_by_block(betas, models)
#                       per model x block, over the active predictors: mean
#                       percent bias, mean absolute bias, RMSE, share of
#                       estimates exactly zero. Raw mean bias is left out on
#                       purpose: with mixed-sign truths, shrinkage of positive
#                       and negative effects cancels in it.
#   plot_estimates(bias_sum, models)
#   plot_bias(bias_sum, models)
#   plot_pct_bias(bias_sum, models)
#                       one panel per correlation block, predictors on x,
#                       one dodged point + interval per model
#
# REPLACES
#   plot_model_comparison() in the old plot_metrics.R / plot_model_comparisons.R,
#   which hard-coded the model palette and printed three plots from one call.
#   The summary and each plot are now separate so the .qmd can place, caption
#   and size them independently.
#
# INPUT CONTRACT (betas table from parse.r)
#   predictor, Block, BlockRho, True, Model, Estimate, run_number
#   Failed fits arrive with Estimate == NA for every predictor. They are
#   dropped here with a message, never averaged in as zeros.
#
# READING PERCENT BIAS
#   100 * (Estimate - True) / True. Because it is scaled by the signed truth,
#   a negative value means the estimate is shrunk toward zero whatever the sign
#   of the true coefficient, and a positive value means it is inflated away
#   from zero. Raw bias does not have this property: shrinkage of a negative
#   coefficient shows up as positive bias.

summarise_bias <- function(betas, models, active_only = TRUE) {
  models <- resolve_models(models, betas$Model)
  d <- betas[betas$Model %in% models$model, ]
  if (active_only) d <- d[d$True != 0, ]

  failed <- is.na(d$Estimate)
  if (any(failed)) {
    f <- unique(d[failed, c("Model", "run_number")])
    msg <- vapply(split(f$run_number, f$Model), length, integer(1))
    message("Dropping failed fits (NA estimates), runs per model: ",
            paste(sprintf("%s = %d", names(msg), msg), collapse = ", "))
    d <- d[!failed, ]
  }

  d$Bias     <- d$Estimate - d$True
  d$Pct_Bias <- ifelse(d$True != 0, 100 * d$Bias / d$True, NA_real_)

  out <- d |>
    dplyr::group_by(predictor, Block, BlockRho, True, Model) |>
    dplyr::summarise(
      n_runs      = dplyr::n(),
      Avg_Est     = mean(Estimate),
      SD_Est      = stats::sd(Estimate),
      Avg_Bias    = mean(Bias),
      SE_Bias     = stats::sd(Bias) / sqrt(dplyr::n()),
      Avg_Pct_Bias = mean(Pct_Bias),
      SE_Pct_Bias = stats::sd(Pct_Bias) / sqrt(dplyr::n()),
      Prop_Zero   = mean(Estimate == 0),
      .groups = "drop") |>
    as.data.frame()

  out <- label_models(out, models)
  out$predictor <- factor(out$predictor,
                          levels = stringr::str_sort(unique(out$predictor),
                                                     numeric = TRUE))
  out$Block <- .block_factor(out$Block, out$BlockRho)
  attr(out, "models") <- models
  out
}

bias_by_block <- function(betas, models) {
  models <- resolve_models(models, betas$Model)
  d <- betas[betas$Model %in% models$model & betas$True != 0 &
               !is.na(betas$Estimate), ]
  d$Bias <- d$Estimate - d$True
  out <- d |>
    dplyr::group_by(Model, Block, BlockRho) |>
    dplyr::summarise(
      `Mean % bias`       = mean(100 * Bias / True),
      `Mean |bias|`       = mean(abs(Bias)),
      RMSE                = sqrt(mean(Bias^2)),
      `% estimated as 0`  = 100 * mean(Estimate == 0),
      .groups = "drop") |>
    as.data.frame()
  out <- label_models(out, models)
  out$Block <- .block_factor(out$Block, out$BlockRho)
  out <- out[order(out$Block, out$Model), setdiff(names(out), "BlockRho")]
  rownames(out) <- NULL
  out
}

# Blocks labelled with their rho and ordered by rho, e.g. "strong (rho = 0.7)"
.block_factor <- function(block, rho) {
  key <- unique(data.frame(block = as.character(block), rho = rho))
  key <- key[order(key$rho), ]
  lab <- sprintf("%s (rho = %s)", key$block, as.character(key$rho))
  factor(lab[match(as.character(block), key$block)], levels = lab)
}

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

.bias_base <- function(s, y, ymin, ymax, models, ylab, title, caption,
                       hline = TRUE) {
  dodge <- ggplot2::position_dodge(width = 0.7)
  p <- ggplot2::ggplot(s, ggplot2::aes(x = predictor, y = {{ y }},
                                       colour = Model))
  if (hline)
    p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                                 colour = "grey50")
  p +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = {{ ymin }}, ymax = {{ ymax }}),
                           position = dodge, width = 0) +
    ggplot2::geom_point(position = dodge, size = 1.6) +
    ggplot2::facet_grid(. ~ Block, scales = "free_x", space = "free_x") +
    scale_colour_models(models) +
    ggplot2::labs(x = "Active predictor", y = ylab, title = title,
                  caption = caption, colour = NULL) +
    .theme_sim()
}

plot_estimates <- function(bias_sum, models = attr(bias_sum, "models")) {
  truth <- unique(bias_sum[, c("predictor", "Block", "True")])
  .bias_base(bias_sum, Avg_Est, Avg_Est - SD_Est, Avg_Est + SD_Est, models,
             ylab = "Mean estimate",
             title = "Coefficient estimates for the active predictors",
             caption = "Points: mean over runs. Bars: +/- 1 SD over runs. Black dash: true value.",
             hline = TRUE) +
    ggplot2::geom_errorbar(data = truth,
                           ggplot2::aes(x = predictor, ymin = True, ymax = True),
                           inherit.aes = FALSE, width = 0.9, linewidth = 0.6,
                           colour = "black")
}

plot_bias <- function(bias_sum, models = attr(bias_sum, "models")) {
  .bias_base(bias_sum, Avg_Bias, Avg_Bias - 1.96 * SE_Bias,
             Avg_Bias + 1.96 * SE_Bias, models,
             ylab = "Mean bias (estimate - true)",
             title = "Bias for the active predictors",
             caption = "Bars: +/- 1.96 Monte Carlo SE of the mean bias.")
}

plot_pct_bias <- function(bias_sum, models = attr(bias_sum, "models")) {
  .bias_base(bias_sum, Avg_Pct_Bias, Avg_Pct_Bias - 1.96 * SE_Pct_Bias,
             Avg_Pct_Bias + 1.96 * SE_Pct_Bias, models,
             ylab = "Mean percent bias (%)",
             title = "Percent bias for the active predictors",
             caption = paste("Percent bias = 100 x (estimate - true) / true;",
                             "negative = shrunk toward zero.",
                             "Bars: +/- 1.96 Monte Carlo SE.")) +
    ggplot2::geom_hline(yintercept = -100, linetype = "dotted",
                        colour = "grey60")
}

.theme_sim <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(angle = 90, vjust = 0.5,
                                               hjust = 1, size = 7),
      legend.position  = "bottom",
      panel.grid.major.x = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(face = "bold"),
      plot.caption     = ggplot2::element_text(colour = "grey40", hjust = 0))
}
