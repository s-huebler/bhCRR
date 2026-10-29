# sim/analysis/R/selection.r -- variable selection (confusion matrix) metrics
#
# OWNS
#   selection_metrics(performance, models)
#                       per model x run: tp, fp, tn, fn and the derived rates
#                       (precision, recall, specificity, F1, MCC) plus the
#                       Wolbers C-index; failed fits removed
#   summarise_selection(sel)
#                       per model: mean (SD) of every metric over runs
#   mean_confusion(sel) per model: the 2 x 2 confusion matrix averaged over runs
#   plot_confusion(sel) one tile panel per model showing mean_confusion()
#   plot_selection(sel, metrics = ...)
#                       distribution over runs of the chosen metrics
#   selection_by_block(by_block, performance, models)
#                       per model x block x run: recall and false positive rate
#   plot_selection_by_block(sbb)
#
# SELECTION RULE
#   A predictor is "selected" when its coefficient is nonzero (score.r uses
#   |coef| > 0). Truth is beta1 != 0. These counts are precomputed by score.r;
#   this file only derives rates from them.
#
# FAILED FITS
#   A model-run that did not converge arrives with tp = fp = tn = fn = 0 and
#   NA rates. It is dropped (see .usable_runs() in scenario_summary.r) so that
#   it does not drag the averages toward zero. run_health() reports how many.
#
# RATE CONVENTIONS
#   precision = tp / (tp + fp)   NA when nothing was selected
#   recall    = tp / (tp + fn)   true positive rate
#   F1        = 2 tp / (2 tp + fp + fn)
#   MCC       = (tp tn - fp fn) / sqrt((tp+fp)(tp+fn)(tn+fp)(tn+fn)), NA on a
#               zero denominator. MCC is the single number least distorted by
#               the heavy class imbalance at large p (30 actives vs 970 nulls),
#               where specificity sits near 1 for every method.

.metric_labels <- c(
  tp = "True positives", fp = "False positives",
  fn = "False negatives", tn = "True negatives",
  n_selected = "Selected",
  precision = "Precision", recall = "Recall (TPR)",
  specificity = "Specificity", f1 = "F1", mcc = "MCC",
  wolbers_c = "C-index (Wolbers)")

selection_metrics <- function(performance, models) {
  models <- resolve_models(models, performance$Model)
  d <- performance[performance$Model %in% models$model, ]
  ok <- .usable_runs(d)
  if (any(!ok)) {
    n <- table(d$Model[!ok])
    message("Dropping failed fits, runs per model: ",
            paste(sprintf("%s = %d", names(n), n), collapse = ", "))
  }
  d <- d[ok, ]

  tp <- as.numeric(d$tp); fp <- as.numeric(d$fp)
  tn <- as.numeric(d$tn); fn <- as.numeric(d$fn)
  den <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))

  out <- data.frame(
    Model       = d$Model,
    run_number  = d$run_number,
    tp = tp, fp = fp, tn = tn, fn = fn,
    n_selected  = tp + fp,
    precision   = ifelse(tp + fp > 0, tp / (tp + fp), NA_real_),
    recall      = ifelse(tp + fn > 0, tp / (tp + fn), NA_real_),
    specificity = ifelse(tn + fp > 0, tn / (tn + fp), NA_real_),
    f1          = ifelse(2 * tp + fp + fn > 0, 2 * tp / (2 * tp + fp + fn),
                         NA_real_),
    mcc         = ifelse(den > 0, (tp * tn - fp * fn) / den, NA_real_),
    wolbers_c   = d$wolbers_c,
    stringsAsFactors = FALSE)
  out <- label_models(out, models)
  attr(out, "models") <- models
  out
}

summarise_selection <- function(sel,
                                metrics = c("tp", "fp", "fn", "n_selected",
                                            "precision", "recall", "specificity",
                                            "f1", "mcc", "wolbers_c"),
                                digits = 3) {
  fmt <- function(x) {
    x <- x[!is.na(x)]
    if (!length(x)) return("")
    sprintf("%s (%s)", format(signif(mean(x), digits), trim = TRUE),
            format(signif(stats::sd(x), 2), trim = TRUE))
  }
  out <- do.call(rbind, lapply(split(sel, sel$Model, drop = TRUE), function(d)
    data.frame(Model = d$Model[1L], Runs = nrow(d),
               t(vapply(metrics, function(m) fmt(d[[m]]), "")),
               check.names = FALSE, stringsAsFactors = FALSE)))
  names(out) <- c("Model", "Runs", .metric_labels[metrics])
  rownames(out) <- NULL
  out
}

mean_confusion <- function(sel) {
  agg <- stats::aggregate(cbind(tp, fp, fn, tn) ~ Model, data = sel, FUN = mean)
  long <- rbind(
    data.frame(Model = agg$Model, Truth = "Active",   Selected = "Selected",     n = agg$tp),
    data.frame(Model = agg$Model, Truth = "Inactive", Selected = "Selected",     n = agg$fp),
    data.frame(Model = agg$Model, Truth = "Active",   Selected = "Not selected", n = agg$fn),
    data.frame(Model = agg$Model, Truth = "Inactive", Selected = "Not selected", n = agg$tn))
  long$Truth    <- factor(long$Truth, levels = c("Active", "Inactive"))
  long$Selected <- factor(long$Selected, levels = c("Not selected", "Selected"))
  # Share of the truth class, so the colour reads the same for the 30-row and
  # the 970-row column.
  tot <- stats::ave(long$n, long$Model, long$Truth, FUN = sum)
  long$share <- long$n / tot
  long$correct <- (long$Truth == "Active") == (long$Selected == "Selected")
  long
}

plot_confusion <- function(sel) {
  cm <- mean_confusion(sel)
  cm$fill <- ifelse(cm$correct, cm$share, -cm$share)
  ggplot2::ggplot(cm, ggplot2::aes(x = Truth, y = Selected, fill = fill)) +
    ggplot2::geom_tile(colour = "white", linewidth = 1) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f\n(%.1f%%)", n,
                                                    100 * share)),
                       size = 3) +
    ggplot2::facet_wrap(~ Model, nrow = 1) +
    ggplot2::scale_fill_gradient2(low = "#D6604D", mid = "white",
                                  high = "#4393C3", midpoint = 0,
                                  limits = c(-1, 1), guide = "none") +
    ggplot2::labs(x = "Truth", y = NULL,
                  title = "Mean confusion matrix over runs",
                  caption = paste("Cells: mean count over usable runs, with the",
                                  "percentage of that truth column. Blue = correct,",
                                  "red = error; shade = share of the column.")) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(face = "bold"),
                   plot.caption = ggplot2::element_text(colour = "grey40",
                                                        hjust = 0))
}

plot_selection <- function(sel,
                           metrics = c("recall", "precision", "f1", "mcc",
                                       "fp", "wolbers_c"),
                           models = attr(sel, "models")) {
  long <- do.call(rbind, lapply(metrics, function(m)
    data.frame(Model = sel$Model, metric = .metric_labels[[m]],
               value = sel[[m]], stringsAsFactors = FALSE)))
  long$metric <- factor(long$metric, levels = .metric_labels[metrics])
  long <- long[!is.na(long$value), ]
  ggplot2::ggplot(long, ggplot2::aes(x = Model, y = value, fill = Model)) +
    ggplot2::geom_boxplot(outlier.size = 0.8, width = 0.65, alpha = 0.85) +
    ggplot2::facet_wrap(~ metric, scales = "free_y", ncol = 3) +
    scale_fill_models(models) +
    ggplot2::labs(x = NULL, y = NULL, fill = NULL,
                  title = "Selection and prediction metrics across runs") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   legend.position = "none",
                   strip.text = ggplot2::element_text(face = "bold"))
}

# ---------------------------------------------------------------------------
# By correlation block
# ---------------------------------------------------------------------------

selection_by_block <- function(by_block, performance, models) {
  models <- resolve_models(models, performance$Model)
  ok <- performance[.usable_runs(performance), c("Model", "run_number")]
  d <- merge(by_block, ok, by = c("Model", "run_number"))
  d <- d[d$Model %in% models$model, ]
  d$recall  <- ifelse(d$tp + d$fn > 0, d$tp / (d$tp + d$fn), NA_real_)
  d$fp_rate <- ifelse(d$fp + d$tn > 0, d$fp / (d$fp + d$tn), NA_real_)
  d <- label_models(d, models)
  d$Block <- .block_factor(d$Block, d$BlockRho)
  attr(d, "models") <- models
  d
}

summarise_selection_by_block <- function(sbb) {
  out <- sbb |>
    dplyr::group_by(Block, Model) |>
    dplyr::summarise(
      Runs = dplyr::n(),
      `Mean TP` = mean(tp), `Mean FN` = mean(fn), `Mean FP` = mean(fp),
      `Recall` = mean(recall, na.rm = TRUE),
      `FP rate (%)` = 100 * sum(fp) / sum(fp + tn),
      .groups = "drop") |>
    as.data.frame()
  out
}

plot_selection_by_block <- function(sbb, models = attr(sbb, "models")) {
  long <- rbind(
    data.frame(Model = sbb$Model, Block = sbb$Block,
               metric = "Recall (TPR) among actives", value = sbb$recall),
    data.frame(Model = sbb$Model, Block = sbb$Block,
               metric = "False positive rate among nulls (%)",
               value = 100 * sbb$fp_rate))
  long <- long[!is.na(long$value), ]
  long$metric <- factor(long$metric, levels = unique(long$metric))
  s <- long |>
    dplyr::group_by(Model, Block, metric) |>
    dplyr::summarise(mean = mean(value),
                     se = stats::sd(value) / sqrt(dplyr::n()),
                     .groups = "drop")
  dodge <- ggplot2::position_dodge(width = 0.6)
  ggplot2::ggplot(s, ggplot2::aes(x = Block, y = mean, colour = Model)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - 1.96 * se,
                                        ymax = mean + 1.96 * se),
                           position = dodge, width = 0) +
    ggplot2::geom_point(position = dodge, size = 2.2) +
    ggplot2::facet_wrap(~ metric, scales = "free_y") +
    scale_colour_models(models) +
    ggplot2::labs(x = "Correlation block", y = NULL, colour = NULL,
                  title = "Selection by correlation block",
                  caption = "Points: mean over runs. Bars: +/- 1.96 Monte Carlo SE.") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom",
                   panel.grid.major.x = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(face = "bold"),
                   plot.caption = ggplot2::element_text(colour = "grey40",
                                                        hjust = 0))
}
