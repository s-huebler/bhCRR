# sim/timing/plot_timing.r
#
# Visualise ClinTALL fit-time scaling results.
# Reads ONE collected CSV; tag is inferred from the filename if not set.
#
# Usage (from repo root):
#   Rscript sim/timing/plot_timing.r
#   TIMING_CSV=sim/timing/results/clintall_timing_<tag>.csv Rscript sim/timing/plot_timing.r
#
# Outputs (next to the CSV):
#   clintall_timing_<tag>.pdf
#   clintall_timing_<tag>.png

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# ---------------------------------------------------------------------------
# 1. Locate CSV
# ---------------------------------------------------------------------------
csv_env <- Sys.getenv("TIMING_CSV", unset = "")
if (nzchar(csv_env)) {
  csv_path <- csv_env
} else {
  hits <- list.files("sim/timing/results",
                     pattern = "^clintall_timing_.*\\.csv$",
                     full.names = TRUE)
  hits <- hits[!grepl("calib", hits)]   # exclude any *_calib CSVs
  if (length(hits) == 0L)
    stop("No clintall_timing_*.csv found in sim/timing/results/. ",
         "Set TIMING_CSV= or place the file there.")
  if (length(hits) > 1L) {
    message("Multiple CSVs found; using the most recently modified:")
    hits <- hits[order(file.mtime(hits), decreasing = TRUE)]
    message("  ", hits[1])
  }
  csv_path <- hits[1]
}

cat("[plot] Reading:", csv_path, "\n")
d <- read.csv(csv_path, stringsAsFactors = FALSE)
stopifnot(nrow(d) > 0L)

out_stem <- sub("\\.csv$", "", csv_path)

# ---------------------------------------------------------------------------
# 2. Derive subtitle metadata from the CSV (never hard-coded)
# ---------------------------------------------------------------------------
s0_val      <- unique(d$s0)[1]
s1_val      <- unique(d$s1)[1]
init_meth   <- unique(d$init_method)[1]
cpu_val     <- unique(d$cpu_model[!is.na(d$cpu_model)])[1]
host_val    <- unique(d$host[!is.na(d$host)])[1]
run_dates   <- as.Date(substr(d$timestamp, 1, 10))
run_date    <- format(min(run_dates, na.rm = TRUE), "%Y-%m-%d")
n_val       <- unique(d$n)[1]

subtitle_str <- sprintf(
  "n = %d  |  s0 = %g, s1 = %g  |  init: %s  |  %s on %s  |  run %s",
  n_val, s0_val, s1_val, init_meth,
  ifelse(is.na(cpu_val), "unknown CPU", cpu_val),
  ifelse(is.na(host_val), "unknown host", host_val),
  run_date
)

# ---------------------------------------------------------------------------
# 3. Summaries
# ---------------------------------------------------------------------------
# Mark censored rows (TIMEOUT)
d$censored <- !is.na(d$err) & d$err == "TIMEOUT"

# Numeric time columns — coerce NAs that came from TIMEOUT
d$t_total <- as.numeric(d$t_total)
d$t_init  <- as.numeric(d$t_init)
d$t_em    <- as.numeric(d$t_em)

# Medians across reps (uncensored only for the median line)
med_df <- do.call(rbind, lapply(split(d, d$p), function(sub) {
  unc <- sub[!sub$censored & !is.na(sub$t_total), ]
  data.frame(
    p          = sub$p[1],
    med_total  = if (nrow(unc) > 0) median(unc$t_total) else NA_real_,
    med_init   = if (nrow(unc) > 0) median(unc$t_init, na.rm = TRUE) else NA_real_,
    med_em     = if (nrow(unc) > 0) median(unc$t_em,   na.rm = TRUE) else NA_real_,
    med_nnz    = if (nrow(unc) > 0) median(unc$nnz,    na.rm = TRUE) else NA_real_,
    med_iter   = if (nrow(unc) > 0) median(unc$iterations, na.rm = TRUE) else NA_real_
  )
}))
med_df <- med_df[order(med_df$p), ]

any_zero_nnz <- any(!is.na(d$nnz) & d$nnz == 0L)

# ---------------------------------------------------------------------------
# 4. Shared theme
# ---------------------------------------------------------------------------
theme_timing <- function() {
  theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold"),
      legend.position  = "bottom"
    )
}

# ---------------------------------------------------------------------------
# 5. Panel A: t_total (seconds) vs p
# ---------------------------------------------------------------------------
# Long-form for t_init / t_em overlay
d_long <- rbind(
  data.frame(p = d$p, t = d$t_total, rep = d$rep,
             series = "t_total", censored = d$censored),
  data.frame(p = d$p, t = d$t_init,  rep = d$rep,
             series = "t_init",  censored = d$censored),
  data.frame(p = d$p, t = d$t_em,    rep = d$rep,
             series = "t_em",    censored = d$censored)
)
d_long <- d_long[!is.na(d_long$t), ]

med_long <- rbind(
  data.frame(p = med_df$p, med = med_df$med_total, series = "t_total"),
  data.frame(p = med_df$p, med = med_df$med_init,  series = "t_init"),
  data.frame(p = med_df$p, med = med_df$med_em,    series = "t_em")
)
med_long <- med_long[!is.na(med_long$med), ]

series_cols   <- c(t_total = "#2166ac", t_init = "#d6604d", t_em = "#4dac26")
series_labels <- c(t_total = "t_total (EM + init)",
                   t_init  = "t_init (LASSO BIC)",
                   t_em    = "t_em (EM only)")

p_A <- ggplot() +
  # uncensored points
  geom_point(
    data = d_long[!d_long$censored, ],
    aes(x = p, y = t, colour = series),
    alpha = 0.55, size = 1.8, shape = 16
  ) +
  # censored as up-pointing triangles at their time limit
  geom_point(
    data = d_long[d_long$censored & d_long$series == "t_total", ],
    aes(x = p, y = t),
    shape = 24, size = 2.5, fill = "grey40", colour = "grey20"
  ) +
  # median lines
  geom_line(
    data = med_long,
    aes(x = p, y = med, colour = series),
    linewidth = 0.9
  ) +
  scale_colour_manual(values = series_cols, labels = series_labels,
                      name = NULL) +
  scale_x_continuous(
    labels = function(x) ifelse(x >= 1000, paste0(x / 1000, "k"), x)
  ) +
  labs(
    title    = "fit_ssl_psdh wall time vs number of gene predictors",
    subtitle = subtitle_str,
    x        = "p (gene predictors)",
    y        = "elapsed time (s)"
  ) +
  theme_timing()

# ---------------------------------------------------------------------------
# 6. Panel B: log–log with fitted slope
# ---------------------------------------------------------------------------
ll_df <- med_df[!is.na(med_df$med_total) & med_df$med_total > 0, ]
lm_fit <- lm(log(med_total) ~ log(p), data = ll_df)
slope  <- round(coef(lm_fit)[["log(p)"]], 3)
r2     <- round(summary(lm_fit)$r.squared, 3)

ll_df$fitted <- exp(predict(lm_fit))

p_B <- ggplot(ll_df, aes(x = p, y = med_total)) +
  geom_line(aes(y = fitted), colour = "grey50", linetype = "dashed",
            linewidth = 0.8) +
  geom_point(colour = series_cols["t_total"], size = 2.2) +
  geom_line(colour  = series_cols["t_total"], linewidth = 0.8) +
  annotate("text", x = min(ll_df$p) * 1.4, y = max(ll_df$med_total) * 0.6,
           label = sprintf("slope = %.3f\nR² = %.3f", slope, r2),
           hjust = 0, size = 3.5, colour = "grey30") +
  scale_x_log10(labels = function(x) ifelse(x >= 1000, paste0(x/1000,"k"), x)) +
  scale_y_log10() +
  labs(
    title = sprintf("Log-log: empirical scaling exponent = %.3f", slope),
    x     = "p (log scale)",
    y     = "median t_total (s, log scale)"
  ) +
  theme_timing() +
  theme(legend.position = "none")

# ---------------------------------------------------------------------------
# 7. Panel C: EM iterations and nnz vs p
# ---------------------------------------------------------------------------
# nnz per rep (points) + median line; flag nnz = 0
d_c <- d[!is.na(d$nnz), ]
nnz_zero_p <- unique(d_c$p[d_c$nnz == 0L])

nnz_med <- med_df[!is.na(med_df$med_nnz), ]
iter_med <- med_df[!is.na(med_df$med_iter), ]

# Dual-series within one panel (nnz left axis, iterations right)
# Use a single panel with two geom layers scaled by a common factor
iter_scale <- max(d_c$nnz, na.rm = TRUE) /
              max(d_c$iterations, na.rm = TRUE)

p_C <- ggplot() +
  # nnz points + median line
  geom_point(data = d_c,
             aes(x = p, y = nnz),
             colour = "#762a83", alpha = 0.55, size = 1.8) +
  geom_line(data = nnz_med,
            aes(x = p, y = med_nnz),
            colour = "#762a83", linewidth = 0.9) +
  # iterations points + median line (scaled to nnz axis)
  geom_point(data = d_c,
             aes(x = p, y = iterations * iter_scale),
             colour = "#1b7837", alpha = 0.55, size = 1.8, shape = 17) +
  geom_line(data = iter_med,
            aes(x = p, y = med_iter * iter_scale),
            colour = "#1b7837", linewidth = 0.9, linetype = "dashed") +
  # flag zero-nnz p values if any
  { if (any_zero_nnz)
      geom_vline(xintercept = nnz_zero_p, colour = "red",
                 linetype = "dotted", linewidth = 0.8)
    else list() } +
  scale_x_continuous(
    labels = function(x) ifelse(x >= 1000, paste0(x / 1000, "k"), x)
  ) +
  scale_y_continuous(
    name = "nnz (non-zero coefficients)",
    sec.axis = sec_axis(~ . / iter_scale,
                        name = "EM iterations",
                        breaks = seq(1, 10))
  ) +
  labs(
    title = "Model complexity vs p",
    x     = "p (gene predictors)",
    caption = if (any_zero_nnz)
      paste("Red dotted line: p with nnz = 0:", paste(nnz_zero_p, collapse = ", "))
    else
      "Purple circles: nnz (left axis)  |  Green triangles: EM iterations (right axis)"
  ) +
  theme_timing() +
  theme(
    axis.title.y.right = element_text(colour = "#1b7837"),
    axis.text.y.right  = element_text(colour = "#1b7837"),
    axis.title.y.left  = element_text(colour = "#762a83"),
    axis.text.y.left   = element_text(colour = "#762a83")
  )

# ---------------------------------------------------------------------------
# 8. Assemble and save
# ---------------------------------------------------------------------------
fig <- (p_A / (p_B | p_C)) +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(face = "bold"))
  )

pdf_path <- paste0(out_stem, ".pdf")
png_path <- paste0(out_stem, ".png")

ggsave(pdf_path, fig, width = 12, height = 10)
ggsave(png_path, fig, width = 12, height = 10, dpi = 150)

cat("[plot] Saved:", pdf_path, "\n")
cat("[plot] Saved:", png_path, "\n")
cat(sprintf("[plot] Empirical scaling exponent: %.3f  (R² = %.3f)\n", slope, r2))
if (any_zero_nnz) {
  cat("[plot] WARNING: nnz = 0 at p =", paste(nnz_zero_p, collapse = ", "), "\n")
} else {
  cat("[plot] nnz > 0 at all p values.\n")
}
