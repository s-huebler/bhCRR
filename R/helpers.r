


plot_score_heatmap <- function(df) {

  # Ensure the dataframe is valid
  if (!all(c("s0", "s1", "score_mean") %in% colnames(df))) {
    stop("Dataframe must contain 's0', 's1', and 'score_mean' columns.")
  }

  # Create the heatmap using ggplot2
  p <- ggplot(data = df, aes(x = as.factor(s0), y = as.factor(s1), fill = score_mean)) +
    geom_tile(color = "white", size = 0.5) + # Adds white borders between tiles
    scale_fill_viridis_c(option = "plasma", name = "Score Mean") + # Applies a colorblind-friendly continuous color scale
    labs(
      title = "Heatmap of Score Mean by s0 and s1",
      x = "s0",
      y = "s1"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(), # Removes background grid lines
      axis.text.x = element_text(angle = 45, hjust = 1) # Tilts x-axis labels if there are many
    )

  return(p)
}


#' Permute the columns of each row independently
#'
#' For each row of \code{df}, independently shuffles the column values.
#' Used to break between-column correlations while preserving the marginal
#' distribution of each row.
#'
#' @param df A numeric data frame (or object coercible to a numeric matrix).
#' @param seed Optional integer seed passed to \code{\link{set.seed}} for
#'   reproducibility.  \code{NULL} (default) leaves the RNG state unchanged.
#'
#' @return A data frame with the same dimensions and column/row names as
#'   \code{df}, with each row's values independently permuted across columns.
#'
#' @examples
#' set.seed(1)
#' df <- as.data.frame(matrix(1:12, nrow = 3))
#' permute_rows(df, seed = 42)
#'
#' @export
permute_rows <- function(df, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  m <- as.matrix(df)
  permuted_mat <- t(apply(m, 1, function(row) row[sample.int(length(row))]))
  colnames(permuted_mat) <- colnames(m)
  rownames(permuted_mat) <- rownames(m)
  as.data.frame(permuted_mat)
}
