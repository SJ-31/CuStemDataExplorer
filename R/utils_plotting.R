# Miscellaneous plotting functions used by multiple modules

#' Plot PCA
#'
#' @param pr_obj Object of type `prcomp`
#' @param expr_tbs Expression tibbles, list with keys `meta`, `expr`
#' @param color_by Column in expr_tbs$meta used to color points
plot_pca <- function(
  pr_obj,
  expr_tbs,
  color_by,
  key_group,
  tumor_types = NULL,
  cohorts = NULL
) {
  if (nchar(color_by) == 0) {
    return()
  }
  box::use(ggplot2[ggplot, aes])
  joined <- tibble::as_tibble(pr_obj$x, rownames = "sample") |>
    dplyr::inner_join(
      expr_tbs$meta,
      by = dplyr::join_by(sample)
    )
  n_labels <- length(unique(joined[[color_by]]))
  if (!is.null(tumor_types)) {
    joined <- joined |> dplyr::filter(tumor_type %in% tumor_types)
  }
  if (!is.null(cohorts)) {
    joined <- joined |> dplyr::filter(cohort %in% cohorts)
  }

  summary_df <- as.data.frame(summary(pr_obj)$importance)
  pc1_var <- summary_df$PC1[2]
  pc2_var <- summary_df$PC2[2]

  plot <- ggplot(joined, aes(x = PC1, y = PC2, color = !!as.symbol(color_by))) +
    ggplot2::geom_point() +
    ggplot2::xlab(paste0("PC1 (", pc1_var, ")")) +
    ggplot2::ylab(paste0("PC2 (", pc2_var, ")"))
  plot |>
    add_palette_from_cache(
      key = color_by,
      key_group = key_group,
      min_length = n_labels,
      fill = FALSE
    )
}
