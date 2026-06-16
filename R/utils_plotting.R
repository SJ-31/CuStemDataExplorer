# Miscellaneous plotting functions used by multiple modules

#' Plot PCA
#'
#' @param pr_obj Object of type `prcomp`
#' @param expr_tbs Expression tibbles, list with keys `meta`, `expr`
#' @param color_by Column in expr_tbs$meta used to color points
plot_pca <- function(pr_obj, expr_tbs, color_by) {
  box::use(ggplot2[ggplot, aes])
  joined <- tibble::as_tibble(pr_obj$x, rownames = "sample") |>
    dplyr::inner_join(
      dplyr::select(expr_tbs$meta, dplyr::all_of(c("sample", color_by))),
      by = dplyr::join_by(sample)
    )
  summary_df <- as.data.frame(summary(pr_obj)$importance)
  pc1_var <- summary_df$PC1[2]
  pc2_var <- summary_df$PC2[2]

  ggplot(joined, aes(x = PC1, y = PC2, color = !!as.symbol(color_by))) +
    ggplot2::geom_point() +
    ggplot2::xlab(paste0("PC1 (", pc1_var, ")")) +
    ggplot2::ylab(paste0("PC2 (", pc2_var, ")"))
}
