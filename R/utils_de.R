# TODO: these probably work for edgeR too, would be useful to keep

deseq_ovr_contrasts <- function(factor, dds) {
  box::use(stringr[str_starts])
  stopifnot("(Intercept)" %notin% DESeq2::resultsNames(dds))
  all_levels <- DESeq2::resultsNames(dds)
  levels <- all_levels[str_starts(all_levels, factor)]

  removed <- stringr::str_remove(levels, factor)

  n_levels <- length(levels)
  lapply(levels, \(l) {
    contrast <- rep(0, length(all_levels))
    contrast[str_starts(all_levels, factor)] <- 1 / (n_levels - 1)
    contrast[which(all_levels == l)] <- -1
    DESeq2::results(dds, contrast = contrast)
  }) |>
    `names<-`(glue::glue("{removed} vs. Rest"))
}

deseq_pairwise_contrasts <- function(factor, dds) {
  levels <- DESeq2::resultsNames(dds)
  levels <- levels[stringr::str_starts(levels, factor)] |>
    stringr::str_remove(factor)
  combos <- combn(unique(levels), 2)
  names <- apply(combos, 2, \(x) glue::glue("{x[1]} vs. {x[2]}")) |> unlist()
  apply(combos, 2, \(col) {
    DESeq2::results(dds, c(factor, col[1], col[2]))
  }) |>
    `names<-`(names)
}

deseq_wrapper <- function(obj, kws) {
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = tibble::column_to_rownames(obj$expr, var = "gene_id"),
    colData = obj$meta,
    # TODO: include batch as well
    # TODO: read up about the different ways to parameterize, including
    # housekeeping genes
    design = ~ 0 + tumor_type
  )
  dds <- DESeq2::DESeq(dds)
  ovr <- deseq_ovr_contrasts("tumor_type", dds)
  pairwise <- deseq_pairwise_contrasts("tumor_type", dds)
  list(dds = dds, ovr = ovr, pairwise = pairwise)
}

do_de <- function(obj, how = "deseq2", kws) {
  if (how == "deseq2") {
    deseq_wrapper(obj, kws)
  } else {
    stop("Not implemented yet")
  }
}


volcano_plot <- function(obj, key_group, p_threshold = 0.05, key = "volcano") {
  if (!is.data.frame(obj)) {
    tb <- obj |>
      as.data.frame() |>
      tibble::rownames_to_column(var = "gene_id")
  } else {
    tb <- obj
  }
  if ("is_sig" %notin% colnames(tb)) {
    tb$is_sig <- tb$padj < p_threshold
  }
  plot <- ggplot2::ggplot(
    tb,
    ggplot2::aes(x = log2FoldChange, y = -log(padj), color = is_sig)
  ) +
    ggplot2::geom_point() +
    ggplot2::guides(color = ggplot2::guide_legend(title = "Significant"))
  add_palette_from_cache(plot, key = key, key_group = key_group, min_length = 2)
}

de_scatter_plot <- function(
  expr_tbs,
  chosen_genes,
  contrast,
  key_group,
  key = "de_scatter",
  factor = "tumor_type"
) {
  splits <- stringr::str_split_1(contrast, "vs.") |>
    purrr::map_chr(stringr::str_trim)
  if (length(chosen_genes) > 20) {
    return("Can plot at most 20 genes")
  }
  first <- splits[1]
  if (splits[2] == "Rest") {
    filter_by <- unique(expr_tbs$meta[[factor]])
  } else {
    filter_by <- c(first, splits[2])
  }
  long <- expr_tbs$expr |>
    dplyr::filter(gene_id %in% chosen_genes) |>
    tidyr::pivot_longer(-gene_id, names_to = "sample") |>
    dplyr::left_join(expr_tbs$meta, by = dplyr::join_by(sample)) |>
    dplyr::filter(!!as.symbol(factor) %in% filter_by)

  # TODO: wanna have the reference level be the first one
  # TODO: center the jitter points on each box and decrease their widths
  plot <- ggplot2::ggplot(
    long,
    ggplot2::aes(x = gene_id, y = value, color = !!as.symbol(factor))
  ) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_jitter() +
    ggplot2::ylab("Normalized Expression") +
    ggplot2::xlab("Gene")

  add_palette_from_cache(
    plot,
    key = key,
    key_group = key_group,
    min_length = length(unique(expr_tbs$meta[[factor]]))
  )
}
