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
    results(dds, contrast = contrast)
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
    countData = column_to_rownames(obj$expr, var = "gene_id"),
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
    deseq2_wrapper(obj, kws)
  } else {
    stop("Not implemented yet")
  }
}
