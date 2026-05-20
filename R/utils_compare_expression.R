#' @importFrom rlang %||%

merge_expr_tbs <- function(tb_list) {
  tb_list <- purrr::discard(tb_list, is.null)
  if (length(tb_list) > 1) {
    purrr::reduce(tb_list, \(x, y) {
      list(
        expr = dplyr::left_join(
          x$expr,
          y$expr,
          by = dplyr::join_by(gene_id)
        ),
        meta = dplyr::bind_rows(x$meta, y$meta)
      )
    })
  } else if (length(tb_list) == 1) {
    tb_list[[1]]
  } else {
    NULL
  }
}

#' Read an h5ad file and pseudobulk it
#'
#' @return List of two tibbles, `expr` and `meta`
#'
read_anndata_pb <- function(
  file,
  cache = TRUE,
  convert_names_to = "symbol",
  sample_col = "sample"
) {
  box::use(SingleCellExperiment[colData])
  adata <- anndataR::read_h5ad(file, as = "SingleCellExperiment")
  bulked <- scrapper::aggregateAcrossCells.se(
    adata,
    list(sample = colData(adata)[[sample_col]]),
    assay.type = "X"
  )
  sample_names <- colData(bulked)[[sample_col]]
  expr <- SummarizedExperiment::assay(bulked, "sums") |>
    as.data.frame() |>
    `colnames<-`(sample_names) |>
    tibble::rownames_to_column(var = "gene_id") |>
    tidyr::as_tibble()
  meta <- colData(bulked) |>
    tidyr::as_tibble() |>
    dplyr::select(dplyr::any_of(
      c(sample_col, "cohort", "tumor_type", "treatment", "patient")
    ))
  if ("cohort" %notin% colnames(meta)) {
    meta$cohort <- NA_character_
  }
  if ("tumor_type" %notin% colnames(meta)) {
    meta$tumor_type <- NA_character_
  }
  list(expr = expr, meta = meta)
}

#' Helper function to read expression data from yaml configuration
#'
#' @return
#'
read_expression_spec <- function(
  file,
  convert_names_to = "symbol",
  allowed_exts = c("csv", "tsv", "h5ad"),
  blacklist = NULL
) {
  all_spec <- yaml::read_yaml(file)

  read_helper <- function(spec, ttype) {
    cohort_val <- spec$cohort %||% "unassigned"
    checkmate::assert_list(spec)
    checkmate::assert_names(
      names(spec),
      subset.of = c(
        "cohort",
        "patient",
        "treatment",
        "counts",
        "gene_name_format",
        "gene_col"
      )
    )
    checkmate::assert_file_exists(spec$counts)
    meta <- NULL
    ext <- tools::file_ext(spec$counts)
    if (ext %notin% allowed_exts) {
      return(NULL)
    } else if (ext == "csv") {
      tb <- suppressMessages(readr::read_csv(spec$counts))
    } else if (ext == "tsv") {
      tb <- suppressMessages(readr::read_tsv(spec$counts))
    } else if (ext == "h5ad") {
      from_anndata <- read_anndata_pb(spec$counts)
      tb <- from_anndata$expr
      meta <- from_anndata$meta
      meta <- meta |>
        dplyr::mutate(
          tumor_type = dplyr::replace_values(tumor_type, NA ~ ttype),
          cohort = dplyr::replace_values(cohort, NA ~ cohort_val),
        )
    } else {
      stop("Extension not supported yet")
    }
    gene_col <- spec$gene_col %||% "gene_id"
    gene_name_format <- spec$gene_name_format %||% "ensembl"
    pointblank::col_exists(tb, gene_col)
    tb <- dplyr::rename(tb, gene_id = gene_col)
    tb$gene_id <- recode_genes(
      tb$gene_id,
      to = convert_names_to,
      from = gene_name_format
    )
    tb <- dplyr::distinct(tb, gene_id, .keep_all = TRUE) |>
      dplyr::filter(gene_id %notin% blacklist)
    samples <- colnames(tb) |> purrr::discard(\(x) x == gene_col)
    if (is.null(meta)) {
      meta <- tibble::tibble(
        sample = samples,
        cohort = cohort_val,
        tumor_type = ttype
      )
    }
    list(expr = tb, meta = meta)
  }

  lapply(names(all_spec), \(tumor_type) {
    tt_list <- all_spec[[tumor_type]]
    lapply(tt_list, \(t) read_helper(t, ttype = tumor_type)) |>
      merge_expr_tbs()
  }) |>
    merge_expr_tbs()
}

read_all_expr <- function(
  cfg,
  allowed_exts = c("csv", "tsv", "h5ad"),
  blacklist = c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")
) {
  dir <- cfg$spec_directory
  checkmate::assert_directory_exists(dir)
  convert_names_to <- cfg$name_format %||% "symbol"

  combined <- lapply(
    list.files(dir, pattern = ".yml|yaml$", full.names = TRUE),
    \(f) {
      read_expression_spec(
        f,
        convert_names_to,
        allowed_exts = allowed_exts,
        blacklist = blacklist
      )
    }
  ) |>
    merge_expr_tbs()
}

#' Helper function to remove samples from expr_tbs
#'
#' @description
#' @param expr_tbs List with two elements: expr containing the expression
#' tibble and meta containing the sample metadata
#' @param fn Filter function applied to metadata
#' e.g. \(x) filter(x, x$tumor_type == "HCC")
filter_expr_tbs <- function(expr_tbs, fn) {
  meta <- fn(expr_tbs$meta)
  samples_keep <- meta$sample
  kept <- dplyr::select(expr_tbs$expr, "gene_id", dplyr::all_of(samples_keep))
  list(
    expr = kept,
    meta = meta
  )
}


#' Compute TMM normalization factors from expression with TMM
#' convert to (log) CPM
#'
tmm_normalize <- function(obj, log = TRUE) {
  dge <- local({
    tmp <- tibble::column_to_rownames(obj$expr, var = "gene_id")
    edgeR::DGEList(counts = tmp)
  })
  dge <- edgeR::normLibSizes(dge)
  counts <- edgeR::cpm(dge, log = log) |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "gene_id") |>
    tibble::as_tibble()
  list(expr = counts, meta = obj$meta)
}


do_heatmap <- function(
  expr_tbs,
  genes,
  cfg,
  tumor_types = NULL,
  cohorts = NULL
) {
  box::use(
    ggplot2[aes, theme, element_blank, geom_tile, theme_void, ggplot],
    paletteer[scale_color_paletteer_c, scale_fill_paletteer_d],
    dplyr[filter]
  )
  checkmate::assert_list(cfg)

  long <- expr_tbs$expr |>
    filter(gene_id %in% genes) |>
    tidyr::pivot_longer(-gene_id, names_to = "sample") |>
    dplyr::left_join(expr_tbs$meta, by = dplyr::join_by(sample))
  meta <- expr_tbs$meta
  if (!is.null(tumor_types)) {
    long <- filter(long, tumor_type %in% tumor_types)
    meta <- filter(meta, tumor_type %in% tumor_types)
  }
  if (!is.null(cohorts)) {
    long <- filter(long, cohort %in% cohorts)
    meta <- filter(meta, cohort %in% cohorts)
  }

  expr_palette <- cfg$expression %||% "ggthemes::Red-Gold"

  n_cohorts <- length(unique(long$cohort))
  n_tumor_types <- length(unique(long$tumor_type))
  cohort_palette <- random_palette_d(n_cohorts)
  ttype_palette <- random_palette_d(n_tumor_types)

  top_theming <- theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  )
  bot_theming <- theme(
    axis.text.x = ggplot2::element_text(angle = 90),
    axis.title.x = ggplot2::element_text(size = 15)
  )

  expr_plot <- ggplot(
    long,
    aes(x = sample, y = gene_id, fill = value)
  ) +
    geom_tile(width = 0.90) +
    theme(panel.grid = element_blank()) +
    paletteer::scale_fill_paletteer_c(expr_palette) +
    ggplot2::guides(fill = ggplot2::guide_legend("Normalized expression")) +
    ggplot2::ylab("Gene")

  if (n_cohorts <= 1 && n_tumor_types <= 1) {
    return(expr_plot + bot_theming + ggplot2::xlab("Sample"))
  }
  if (n_cohorts > 1) {
    cohort_labels <- ggplot(meta, aes(x = sample, fill = cohort, y = "1")) +
      geom_tile() +
      theme_void() +
      scale_fill_paletteer_d(cohort_palette) +
      ggplot2::guides(fill = ggplot2::guide_legend("Cohort"))
  }
  if (n_tumor_types > 1) {
    ttype_labels <- ggplot(meta, aes(x = sample, fill = tumor_type, y = "1")) +
      geom_tile() +
      theme_void() +
      scale_fill_paletteer_d(ttype_palette) +
      ggplot2::guides(fill = ggplot2::guide_legend("Tumor type"))
  }

  if (n_tumor_types > 1 && n_cohorts <= 1) {
    patchwork::wrap_plots(
      expr_plot + top_theming,
      ttype_labels + bot_theming + ggplot2::xlab("Sample"),
      nrow = 2,
      heights = c(0.95, 0.05),
      guides = "collect"
    )
  } else if (n_cohorts > 1 && n_tumor_types <= 1) {
    patchwork::wrap_plots(
      expr_plot + top_theming,
      cohort_labels + bot_theming + ggplot2::xlab("Sample"),
      nrow = 2,
      heights = c(0.95, 0.05),
      guides = "collect"
    )
  } else {
    patchwork::wrap_plots(
      expr_plot + top_theming,
      cohort_labels + top_theming,
      ttype_labels + bot_theming + ggplot2::xlab("Sample"),
      nrow = 3,
      heights = c(0.8, 0.05, 0.05),
      guides = "collect"
    )
  }
}
