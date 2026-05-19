#' @importFrom rlang %||%

merge_expr_tbs <- function(tb_list) {
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
  } else {
    tb_list[[1]]
  }
}


#' Helper function to read expression data from yaml configuration
#'
#' @description
#'
read_expression_spec <- function(file, convert_names_to = "symbol") {
  all_spec <- yaml::read_yaml(file)

  read_helper <- function(spec, ttype) {
    checkmate::assert_list(spec)
    checkmate::assert_names(
      names(spec),
      subset.of = c(
        "cohort",
        "counts",
        "gene_name_format",
        "gene_col"
      )
    )
    checkmate::assert_file_exists(spec$counts)
    if (stringr::str_ends(spec$counts, "csv")) {
      tb <- readr::read_csv(spec$counts)
    } else {
      tb <- readr::read_tsv(spec$counts)
    }
    gene_col <- spec$gene_col %||% "gene_id"
    gene_name_format <- spec$gene_name_format %||% "ensembl"
    cohort <- spec$cohort %||% "unassigned"
    pointblank::col_exists(tb, gene_col)
    tb <- dplyr::rename(tb, gene_id = gene_col)
    tb$gene_id <- recode_genes(
      tb$gene_id,
      to = convert_names_to,
      from = gene_name_format
    )
    samples <- colnames(tb) |> purrr::discard(\(x) x == gene_col)
    meta <- tibble::tibble(
      sample = samples,
      cohort = cohort,
      tumor_type = ttype
    )

    list(expr = tb, meta = meta)
  }

  lapply(names(all_spec), \(tumor_type) {
    tt_list <- all_spec[[tumor_type]]
    lapply(tt_list, \(t) read_helper(t, ttype = tumor_type)) |>
      merge_expr_tbs()
  }) |>
    merge_expr_tbs()
}

read_all_expr <- function(cfg) {
  dir <- cfg$spec_directory
  checkmate::assert_directory_exists(dir)
  convert_names_to <- cfg$name_format %||% "symbol"

  combined <- lapply(
    list.files(dir, pattern = ".yml|yaml$", full.names = TRUE),
    \(f) read_expression_spec(f, convert_names_to)
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


tmm_normalize <- function(obj) {}

# TODO: would also like a single comparison

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
