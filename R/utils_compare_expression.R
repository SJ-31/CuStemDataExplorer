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
    samples <- colnames(tb) |> purrr::discard(\(x) x == gene_col)
    meta <- tibble::tibble(
      samples = samples,
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

read_all_expr <- function() {
  dir <- get_golem_config("expression_viewer")$spec_directory
  checkmate::assert_directory_exists(dir)
  convert_names_to <- get_golem_config("expression_viewer")$name_format %||%
    "symbol"

  lapply(
    list.files(dir, pattern = ".yml|yaml$"),
    \(f) read_expression_spec(f, convert_names_to)
  ) |>
    merge_expr_tbs()
}
