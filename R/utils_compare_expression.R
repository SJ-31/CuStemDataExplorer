merge_expr_tbs <- function(x, y) {
  list(
    expr = dplyr::left_join(
      x$expr,
      y$expr,
      by = join_by(gene_id)
    ),
    meta = dplyr::bind_rows(x$meta, y$meta)
  )
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
    tb <- ifelse(
      stringr::str_ends(spec$counts, "csv"),
      readr::read_csv(spec$counts),
      readr::read_tsv(spec$counts)
    )
    gene_col <- spec$gene_name_format %||% "gene_id"
    cohort <- spec$gene_name_format %||% "unassigned"
    pointblank::col_exists(tb, id)
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
    result <- lapply(tt_list, \(t) read_helper(t, ttype = tumor_type)) |>
      purrr::reduce(merge_expr_tbs)
  }) |>
    purrr::reduce(merge_expr_tbs)
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
    purrr::reduce(merge_expr_tbs)
}
