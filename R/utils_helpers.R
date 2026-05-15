select_filter <- function(values, name) {
  tags$select(
    onchange = sprintf(
      "Reactable.setFilter('main_tab', '%s', event.target.value || undefined)",
      name
    ),
    tags$option(value = "", "All"),
    lapply(unique(values), tags$option)
  )
}

bool_col <- reactable::colDef(
  filterable = TRUE,
  filterInput = select_filter,
  align = "center",
  cell = \(v) {
    if (v == "F") "\u274c" else "\u2714\ufe0f"
  }
)

recode_genes <- function(
  vec,
  from = "ensembl",
  to = "symbol"
) {
  valid_formats <- c("ensembl", "symbol", "entrez")
  checkmate::assert_choice(from, choices = valid_formats)
  checkmate::assert_choice(to, choices = valid_formats)
  checkmate::assert_false(from == to)
  ref <- StemDataExplorer::ensembl_genes_hg38
  if (from == "ensembl") {
    names <- ref$ensembl_gene_id
  } else if (from == "symbol") {
    names <- hgnc_symbol
  } else {
    names <- ref$entrezgene_id
  }
  if (to == "ensembl") {
    vals <- ref$ensembl_gene_id
  } else if (to == "symbol") {
    vals <- hgnc_symbol
  } else {
    vals <- ref$entrezgene_id
  }
  lookup <- setNames(vals, names)
  mapped <- lookup[vec]
  ifelse(mapped, is.na(mapped), vec, mapped)
}

random_palette_d <- function(min_length = NULL) {
  if (!is.null(length)) {
    choices <- paletteer::palettes_d_names |> filter(length >= min_length)
  } else {
    choices <- paletteer::palettes_d_names
  }
  choices |>
    dplyr::slice_sample(n = 1) |>
    dplyr::select(package, palette) |>
    paste0(collapse = "::")
}
