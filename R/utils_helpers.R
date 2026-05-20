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
  ref <- ensembl_genes_hg38
  if (from == "ensembl") {
    names <- ref$ensembl_gene_id
  } else if (from == "symbol") {
    names <- ref$hgnc_symbol
  } else {
    names <- ref$entrezgene_id
  }
  if (to == "ensembl") {
    vals <- ref$ensembl_gene_id
  } else if (to == "symbol") {
    vals <- ref$hgnc_symbol
  } else {
    vals <- ref$entrezgene_id
  }
  lookup <- setNames(vals, names)
  mapped <- lookup[vec]
  ifelse(is.na(mapped), vec, mapped)
}

random_palette <- function(min_length = NULL, continuous = FALSE) {
  if (continuous) {
    choices <- paletteer::palettes_c_names
  } else if (!is.null(min_length)) {
    choices <- paletteer::palettes_d_names |>
      dplyr::filter(length >= min_length)
  } else {
    choices <- paletteer::palettes_d_names
  }
  choices |>
    dplyr::slice_sample(n = 1) |>
    dplyr::select(package, palette) |>
    paste0(collapse = "::")
}

random_palette_d <- function(min_length = NULL) {
  random_palette(min_length = min_length, continuous = FALSE)
}

random_palette_c <- function() {
  random_palette(continuous = TRUE)
}

palette_from_cache <- function(key, min_length = NULL, discrete = TRUE) {
  if (discrete) {
    fn <- \() random_palette_d(min_length = min_length)
  } else {
    fn <- random_palette_c
  }
  if (!is.null(min_length)) {
    key <- glue::glue("palette_{key}_{min_length}")
  } else {
    key <- glue::glue("palette_{key}")
  }
  if (exists("CACHE")) {
    lookup <- CACHE$get(key)
    if (fastmap::is.key_missing(lookup)) {
      palette <- fn()
      CACHE$set(key, palette)
      palette
    } else {
      lookup
    }
  } else {
    fn()
  }
}

#' Retrieve the most recently saved resource `rname` from the
#' global cache
#'
from_bfc <- function(rname) {
  BiocFileCache::bfcquery(BFC, rname) |>
    dplyr::arrange(dplyr::desc(create_time)) |>
    head(n = 1) |>
    purrr::pluck("rpath") |>
    readRDS()
}
