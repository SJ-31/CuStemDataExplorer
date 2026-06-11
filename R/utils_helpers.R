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
    if (v) "\u274c" else "\u2714\ufe0f"
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
  library(ggthemes)
  library(ggsci)
  if (continuous) {
    choices <- paletteer::palettes_c_names
  } else if (!is.null(min_length)) {
    choices <- paletteer::palettes_d_names |>
      dplyr::filter(length >= min_length)
    if (nrow(choices) == 0) {
      custom <- c()
      while (length(custom) < min_length) {
        pal_group <- sample(names(paletteer::palettes_d), size = 1)
        colors <- paletteer::palettes_d[[pal_group]] |>
          sample(size = 1) |>
          unlist(use.names = FALSE)
        custom <- c(custom, colors) |> unique()
      }
      return(custom)
    }
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

#' Retrieve a palette name from the in-memory cache for plotting consistency
#'
#' @description If the key for the palette isn't found in the cache, a new one
#' is generated and added to the cache. If the cache doesn't exist, simply
#' return a random palette
palette_from_cache <- function(
  key,
  key_group,
  min_length = NULL,
  discrete = TRUE
) {
  if (discrete) {
    fn <- \() random_palette_d(min_length = min_length)
  } else {
    fn <- random_palette_c
  }
  if (!is.null(min_length)) {
    key <- glue::glue("{key}_{min_length}")
  }

  if (exists("CACHE")) {
    if (fastmap::is.key_missing(CACHE$get(key_group))) {
      CACHE$set(key_group, list())
    }
    lookup <- CACHE$get(key_group)
    if (key %notin% names(lookup)) {
      palette <- fn()
      lookup[[key]] <- palette
      CACHE$set(key_group, lookup)
      palette
    } else {
      lookup[[key]]
    }
  } else {
    fn()
  }
}


add_palette_from_cache <- function(
  plot,
  key,
  key_group,
  min_length = NULL,
  fill = FALSE,
  discrete = TRUE
) {
  pal <- palette_from_cache(
    key,
    key_group = key_group,
    min_length = min_length,
    discrete = discrete
  )
  if (length(pal) > 1 && fill) {
    plot + ggplot2::scale_fill_discrete(pal)
  } else if (length(pal) > 1) {
    plot + ggplot2::scale_color_discrete(pal)
  } else if (discrete && fill) {
    plot + paletteer::scale_fill_paletteer_d(pal)
  } else if (discrete) {
    plot + paletteer::scale_color_paletteer_d(pal)
  } else if (fill) {
    plot + paletteer::scale_fill_paletteer_c(pal)
  } else {
    plot + paletteer::scale_color_paletteer_c(pal)
  }
}

#' Retrieve the most recently saved resource `rname` from the
#' global cache
#'
from_bfc <- function(q, as_row = FALSE, bfc = BFC, read = TRUE) {
  query_res <- BiocFileCache::bfcquery(bfc, q) |>
    dplyr::filter(rname == q)
  if (nrow(query_res) > 0) {
    row <- query_res |>
      dplyr::arrange(dplyr::desc(create_time)) |>
      head(n = 1)
    if (as_row) {
      row
    } else if (read && stringr::str_ends(row$rpath, ".rds")) {
      readRDS(row$rpath)
    } else if (read && stringr::str_ends(row$rpath, ".db")) {
      duckdb::dbConnect(duckdb::duckdb(), dbdir = row$rpath)
    } else {
      row$rpath
    }
  }
}

set_logger <- function() {
  file <- get_golem_config("log")
  logger::log_appender(logger::appender_tee(file))
  logger::log_info("Logger set up")
}

dispatch_normalize <- function(obj) {
  cfg <- local({
    tmp <- get_golem_config("expression_viewer") %||% list()
    tmp$normalization
  })
  method <- cfg$method %||% "vst"
  if (cfg$method == "tmm") {
    tmm_normalize(obj)
  } else {
    kws <- cfg$kws %||% list()
    deseq2_normalize(obj, method = method, kws = kws)
  }
}


get_validate_cache <- function() {
  cache_path <- get_golem_config("cache")
  bfc <- withCallingHandlers(
    expr = BiocFileCache::BiocFileCache(cache_path),
    error = \(e) {
      logger::log_fatal(
        "Failed to read from passed cache `{cache_path}`\n--- TRACEBACK ---\n {e}"
      )
      logger::log_fatal("Current configuration: {get_golem_config()}")
      cache_exists <- dir.exists(cache_path)
      logger::log_info("Cache existence: {cache_exists}")
      if (cache_exists) {
        logger::log_info("Cache contents: {list.files(cache_path)}")
      }
      NULL
    }
  )
  if (is.null(bfc)) {
    stop("Cache is null")
  }
  cached <- BiocFileCache::bfcinfo(bfc)$rname |> unique()
  vals <- vapply(
    get_cache_spec(NULL),
    \(spec) {
      key <- spec$name
      usage <- spec$usage
      if (key %notin% cached) {
        logger::log_warn(
          "The required key {key} is missing from the cache Run inst/set_up.R\n
--- Module(s)/feature(s) `{usage}` will not be available"
        )
      }
      key %in% cached
    },
    FUN.VALUE = logical(1)
  )
  if (all(!vals)) {
    stop("Cache contains none of the required keys. Run inst/set_up.R")
  }
  bfc
}
