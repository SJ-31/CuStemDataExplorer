cache_pca <- function(key, cache, cache_key, remove = TRUE) {
  list(
    name = paste0(key, "_pca"),
    fn = \() {
      data <- cache$get(cache_key)$expr
      as_mat <- tibble::column_to_rownames(data, var = "gene_id") |>
        as.matrix() |>
        t()
      pca <- prcomp(as_mat)
      if (remove) {
        cache$remove(cache_key)
      }
      pca
    }
  )
}

#' List of objects to cache and how
#'
#' @description
#' Each element of this list must be a list with three keys:
#' name: the name of the key used to cache the object
#' fn: a function of no arguments that returns the object to cache
#' usage: string listing the modules, features the object is required by.
#'      for logging and debugging purposes when the cache is empty
get_cache_spec <- function(cfg = NULL, cache, bfc) {
  list(
    list(
      name = "sheets",
      fn = \() get_sheets(cfg),
      usage = "Clinical, metadata, and sample tabs"
    ),
    list(
      name = "raw_bulk_expression",
      fn = \() {
        read_all_expr(
          cfg$expression_viewer,
          allowed_exts = c("csv", "tsv")
        )
      },
      usage = "DE analysis",
      cache = "raw_bulk_expr"
    ),
    list(
      name = "bulk_de",
      fn = \() {
        cache_de_analysis(
          cfg$bulk_de %||% list(),
          cache = cache,
          bfc = bfc,
          cache_key = "raw_bulk_expr",
          prefix = "bulk_de"
        )
      },
      usage = "DE viewer"
    ),
    list(
      name = "bulk_expression",
      fn = \() {
        norm <- dispatch_normalize(cache$get("raw_bulk_expr"))
        cache$remove("raw_bulk_expr")
        norm
      },
      usage = "Bulk expression comparison",
      cache_key = "bulk_expr"
    ),
    cache_pca("bulk_expression", cache = cache, cache_key = "bulk_expr"),
    list(
      name = "sc_pseudobulk_expression_raw",
      fn = \() read_all_expr(cfg$expression_viewer, allowed_exts = "h5ad"),
      cache = "single_cell_expr",
      usage = "Single-cell DE analysis"
    ),
    list(
      name = "sc_pseudobulk_expression",
      fn = \() {
        norm <- cache$get("single_cell_expr") |> dispatch_normalize()
        cache$remove("single_cell_expr")
        norm
      },
      usage = "Single-cell pseudobulk expression comparison",
      cache_key = "sc_psbulk_expr"
    ),
    cache_pca(
      "sc_pseudobulk_expression",
      cache = cache,
      cache_key = "sc_psbulk_expr"
    ),
    cache_exome(cfg, cache = cache, bfc = bfc, cfg_key = "exome"),
    cache_exome(cfg, cache = cache, bfc = bfc, cfg_key = "exome_sv")
  )
}

cache_exome <- function(cfg, cache, bfc, cfg_key) {
  list(
    name = cfg_key,
    fn = \() {
      data <- read_variant_spec_all(cfg[[cfg_key]])

      sapply(
        names(data),
        \(name) {
          key <- paste0(cfg_key, "_", name)
          db_path <- BiocFileCache::bfcnew(bfc, key, ext = ".db")
          cur <- data[[name]]
          con <- duckdb::dbConnect(duckdb::duckdb(), dbdir = db_path)
          for (ttype in names(cur)) {
            duckdb::dbWriteTable(con, ttype, cur[[ttype]])
          }
          db_path
        },
        simplify = FALSE,
        USE.NAMES = TRUE
      )
    }
  )
}

cache_de_analysis <- function(cfg, cache, bfc, cache_key, prefix) {
  result <- do_de(
    obj = cache$get(cache_key),
    how = cfg$how %||% "deseq2",
    kws = cfg$kws %||% list()
  )

  dds_name <- paste0(prefix, "::dds")
  dds_path <- BiocFileCache::bfcnew(bfc, dds_name, ext = ".rds")
  saveRDS(result$dds, dds_path)
  mapping <- list(dds = dds_name)

  cache_into <- function(subset) {
    name <- paste0(prefix, "::", subset)
    db_file <- BiocFileCache::bfcnew(bfc, name, ext = ".db")
    con <- duckdb::dbConnect(duckdb::duckdb(), dbdir = db_file)
    for (n in names(result[[subset]])) {
      df <- result[[subset]][[n]] |>
        as.data.frame() |>
        tibble::rownames_to_column(var = "gene")
      duckdb::dbWriteTable(con, n, df)
    }
    mapping[[subset]] <<- name
  }
  cache_into("ovr")
  cache_into("pairwise")

  mapping
}
