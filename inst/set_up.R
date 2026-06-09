# Gather data resources and store in cache

devtools::load_all()

main <- function(cfg, force_update = NULL) {
  box::use(BiocFileCache[BiocFileCache, bfcnew, bfcquery], glue[glue])
  bfc <- BiocFileCache(cfg$cache)
  cache <- cachem::cache_mem()
  specs <- get_cache_spec(cfg, cache, bfc)
  for (spec in specs) {
    if (nrow(bfcquery(bfc, spec$name, exact = TRUE)) > 0) {
      message(glue("Object `{spec$name}` found in cache already, skipping..."))
      obj <- from_bfc(q = spec$name, bfc = bfc)
    } else {
      obj <- spec$fn()
      file_path <- bfcnew(bfc, spec$name, ext = ".rds")
      saveRDS(obj, file_path)
    }
    if (!is.null(spec$cache)) {
      cache$set(spec$cache, obj)
    }
  }
}

if (sys.nframe() == 0) {
  library(optparse)
  parser <- OptionParser()
  parser <- add_option(
    parser,
    c("-f", "--file"),
    type = "character",
    help = "Configuration file",
    default = app_sys("golem-config.yml")
  )
  parser <- add_option(
    parser,
    c("-c", "--config"),
    type = "character",
    help = "Additional configuration to read from config",
    default = "default"
  )
  parser <- add_option(
    parser,
    c("-u", "--force_update"),
    action = "append",
    type = "character",
    help = "Force update the following items in the cache",
    default = NULL
  )
  parser <- add_option(
    parser,
    c("-r", "--remove"),
    action = "store_true",
    type = "logical",
    help = "Remove old cache and create a new one",
    default = FALSE
  )
  args <- parse_args(parser)
  cfg <- config::get(file = args$file, config = args$config)
  if (args$remove) {
    BiocFileCache::removebfc(
      BiocFileCache::BiocFileCache(cfg$cache),
      ask = FALSE
    )
  } else if (length(args$force_update) > 0) {
    bfc <- BiocFileCache::BiocFileCache(cfg$cache)
    for (to_remove in args$force_update) {
      q <- BiocFileCache::bfcquery(bfc, to_remove, exact = TRUE)
      if (nrow(q) > 0) {
        BiocFileCache::bfcremove(bfc, q$rid)
        sprintf("Removing key `%s`", q)
      } else {
        warning(sprintf(
          "Cannot force-update `%s`, doesn't exist in cache",
          q
        ))
      }
    }
  }

  main(cfg)
}
