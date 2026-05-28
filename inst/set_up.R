# Gather data resources and store in cache

devtools::load_all()

main <- function(cfg) {
  box::use(BiocFileCache[BiocFileCache, bfcnew, bfcquery], glue[glue])
  bfc <- BiocFileCache(cfg$cache)
  specs <- get_cache_spec(cfg)
  for (spec in specs) {
    if (nrow(bfcquery(bfc, spec$name)) > 0) {
      message(glue("Object `{spec$name}` found in cache already, skipping..."))
    } else {
      obj <- spec$fn()
      file_path <- bfcnew(bfc, spec$name, ext = ".rds")
      saveRDS(obj, file_path)
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
  }
  main(cfg)
}
